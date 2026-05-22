# v3.1.0 — Correctness, performance, and configuration cleanup

This release fixes a class of correctness bugs around `-resume`,
removes a real race condition on shared reference index files, makes
the pipeline robust to transient cluster failures, and consolidates
configuration and help output. There are no changes to the public
input/output schema.

## Highlights

- `-resume` now actually works after every step. The combination of
  `publishDir mode: 'move'` + `file("results/...")` re-staging was
  invalidating the work-directory cache. All publish modes are now
  `symlink` (intermediates that are not user-facing simply skip
  `publishDir`), and downstream channels read directly from the
  upstream process output.
- Removed a `.fai` / `.dict` race condition. Six calling processes
  each guarded an inline `samtools faidx` / `gatk CreateSequenceDictionary`
  against the shared reference directory. Under concurrent task launches
  this could leave the index in a corrupted state. A new
  `PREPARE_REFERENCE` process generates indices once and emits a value
  channel that all callers consume as staged input.
- The pipeline no longer aborts on a single transient failure. Tasks
  retry up to 3 times; if a task still fails, the run switches to
  `finish` so other in-flight branches can produce output instead of
  being killed.
- GenomicsDBImport tuning is now correct. JVM `-Xmx` scales with
  `task.memory`, `--batch-size` defaults to 50 (GATK recommendation,
  was 1 which disabled batch IO), and the conflicting
  `-XX:ActiveProcessorCount=1` is gone so `--reader-threads` actually
  uses the allocated CPUs.

## What changed

### P0 — correctness
- **Resume + publish layout (#§2.1).** Replaced `mode: 'move'` with
  `mode: 'symlink'` across markdup, calling, QC, vcf_filtering and
  download modules. Dropped the `.map { ... file("results/...") }`
  rewrap in `emit:` blocks. Added `--publish_cram` (default `true`)
  for cohorts that want to skip publishing markdup CRAM.
- **Reference indexing race (#§2.2).** Added
  `modules/calling/prepare_reference.nf`. Six calling processes
  (`gatk_haplotypecaller`, `gatk_genomicsdbimport`, `gatk_bqsr`,
  `gatk_genotypegvcfs`, `gatk_filter`, `bcftools_joint_call`) and the
  `variant_calling` subworkflow updated to consume
  `tuple(ref_fa, ref_fai, ref_dict)`.
- **Error strategy (#§2.3).** `errorStrategy = 'terminate'` →
  `{ task.attempt <= 3 ? 'retry' : 'finish' }`, `maxRetries = 3`.

### P1 — performance
- **GenomicsDBImport (#§2.6).** `-Xmx` scales with `task.memory - 4g`
  (floor 2g). New params `gatk_genomicsdbimport_batch_size` (default
  50) and `gatk_genomicsdbimport_consolidate` (default false). Removed
  `-XX:ActiveProcessorCount=1`.
- **Thread environment cleanup (#§2.7).** Dropped redundant
  `OMP_NUM_THREADS` / `MKL_NUM_THREADS` / `OPENBLAS_NUM_THREADS` /
  `POLARS_MAX_THREADS` / `RAYON_NUM_THREADS` exports from 9 processes.
  None of GATK (JVM), bcftools, samtools, or MultiQC read these
  variables — they used their own `--threads` / `-@` flags, so the
  env exports only risked thread oversubscription.
- **SLURM throughput (#§2.8).** `executor.queueSize` raised from 4 to
  50. Added `executor.submitRateLimit = '50/sec'`.

### P2 — defaults and small wins
- `cleanup_GVCFS` default flipped to `true`. Per-interval gVCFs are
  redundant once GenomicsDBImport finishes (#§2.9).
- `bcftools concat -a` removed. The interval list is non-overlapping
  by construction, so `-a` only added a sort-dedup pass (#§2.10).
- `fasterq-dump` postprocess switched from `gzip` to
  `pigz -p task.cpus` (#§2.11).
- New `--skip_raw_fastqc` toggle to skip Raw FastQC + Raw MultiQC when
  running `--run_step all` (#§2.12).
- `withName` per-process resource templates in the SLURM profile now
  carry the recommended cpus/memory values (FASTP 8/8, ALIGN_MARKDUP
  16/32, GATK_BQSR 4/16, HaplotypeCaller 4/16, GenomicsDBImport 2/32,
  GenotypeGVCFs 2/16, BCFTOOLS_CALL 8/16). Blocks remain commented out
  by default; users uncomment what they need (#§2.13).

### P3 — configuration cleanup
- `vcffilter.config` merged into `nextflow.config`. Fixes a latent
  duplicate-definition bug where `enable_vcffilter` and
  `vcffilter_stop_after` were defined in both files and the include
  order silently picked one.

### Documentation
- `--help` rewritten as a concise 45-line quick reference focused on
  required parameters and common options.
- `--help_full` lists every parameter the pipeline understands,
  including the full phase1/phase2 VCF-filtering tables and the
  cluster-tuning knobs.
- Added `tests/SRA_list_human_test.tsv` — empty SRA accession list
  template with usage notes for human-WGS smoke tests.

## New parameters

| Param | Default | Purpose |
| --- | --- | --- |
| `publish_cram` | `true` | Publish markdup CRAM to `results/03_markdup/`. Set to `false` on large cohorts that only need VCF output. |
| `skip_raw_fastqc` | `false` | Skip Raw FastQC + Raw MultiQC when running `--run_step all`. |
| `gatk_genomicsdbimport_batch_size` | `50` | Pass-through to GenomicsDBImport `--batch-size`. |
| `gatk_genomicsdbimport_consolidate` | `false` | Adds `--consolidate` to GenomicsDBImport. Recommended for large cohorts. |

## Changed defaults (potentially user-visible)

- `cleanup_GVCFS`: `false` → `true`. To keep per-interval gVCFs, pass
  `--cleanup_GVCFS false`.
- `executor.queueSize` (SLURM): `4` → `50`. If your cluster has a
  stricter submit cap, override in a custom `-c` config.
- `process.errorStrategy`: `terminate` → up to 3 retries then
  `finish`. Override to `terminate` if you specifically want
  fail-fast behavior.

## Removed

- `vcffilter.config`. All its parameters are now in `nextflow.config`.
  If you previously copied `vcffilter.config` next to a run, delete
  that copy — it will no longer be picked up.

## Upgrade notes

- No changes to required CLI flags. Existing `nextflow run` commands
  continue to work.
- If you previously relied on `cleanup_GVCFS=false` (e.g. for
  debugging joint-genotyping issues), explicitly set
  `--cleanup_GVCFS false`.
- If you previously bypassed `-resume` because it never re-used cache
  anyway, you should now see real cache hits on reruns.
- If you depended on `terminate` for fail-fast, set
  `process.errorStrategy = 'terminate'` and `maxRetries = 0` in a
  custom config.

## Acknowledgements

Design analysis and prioritization captured in
`federated-rolling-eagle.md` (included in this release).
