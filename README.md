# WGS QC and Variant Calling Pipeline

**Version:** 3.1.0
**Workflow engine:** Nextflow DSL2
**Scope:** short-read whole-genome QC, alignment, joint variant calling, and optional SNP filtering. Supports any reference (tested on human and non-human genomes).

![Pipeline workflow](./Pipeline_figure.png)

---

## Table of contents

1. [Overview](#overview)
2. [Quickstart](#quickstart)
3. [Input formats](#input-formats)
4. [Run modes](#run-modes)
5. [Parameters](#parameters)
6. [Execution profiles](#execution-profiles)
7. [Performance and cluster tuning](#performance-and-cluster-tuning)
8. [Outputs](#outputs)
9. [Step behavior](#step-behavior)
10. [Testing](#testing)
11. [Notes and limitations](#notes-and-limitations)
12. [Repository layout](#repository-layout)

---

## Overview

The pipeline runs the following stages, each one optional via `--run_step`:

```text
Raw_QC -> Trimming_QC -> Aligning -> Calling -> (optional) VCF_Filtering
```

Main features:

- Raw FASTQ QC with FastQC and MultiQC
- Adapter/quality trimming with fastp
- Alignment with bwa-mem2 and duplicate marking with sambamba
- Reference-embedded CRAM output for storage efficiency
- Cohort variant calling with `bcftools` (joint mpileup+call) or `GATK` (HaplotypeCaller -> GenomicsDBImport -> GenotypeGVCFs joint genotyping)
- Optional BQSR for the GATK path
- Optional post-calling VCF filtering for biallelic SNPs (phase1 review + phase2 final)
- Local, SLURM, and AWS Batch execution profiles
- Automatic retry on transient task failures; safe `-resume` across all steps

---

## Quickstart

### 1. Install dependencies

- [Nextflow](https://www.nextflow.io/) >= 25.10.4
- A container runtime: Apptainer/Singularity (default) or Docker

The default container image is `leonardliu0910/wgs_variant_calling:latest`. To use a local image, pass `--sif /path/to/image.sif`.

### 2. Pick an input mode

| Input mode | Flag | Use when |
|---|---|---|
| Raw FASTQ directory | `--input_dir` | You already have FASTQ files on disk |
| SRA accessions | `--SRA_list` | You want the pipeline to download from SRA |
| Trimmed FASTQ directory | `--trimmed_input_dir` | You already trimmed and want to skip Trimming_QC |
| Existing cohort VCF | `--input_vcf` | You only want post-calling VCF filtering |

### 3. Run the full pipeline

```bash
nextflow run main.nf \
  -profile slurm \
  --input_dir /path/to/fastq \
  --ref /path/to/reference.fa \
  --read_type PE \
  --caller bcftools \
  --interval_size 50M
```

This runs Raw_QC -> Trimming_QC -> Aligning -> Calling and writes results under `results/`.

For a quick CLI reference:

```bash
nextflow run main.nf --help        # concise help
nextflow run main.nf --help_full   # every parameter, including phase1/phase2 filter knobs and tuning
```

---

## Input formats

### Paired-end FASTQ (`--read_type PE`)

Filenames must end in `_R1` / `_R2`:

```text
sampleA_R1.fastq.gz
sampleA_R2.fastq.gz
sampleB_R1.fastq.gz
sampleB_R2.fastq.gz
```

### Single-end FASTQ (`--read_type SE`)

```text
sampleA.fastq.gz
sampleB.fq.gz
```

### SRA list (`--SRA_list`)

Tab-separated, at least two columns. Lines starting with `#` are skipped:

```text
# accession    sample_id
SRRxxxxxxx    sampleA
SRRyyyyyyy    sampleB
```

A documented template is provided at [`tests/SRA_list_human_test.tsv`](tests/SRA_list_human_test.tsv).

### Standalone VCF input

`--input_vcf` requires a tabix index next to the file:

```text
/path/to/cohort.vcf.gz
/path/to/cohort.vcf.gz.tbi
```

---

## Run modes

### Raw QC only

```bash
nextflow run main.nf -profile slurm \
  --input_dir /path/to/fastq --ref /path/to/reference.fa \
  --read_type PE --run_step Raw_QC
```

### Start from trimmed FASTQ

```bash
nextflow run main.nf -profile slurm \
  --trimmed_input_dir /path/to/trimmed_fastq \
  --ref /path/to/reference.fa --read_type PE \
  --run_step Aligning --caller bcftools --interval_size 50M
```

### GATK joint calling with BQSR

```bash
nextflow run main.nf -profile slurm \
  --input_dir /path/to/fastq --ref /path/to/reference.fa \
  --read_type PE --run_step all \
  --caller gatk --interval_size 50M \
  --use_bqsr true --bqsr_panel /path/to/known_sites.vcf.gz
```

### Standalone VCF filtering

```bash
nextflow run main.nf -profile slurm \
  --run_step VCF_Filtering --input_vcf /path/to/cohort.vcf.gz \
  --enable_vcffilter true --phase2_max_sample_missing 0.3
```

### Download from SRA

```bash
nextflow run main.nf -profile slurm \
  --SRA_list tests/SRA_list_human_test.tsv \
  --ref /path/to/reference.fa \
  --read_type PE --interval_size 50M
```

---

## Parameters

For the complete list (every advanced flag, every phase1/phase2 threshold), run:

```bash
nextflow run main.nf --help_full
```

### Required and common

| Parameter | Default | Allowed | Description |
|---|---|---|---|
| `--run_step` | `all` | `Raw_QC`,`Trimming_QC`,`Aligning`,`Calling`,`VCF_Filtering`,`all` | Starting step. Continues to end unless `Raw_QC` or `VCF_Filtering`. |
| `--input_dir` | `null` | directory | Raw FASTQ directory. |
| `--SRA_list` | `null` | TSV file | Two-column accession/sample list. Use instead of `--input_dir`. |
| `--trimmed_input_dir` | `null` | directory | Pre-trimmed FASTQ directory; entry point for `Aligning`. |
| `--input_vcf` | `null` | `.vcf.gz` | Standalone VCF filtering input. Requires `.tbi`. |
| `--ref` | `null` | FASTA | Reference. Required for everything except standalone `VCF_Filtering`. |
| `--read_type` | `PE` | `PE`,`SE` | Paired-end or single-end. |
| `--caller` | `bcftools` | `bcftools`,`gatk` | Variant caller backend. |
| `--interval_size` | `null` | e.g. `500K`,`50M`,`1G` | Calling interval block size. Required whenever `Calling` runs. |
| `--enable_vcffilter` | `false` | `true`,`false` | Run post-calling VCF filtering. |
| `--vcffilter_stop_after` | `final` | `phase1`,`final` | Stop after phase1 review or continue to phase2 final. |

### GATK and BQSR

| Parameter | Default | Description |
|---|---|---|
| `--use_bqsr` | `false` | Enable BaseRecalibrator + ApplyBQSR. Requires `--caller gatk`. |
| `--bqsr_panel` | `null` | Known-sites VCF/VCF.GZ, required when `--use_bqsr true`. |
| `--cleanup_GVCFS` | `true` | Delete per-interval gVCFs after GenomicsDBImport succeeds. |
| `--gatk_genomicsdbimport_batch_size` | `50` | GenomicsDBImport `--batch-size`. Raise to 200+ for very large cohorts. |
| `--gatk_genomicsdbimport_consolidate` | `false` | Add `--consolidate` to GenomicsDBImport. Helps downstream GenotypeGVCFs on large cohorts. |

### Output and step toggles

| Parameter | Default | Description |
|---|---|---|
| `--publish_cram` | `true` | Publish markdup CRAM under `results/03_markdup/`. Set `false` if you only need VCF and want to save inodes. |
| `--skip_raw_fastqc` | `false` | Skip Raw FastQC + Raw MultiQC. Only honored with `--run_step all`. |

### Container

| Parameter | Default | Description |
|---|---|---|
| `--image` | `leonardliu0910/wgs_variant_calling:latest` | Docker / registry image. |
| `--sif` | `null` | Local Singularity / Apptainer image. Overrides `--image` when set. |

### VCF filtering thresholds

Defaults are defined in `nextflow.config`. The full list is shown in `--help_full`. Highlights:

| Parameter | Default | Description |
|---|---|---|
| `--phase1_min_gq` | `20` | Mask genotypes with GQ below this value. |
| `--phase1_min_dp` | `5` | Mask genotypes with DP below this value. |
| `--phase1_enable_ab_mask` | `true` | Enable heterozygous allele-balance masking. |
| `--phase1_sample_dp_min_factor` / `--phase1_sample_dp_max_factor` | `0.5` / `2.5` | Per-sample DP factor cutoffs. |
| `--phase2_remove_high_missing_samples` | `true` | Drop samples failing the missingness threshold. |
| `--phase2_max_sample_missing` | `null` | Max sample missingness fraction. User must set when removing high-missing samples. |
| `--phase2_min_qual` | `30` | Minimum site QUAL. |
| `--phase2_max_site_missing` | `0.2` | Maximum site missingness fraction. |
| `--phase2_enable_site_mean_dp_filter` | `true` | Apply mean site-DP filter. |
| `--phase2_site_mean_dp_min_factor` / `--phase2_site_mean_dp_max_factor` | `0.5` / `2.5` | Site mean-DP factor cutoffs. |

### Advanced tool pass-through

These accept verbatim flag strings appended to the underlying tool invocation. Quote them.

| Parameter | Tool |
|---|---|
| `--fastp_parameters` | `fastp` |
| `--bwamem2_parameters` | `bwa-mem2 mem` |
| `--bcftools_mpileup_parameters` | `bcftools mpileup` |
| `--bcftools_call_parameters` | `bcftools call` |
| `--bcftools_concat_parameters` | `bcftools concat` |
| `--gatk_haplotypecaller_parameters` | `GATK HaplotypeCaller` |
| `--gatk_genomicsdbimport_parameters` | `GATK GenomicsDBImport` |
| `--gatk_genotypegvcfs_parameters` | `GATK GenotypeGVCFs` |
| `--gatk_baserecalibrator_parameters` | `GATK BaseRecalibrator` |
| `--gatk_applybqsr_parameters` | `GATK ApplyBQSR` |
| `--gatk_variantfiltration_snp_parameters` | `GATK VariantFiltration` (SNP) |
| `--gatk_variantfiltration_indel_parameters` | `GATK VariantFiltration` (indel) |

---

## Execution profiles

| Profile | Executor | Defaults |
|---|---|---|
| `local` | local | `cpus = 4` global default. Apptainer enabled. |
| `slurm` | SLURM | `cpus = 2`, `memory = 20.GB`, `queueSize = 50`, `submitRateLimit = '50/sec'`. Apptainer enabled. |
| `awsbatch` | AWS Batch | Template only. Edit queue, region, container, and S3 work directory before use. |

Switch with `-profile <name>`. Per-process resource overrides live in `nextflow.config` under each profile's `withName:` blocks (commented-out templates with recommended values — see next section).

---

## Performance and cluster tuning

### Default error handling

Tasks retry up to 3 times on transient failures (`task.attempt <= 3 ? 'retry' : 'finish'`). On the 4th failure, the pipeline switches to `finish` so other running branches still produce output. Override in a custom `-c` config if you need fail-fast.

### SLURM throughput

The SLURM profile keeps up to 50 sbatch submissions in flight (`executor.queueSize = 50`) and rate-limits new submissions to 50/sec. Raise `queueSize` (e.g. 200–500) for large cohorts when your cluster fair-share / `MaxJobs` allows.

### Per-process resources

`nextflow.config` ships commented-out `withName:` templates with recommended values for the heavy processes:

| Process | cpus | memory |
|---|---:|---|
| `FASTP_SE/PE_PROCESS` | 8 | 8 GB |
| `ALIGN_MARKDUP_SE/PE_PROCESS` | 16 | 32 GB |
| `GATK_BQSR_PROCESS` | 4 | 16 GB |
| `GATK_HAPLOTYPECALLER_BY_CHR_PROCESS` | 4 | 16 GB |
| `GATK_GENOMICSDBIMPORT_BY_CHR_PROCESS` | 2 | 32 GB |
| `GATK_GENOTYPEGVCFS_BY_CHR_PROCESS` | 2 | 16 GB |
| `BCFTOOLS_CALL_BY_CHR_PROCESS` | 8 | 16 GB |

Uncomment the blocks you want active. Everything else inherits the global profile defaults.

### GenomicsDBImport scaling

For cohorts > 100 samples consider:

```bash
--gatk_genomicsdbimport_batch_size 200 \
--gatk_genomicsdbimport_consolidate true
```

The process's `-Xmx` already scales with `task.memory - 4 GB` (floor 2 GB).

---

## Outputs

```text
results/
├── 00_sra_download/           # only when --SRA_list is used
├── 01_raw_qc/
│   ├── fastqc/
│   └── multiqc/
├── 02_fastp_qc/
│   ├── fastp/
│   ├── fastqc/
│   └── multiqc/
├── 03_markdup/                # markdup CRAM (skipped when --publish_cram false)
│   ├── SAMPLE.markdup.cram
│   └── SAMPLE.markdup.cram.crai
├── 04_variant_calling/
│   ├── bcftools/              # bcftools caller output
│   └── gatk/                  # gatk caller output (HC -> GDB -> Genotype -> filter)
└── 05_vcf_filtering/          # only when --enable_vcffilter true
    ├── 01_phase1/
    ├── 02_phase2/
    └── final/
```

Run reports:

```text
pipeline_info/
├── timeline.html
├── report.html
├── trace.txt
└── dag.html
```

Final cohort VCF locations:

| Caller | Path |
|---|---|
| bcftools | `results/04_variant_calling/bcftools/cohort.bcftools.vcf.gz` |
| gatk     | `results/04_variant_calling/gatk/filter/cohort.gatk.filtered.vcf.gz` |
| vcffilter | `results/05_vcf_filtering/final/*.vcffilter.phase2.vcf.gz` |

---

## Step behavior

| `--run_step` | Behavior |
|---|---|
| `Raw_QC` | Raw FastQC + Raw MultiQC only, then exit. |
| `Trimming_QC` | Trimming -> Aligning -> Calling (-> VCF_Filtering if enabled). Requires `--input_dir`. |
| `Aligning` | Starts from `--trimmed_input_dir`; otherwise reuses `results/02_fastp_qc/fastp`, falling back to Trimming_QC once if missing. |
| `Calling` | Starts from `results/03_markdup`; falls back to Aligning once if missing. No multi-level fallback. |
| `VCF_Filtering` | Phase1 / phase2 filtering only, from `--input_vcf`. |
| `all` | Full pipeline from raw FASTQ or SRA input. |

---

## Testing

For a smoke test against public data, edit the included template:

```text
tests/SRA_list_human_test.tsv
```

Replace the placeholder accessions with real SRR / ERR / DRR IDs and run:

```bash
nextflow run main.nf -profile slurm \
  --SRA_list tests/SRA_list_human_test.tsv \
  --ref /path/to/human_reference.fa \
  --read_type PE --interval_size 50M --run_step all
```

For a fast smoke test, pick small accessions (< 5 GB) and use a smaller `--interval_size` (e.g. `10M`).

---

## Notes and limitations

- Reference sequence names must not contain `:` or `-`, because those clash with the region parser used for interval generation.
- VCF filtering currently supports biallelic SNPs only — indel-only and mixed SNP+indel inputs are out of scope.
- For `--caller gatk --use_bqsr true`, the known-sites panel must be built against the same reference assembly.
- `-resume` is now fully supported; partial reruns will pick up where they left off as long as `--cleanup_GVCFS` did not delete required intermediates.
- See `nextflow run main.nf --help` for a concise CLI reference and `--help_full` for the complete parameter list.

---

## Repository layout

```text
.
├── main.nf                       # Workflow entry point
├── nextflow.config               # Params + profiles + resources (vcffilter params merged in)
├── subworkflows/                 # Step-level subworkflows
├── modules/
│   ├── alignment/
│   ├── calling/                  # Includes prepare_reference.nf for fai/dict staging
│   ├── download/
│   ├── qc/
│   └── vcf_filtering/
├── vcf_filtering_scripts/        # Python/R helpers used inside vcf_filtering processes
├── tests/
│   └── SRA_list_human_test.tsv   # SRA accession list template
├── Pipeline_figure.png
├── README.md

```
