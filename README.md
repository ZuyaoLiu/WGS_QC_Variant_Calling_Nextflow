# WGS QC and Variant Calling Pipeline

![Pipeline workflow](./Pipeline_figure.png)

A modular **Nextflow DSL2** pipeline for short-read, non-human whole-genome variant calling. The repository is organized around:

- `main.nf` for parameter handling and top-level orchestration
- `subworkflows/` for step-level pipeline composition
- `modules/` for single-responsibility process wrappers
- `vcf_filtering_scripts/` for helper Python and R scripts

## Overview

Main execution flow:

`Raw_QC -> Trimming_QC -> Aligning -> Calling -> optional VCF_Filtering`

Supported entry modes:

- `--run_step Raw_QC`
- `--run_step Trimming_QC`
- `--run_step Aligning`
- `--run_step Calling`
- `--run_step VCF_Filtering`
- `--run_step all`

Key design choices:

- `Raw_QC` is an inspection stop and exits after QC
- alignment emits **reference-embedded CRAM** instead of BAM
- calling is parallelized by `--interval_size`
- optional `vcffilter` currently implements **phase2** only

## Repository Layout

```text
.
├── main.nf
├── nextflow.config
├── vcffilter.config
├── subworkflows/
│   ├── raw_qc.nf
│   ├── trimming_qc.nf
│   ├── alignment.nf
│   ├── variant_calling.nf
│   └── vcf_filtering.nf
├── modules/
│   ├── download/
│   ├── qc/
│   ├── alignment/
│   ├── calling/
│   └── vcf_filtering/
├── vcf_filtering_scripts/
├── test_run/
└── pipeline_info/
```

Structure rules:

- `main.nf` should not contain tool-specific process logic
- `subworkflows/` should express step boundaries and dataflow
- `modules/` should wrap one concrete task or tightly related task group
- helper scripts should stay outside `modules/` unless they are tiny inline shell fragments

## Workflow Modules

### Raw QC

- local FASTQ input or SRA download
- raw `FastQC`
- raw `MultiQC`

### Trimming QC

- `fastp`
- post-trim `FastQC`
- post-trim `MultiQC`

### Alignment

- `bwa-mem2`
- `samtools collate` and `fixmate` for PE
- `sambamba markdup`
- final output: `SAMPLE.markdup.cram` and `SAMPLE.markdup.cram.crai`

### Variant Calling

- `bcftools` cohort calling
- `GATK` HaplotypeCaller / GenomicsDB / GenotypeGVCFs branch
- optional BQSR for GATK
- optional cleanup of per-interval gVCFs after GenomicsDBImport

### VCF Filtering

- controlled by `vcffilter.config`
- current implementation: `phase2`
- current scope: **biallelic SNP only**
- not supported: `indel` or mixed `SNP + indel` filtering

## Inputs

### Global

- `--run_step`: `Raw_QC | Trimming_QC | Aligning | Calling | VCF_Filtering | all`
- `--read_type`: `SE | PE`
- `--ref`: required for `Raw_QC`, `Trimming_QC`, `Aligning`, `Calling`, and `all`
- `--ref` is not required for standalone `VCF_Filtering`

### Step1 Raw_QC

Inputs:

- `--input_dir`: raw FASTQ directory
- `--SRA_list`: two-column SRA list file

Requirement:

- for `--run_step Raw_QC` or `--run_step all`, provide exactly one of `--input_dir` or `--SRA_list`

### Step2 Trimming_QC

Inputs:

- `--input_dir`: raw FASTQ directory

Requirement:

- for `--run_step Trimming_QC`, `--input_dir` is required

### Step3 Aligning

Inputs:

- `--trimmed_input_dir`: trimmed FASTQ directory

Behavior:

- if `--trimmed_input_dir` is not provided, the workflow reuses Step2 outputs from `results/02_fastp_qc/fastp`
- if Step2 outputs are missing and `--run_step Aligning` is used, the workflow auto-runs `Trimming_QC` once

### Step4 Calling

Inputs:

- Step3 CRAM outputs under `results/03_markdup`
- `--trimmed_input_dir`: optional one-step fallback input to `Aligning`
- `--bqsr_panel`: known-sites VCF/VCF.GZ when `--caller gatk --use_bqsr true`

Requirement:

- for any path that executes `Calling`, `--interval_size` is required

### Step5 VCF_Filtering

Inputs:

- Step4 cohort VCF when `--enable_vcffilter true`
- `--input_vcf`: standalone VCF input for `--run_step VCF_Filtering`

Requirement:

- for `--run_step VCF_Filtering`, `--input_vcf` is required
- standalone mode also expects a tabix index at `INPUT_VCF.tbi`

Input naming:

- PE: `sample_R1.fastq.gz` and `sample_R2.fastq.gz`
- SE: `sample.fastq.gz` or `sample.fq.gz`

SRA list format:

- column 1: SRA accession
- column 2: preferred sample name

Reference constraint:

- reference headers must not use `chr:start-end` style names, because they conflict with region parsing in downstream tools

## Running the Pipeline

Check help:

```bash
cd test_run/run
nextflow run ../../main.nf --help
nextflow run ../../main.nf --help_full
```

Run the full pipeline:

```bash
cd test_run/run
nextflow run ../../main.nf \
  -profile local \
  --input_dir ../data \
  --ref ../data/sim_ref_100kb.fa \
  --read_type PE \
  --run_step all \
  --caller bcftools \
  --interval_size 50K
```

Start from alignment with external trimmed reads:

```bash
cd test_run/run
nextflow run ../../main.nf \
  -profile local \
  --trimmed_input_dir ../data \
  --ref ../data/sim_ref_100kb.fa \
  --read_type PE \
  --run_step Aligning \
  --caller bcftools \
  --interval_size 50K
```

Run calling plus vcffilter phase2:

```bash
cd test_run/run
nextflow run ../../main.nf \
  -profile local \
  --trimmed_input_dir ../data \
  --ref ../data/sim_ref_100kb.fa \
  --read_type PE \
  --run_step Calling \
  --caller bcftools \
  --interval_size 50K \
  --enable_vcffilter true
```

Run standalone vcffilter phase2 from an external VCF:

```bash
cd test_run/run
nextflow run ../../main.nf \
  -profile local \
  --run_step VCF_Filtering \
  --input_vcf /path/to/cohort.vcf.gz
```

## Execution Behavior

- `Raw_QC` runs Step 1 only and exits
- `Trimming_QC` runs from trimming to the current pipeline end
- `Aligning` can consume `--trimmed_input_dir` directly
- `Calling` can reuse published Step 3 CRAM outputs
- `VCF_Filtering` can run either from Step4 output or directly from `--input_vcf`
- only one fallback level is allowed when resuming from later steps

## Configuration

### Main pipeline config

Use [nextflow.config](/scratch/zuyao_20260402/WGS_QC_Variant_Calling_Nextflow/nextflow.config) for:

- executor profiles
- process resources
- global defaults
- container settings

Default container image:

- `leonardliu0910/wgs_variant_calling:latest`

Default container runtime:

- `Apptainer`

Container selection rules:

- use `--image` to select a Docker/registry image
- use `--sif` to select a local Singularity or Apptainer image file
- if `--sif` is provided, it overrides `--image` for Singularity/Apptainer execution

### VCF filter config

Use [vcffilter.config](/scratch/zuyao_20260402/WGS_QC_Variant_Calling_Nextflow/vcffilter.config) for:

- `phase1` parameters
- `phase2` parameters
- default filtering thresholds

The file uses short parameter names grouped by comments instead of nested Nextflow maps.

## Outputs

```text
results/
├── 00_sra_download/
├── 01_raw_qc/
├── 02_fastp_qc/
├── 03_markdup/
│   ├── SAMPLE.markdup.cram
│   └── SAMPLE.markdup.cram.crai
├── 04_variant_calling/
│   ├── bcftools/
│   └── gatk/
└── 05_vcf_filtering/
    ├── 02_phase2/
    └── final/
```

Representative final outputs:

- `cohort.bcftools.vcf.gz`
- `cohort.gatk.filtered.vcf.gz`
- `cohort.vcffilter.phase2.vcf.gz`

## Profiles

Supported profiles:

- `local`
- `slurm`
- `awsbatch`

Resource tuning should be done in `nextflow.config`, not by editing module scripts.

Container engines:

- Apptainer is enabled by default
- Docker is supported through `-with-docker`
- Singularity is supported through `-with-singularity`
- Apptainer is supported through `-with-apptainer`
- use `--image` with any engine when pulling from Docker Hub or another registry
- use `--sif` only with `-with-singularity` or `-with-apptainer`
- for local conda-based testing without containers, pass `-c /tmp/wgs_qc_no_container.config`

## Current Limitations

- `vcffilter` currently implements `phase1` review and `phase2` final filtering
- `vcffilter_stop_after phase1` stops after the phase1 review outputs
- `vcffilter` supports **biallelic SNP** only
- the deprecated file `VCF_filtering.config` is kept only as a compatibility marker

## Test Data

The repository includes a small runnable dataset under [test_run](/scratch/zuyao_20260402/WGS_QC_Variant_Calling_Nextflow/test_run).

Use it for smoke tests after structural or module changes.

## BQSR Panel Requirements

If you enable `--caller gatk --use_bqsr true`, you must provide an external known-sites panel with `--bqsr_panel`.

Recommended panel characteristics:

- High-confidence variant sites only.
- Indexed VCF/VCF.GZ built on the same reference FASTA used by the pipeline.
- Stable sites supported by reliable genotype and depth patterns.

For species without a standard community BQSR panel, a practical strategy is:

1. Run the `bcftools` branch first and generate a raw cohort VCF.
2. Apply very strict filtering to retain only highly reliable sites.
3. Use the resulting filtered VCF as `--bqsr_panel` in a second GATK+BQSR run.

Example strict criteria for constructing a candidate BQSR panel  (for reference only, please set you own cut-off instead):

- `QUAL >= 999`
  Site-level variant quality. Higher values indicate stronger evidence that the site is a true variant.
  
- `AF >= 0.1`
  Allele frequency. This avoids very low-frequency sites that may be unstable or weakly supported.
  
- `GQ >= 50`
  Genotype quality. This keeps only highly confident per-sample genotype calls.
  
- Per-sample `DP` within a sample-specific acceptable range
  Depth outside the expected range for a given sample can indicate under-covered or problematic genotypes.
  A practical approach is: compute each sample's mean DP using non-missing genotypes, then mask calls outside an allowed interval such as `0.5x` to `2.5x` that sample mean.
  The in-house script can be found in `bqsr_filter_script/Filter_DP_per_Sample.py`.
  
  An alternative approche is: filter based on meanDP across all samples and sites (recommend for evenly sequenced projects).
  
- Allelic-balance check passed
  This is mainly for heterozygous genotypes. Strong imbalance between REF and ALT support can indicate false heterozygous calls.
  

## FAQ

### Why might I still filter sites by missing fraction when building a BQSR panel?

BQSR itself is applied per sample on each alignment file, so site missingness is not a direct requirement of the BQSR algorithm. However, if you are building your own known-sites panel, missing fraction can still help remove unstable cohort sites.

- A site with many missing genotypes may reflect poor mappability, uneven depth, local complexity, or unreliable genotype calling.
- Such sites are often a poor choice for a high-confidence known-sites panel, even if BQSR does not explicitly require cohort-wide completeness.
- A threshold such as `missing proportion <= 0.1` is therefore best treated as a recommended panel-construction filter, not as a hard requirement of the BQSR algorithm itself.

## Software Scope (Intentional Constraints)

Included only:

- FastQC
- MultiQC
- fastp
- bwa
- bwa-mem2
- samtools
- bcftools
- GATK

For Script `bqsr_filter_script/Filter_DP_per_Sample.py`, PySam needs to be installed.

## Reproducibility Notes

- Every process declares `container`.
- Workflow report/trace/timeline/DAG are written under `pipeline_info/`.
- Intermediate task files are under `work/`.
