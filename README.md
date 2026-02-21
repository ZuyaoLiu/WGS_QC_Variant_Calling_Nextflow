# NGS Variant Calling Nextflow Pipeline

A production-style, containerized, modular **Nextflow DSL2** pipeline for **non-human WGS variant calling**.

This project is designed for strict step-wise execution and reproducibility under:

- **Nextflow**: `25.10.4`
- **Execution runtime**: Apptainer/Singularity
- **Container image**: `WGS_Variant_Calling.sif`

The pipeline follows the required business flow:

- `Raw_QC -> Trimming_QC -> Aligning -> Calling`

and supports both full-run and single-step modes.

## Key Capabilities

- Strict modular structure by step (`modules/step1_*` ... `modules/step4_*`)
- Explicit container declaration in every process
- SE/PE-aware alignment branches
- Caller-aware variant calling branches (`bcftools` / `gatk`)
- Optional GATK BQSR path
- **Per-chromosome/scaffold parallel calling** in Step4, followed by controlled merge
- Standardized output hierarchy under `results/`

## Repository Structure

```text
.
├── main.nf
├── nextflow.config
├── modules/
│   ├── step1_raw_qc/
│   │   ├── fastqc_raw.nf
│   │   └── multiqc_raw.nf
│   ├── step2_fastp_qc/
│   │   ├── fastp.nf
│   │   ├── fastqc_post.nf
│   │   └── multiqc_post.nf
│   ├── step3_alignment/
│   │   ├── bwa_index.nf
│   │   ├── bwa_se.nf
│   │   ├── bwa_mem2_index.nf
│   │   ├── bwa_mem2_pe.nf
│   │   ├── fixmate.nf
│   │   ├── sort.nf
│   │   └── sambamba_markdup.nf
│   └── step4_variant_calling/
│       ├── bcftools_call.nf
│       ├── gatk_bqsr.nf
│       ├── gatk_haplotypecaller.nf
│       ├── gatk_genomicsdbimport.nf
│       ├── gatk_genotypegvcfs.nf
│       ├── gatk_merge_raw.nf
│       └── gatk_filter.nf
└── tests/
    ├── run_smoke_test.sh
    └── README.md
```

## Software Scope (Intentional Constraints)

Included only:

- FastQC
- MultiQC
- fastp
- bwa
- bwa-mem2
- samtools (index/sort/fixmate/faidx helpers)
- sambamba
- bcftools
- GATK


## Execution Modes

### 1) Full run

```bash
nextflow run main.nf --run_step all ...
```

Order:

- Step1 `Raw_QC`
- Step2 `Trimming_QC`
- Step3 `Aligning`
- Step4 `Calling`

### 2) Single-step run

```bash
nextflow run main.nf --run_step Raw_QC ...
nextflow run main.nf --run_step Trimming_QC ...
nextflow run main.nf --run_step Aligning ...
nextflow run main.nf --run_step Calling ...
```

Single-step mode executes only the minimum required upstream logic embedded in that step path.

## Workflow Design Details

## Step1: Raw_QC

- FastQC on raw FASTQ
- MultiQC aggregation

Outputs:

- `results/01_raw_qc/fastqc/`
- `results/01_raw_qc/multiqc/`

## Step2: Trimming_QC

- fastp trimming/QC
- FastQC on post-fastp FASTQ
- MultiQC aggregation

Outputs:

- `results/02_fastp_qc/fastp/`
- `results/02_fastp_qc/fastqc/`
- `results/02_fastp_qc/multiqc/`

## Step3: Aligning

Branching by `--read_type`:

- `SE`:
  - `bwa index`
  - `bwa mem`
  - `samtools sort`
- `PE`:
  - `bwa-mem2 index`
  - `bwa-mem2 mem`
  - `samtools fixmate`
  - `samtools sort`

Common merge point:

- `sambamba markdup`

Outputs:

- `results/03_align/`
- `results/04_markdup/`

## Step4: Calling

Branching by `--caller`:

### A) `bcftools`

Per chromosome/scaffold parallelism:

- build intervals from reference FASTA headers (`>chr`, `>scaffold`, etc.)
- run per-interval bcftools calling in parallel
- merge per-interval VCFs into one final sample VCF

Outputs:

- per-interval: `results/05_variant_calling/bcftools/per_chrom/`
- merged: `results/05_variant_calling/bcftools/SAMPLE.bcftools.vcf.gz`

### B) `gatk`

#### `--use_bqsr false`

- HaplotypeCaller per interval (parallel)
- GenomicsDBImport per interval (parallel)
- GenotypeGVCFs per interval (parallel)
- merge per-interval raw VCF -> one merged raw VCF
- split SNP/INDEL hard filtering
- merge filtered SNP + INDEL into final VCF

#### `--use_bqsr true`

- per-interval bcftools seed call + merge seed VCF
- high-confidence SNP extraction
- **BaseRecalibrator + ApplyBQSR on full sample BAM (single whole-genome recalibration step)**
- then same per-interval GATK calling chain as above
- final merged filtered VCF

GATK outputs:

- per-interval intermediates:
  - `results/05_variant_calling/gatk/haplotypecaller/per_chrom/`
  - `results/05_variant_calling/gatk/genomicsdb/per_chrom/`
  - `results/05_variant_calling/gatk/genotypegvcfs/per_chrom/`
- merged raw VCF:
  - `results/05_variant_calling/gatk/genotypegvcfs/SAMPLE.gatk.raw.vcf.gz`
- final filtered VCF:
  - `results/05_variant_calling/gatk/filter/SAMPLE.gatk.filtered.vcf.gz`

## Input Naming Rules

Sample ID inference:

- PE: strip `_R1.fastq.gz` / `_R2.fastq.gz` or `_R1.fq.gz` / `_R2.fq.gz`
- SE: strip `.fastq.gz` or `.fq.gz`

## Parameters

## Required

- `--input_dir`: FASTQ directory
- `--ref`: reference FASTA

## Core control

- `--run_step`: `Raw_QC|Trimming_QC|Aligning|Calling|all` (default: `all`)
- `--read_type`: `SE|PE` (default: `PE`)
- `--caller`: `bcftools|gatk` (default: `bcftools`)
- `--use_bqsr`: `true|false` (default: `false`)
- `--fastp_parameters`: extra fastp parameter string
- `--help`: print help and exit

## Container

- `--sif`: container image path (default: `WGS_Variant_Calling.sif`)

## Resource controls

- `--threads`: global fallback threads (default: `4`)
- Tool-specific overrides:
  - `--fastqc_cpus`
  - `--fastp_cpus`
  - `--bwa_cpus`
  - `--bwamem2_cpus`
  - `--sambamba_cpus`
  - `--bcftools_cpus`
  - `--gatk_cpus`

## Typical Commands

## Print help

```bash
nextflow run main.nf --help
```

## Full pipeline

```bash
nextflow run main.nf \
  --input_dir /data/fastq \
  --ref /data/ref.fa \
  --read_type PE \
  --run_step all \
  --caller gatk \
  --use_bqsr false \
  --threads 8 \
  --sif WGS_Variant_Calling.sif
```

## Calling only (bcftools)

```bash
nextflow run main.nf \
  --input_dir /data/fastq \
  --ref /data/ref.fa \
  --read_type PE \
  --run_step Calling \
  --caller bcftools \
  --threads 4
```

## Calling only (gatk + BQSR)

```bash
nextflow run main.nf \
  --input_dir /data/fastq \
  --ref /data/ref.fa \
  --read_type PE \
  --run_step Calling \
  --caller gatk \
  --use_bqsr true \
  --threads 4
```



## Operational Notes

- The pipeline parses chromosome/scaffold names directly from FASTA headers.
- Step4 parallelism increases process count; on restricted systems, use conservative settings:
  - `--threads 1`
  - low scheduler parallelism (if configured externally)
- Some GATK subcommands can hit native-thread limits on tightly restricted nodes; the current implementation applies conservative Java/thread environment options in critical processes.

## Reproducibility

- All processes are containerized via `container "${params.sif}"`.
- Process reports/traces/timeline/DAG are enabled in `nextflow.config` under `pipeline_info/`.
- Work directory defaults to `work/`.


