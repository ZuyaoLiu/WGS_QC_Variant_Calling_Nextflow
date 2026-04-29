# WGS QC and Variant Calling Pipeline

**Version:** 3.0.0  
**Workflow engine:** Nextflow DSL2  
**Scope:** short-read, non-human whole-genome QC, alignment, variant calling, and optional SNP filtering.

![Pipeline workflow](./Pipeline_figure.png)

## Overview

This pipeline runs a modular WGS analysis workflow:

```text
Raw_QC -> Trimming_QC -> Aligning -> Calling -> optional VCF_Filtering
```

Main features:

- Raw FASTQ QC with FastQC and MultiQC
- Adapter/quality trimming with fastp
- Alignment with bwa-mem2 and duplicate marking
- Reference-embedded CRAM output
- Cohort variant calling with `bcftools` or `GATK`
- Optional post-calling VCF filtering for biallelic SNPs
- Local, SLURM, and AWS Batch execution profiles

## Quickstart

### 1. Prepare inputs

Paired-end FASTQ files must use this naming pattern:

```text
sampleA_R1.fastq.gz
sampleA_R2.fastq.gz
sampleB_R1.fastq.gz
sampleB_R2.fastq.gz
```

Single-end FASTQ files can use:

```text
sampleA.fastq.gz
sampleB.fq.gz
```

If using SRA input, provide a tab-delimited file:

```text
SRRxxxxxx    sampleA
SRRyyyyyy    sampleB
```

### 2. Run the full pipeline

```bash
nextflow run main.nf \
  -profile local \
  --input_dir /path/to/fastq \
  --ref /path/to/reference.fa \
  --read_type PE \
  --run_step all \
  --caller bcftools \
  --interval_size 5M
```

### 3. Check outputs

Main results are written to:

```text
results/
pipeline_info/
```

For a `bcftools` run, the final cohort VCF is:

```text
results/04_variant_calling/bcftools/cohort.bcftools.vcf.gz
```

## Common Run Modes

Run raw QC only:

```bash
nextflow run main.nf \
  -profile local \
  --input_dir /path/to/fastq \
  --ref /path/to/reference.fa \
  --read_type PE \
  --run_step Raw_QC
```

Start from trimmed FASTQ files:

```bash
nextflow run main.nf \
  -profile local \
  --trimmed_input_dir /path/to/trimmed_fastq \
  --ref /path/to/reference.fa \
  --read_type PE \
  --run_step Aligning \
  --caller bcftools \
  --interval_size 5M
```

Run GATK with BQSR:

```bash
nextflow run main.nf \
  -profile local \
  --input_dir /path/to/fastq \
  --ref /path/to/reference.fa \
  --read_type PE \
  --run_step all \
  --caller gatk \
  --interval_size 5M \
  --use_bqsr true \
  --bqsr_panel /path/to/known_sites.vcf.gz
```

Run standalone VCF filtering:

```bash
nextflow run main.nf \
  -profile local \
  --run_step VCF_Filtering \
  --input_vcf /path/to/cohort.vcf.gz
```

Standalone VCF filtering requires an index next to the input VCF:

```text
/path/to/cohort.vcf.gz.tbi
```

## Parameters

### Required and Common Parameters

| Parameter | Default | Allowed values | Description |
|---|---:|---|---|
| `--run_step` | `all` | `Raw_QC`, `Trimming_QC`, `Aligning`, `Calling`, `VCF_Filtering`, `all` | Start point of the workflow. Later steps continue automatically except `Raw_QC`, which runs QC only. |
| `--input_dir` | `null` | directory | Raw FASTQ directory. Required for `Raw_QC`, `Trimming_QC`, and `all` unless `--SRA_list` is used for raw entry. |
| `--SRA_list` | `null` | file | Two-column tab-delimited SRA accession and sample name file. Use instead of `--input_dir` for raw entry. |
| `--trimmed_input_dir` | `null` | directory | Clean FASTQ directory for starting from `Aligning` or as a one-step fallback for `Calling`. |
| `--input_vcf` | `null` | `.vcf.gz` file | Input VCF for standalone `VCF_Filtering`. Requires `INPUT.vcf.gz.tbi`. |
| `--ref` | `null` | FASTA file | Reference genome. Required except for standalone `VCF_Filtering`. |
| `--read_type` | `PE` | `PE`, `SE` | Paired-end or single-end input mode. |
| `--caller` | `bcftools` | `bcftools`, `gatk` | Variant caller backend. |
| `--interval_size` | `null` | e.g. `500K`, `5M`, `1G`, `5000000` | Interval block size for parallel variant calling. Required for any run that executes `Calling`. |
| `--enable_vcffilter` | `false` | `true`, `false` | Run post-calling VCF filtering after variant calling. |
| `--vcffilter_stop_after` | `final` | `phase1`, `final` | Stop VCF filtering after phase1 review or continue to final phase2 output. |

### GATK and BQSR Parameters

| Parameter | Default | Allowed values | Description |
|---|---:|---|---|
| `--use_bqsr` | `false` | `true`, `false` | Enable GATK BaseRecalibrator and ApplyBQSR. Only valid with `--caller gatk`. |
| `--bqsr_panel` | `null` | indexed VCF/VCF.GZ | Known-sites panel for BQSR. Required when `--use_bqsr true`. |
| `--cleanup_GVCFS` | `false` | `true`, `false` | Remove interval gVCFs after successful GenomicsDBImport. |

### Container Parameters

| Parameter | Default | Description |
|---|---:|---|
| `--image` | `leonardliu0910/wgs_variant_calling:latest` | Docker or registry image used by process containers. |
| `--sif` | `null` | Local Singularity/Apptainer image. Overrides `--image` when set. |

### Advanced Tool Options

These parameters pass extra options directly to the underlying tools. Use quotes.

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
| `--gatk_variantfiltration_snp_parameters` | `GATK VariantFiltration` for SNPs |
| `--gatk_variantfiltration_indel_parameters` | `GATK VariantFiltration` for indels |

### VCF Filtering Parameters

VCF filtering defaults are defined in `vcffilter.config`.

| Parameter | Default | Description |
|---|---:|---|
| `--phase1_remove_mislabeled` | `false` | Remove project-specific mislabeled samples before genotype masking. |
| `--phase1_mislabeled_samples` | `null` | Optional file listing samples to remove when `--phase1_remove_mislabeled true`. |
| `--phase1_min_gq` | `20` | Mask genotypes with GQ below this value. |
| `--phase1_min_dp` | `5` | Mask genotypes with DP below this value. |
| `--phase1_enable_ab_mask` | `true` | Enable heterozygous allele-balance masking. |
| `--phase1_sample_dp_min_factor` | `0.5` | Lower per-sample DP factor for genotype masking. |
| `--phase1_sample_dp_max_factor` | `2.5` | Upper per-sample DP factor for genotype masking. |
| `--phase1_sample_missing_plot_cutoffs_pct` | `5,10,20,40` | Cutoffs used in phase1 sample-missingness plots. |
| `--phase2_remove_high_missing_samples` | `true` | Remove samples listed manually or failing missingness threshold. |
| `--phase2_high_missing_sample_file` | `null` | Optional user-provided sample removal list. |
| `--phase2_max_sample_missing` | `null` | Maximum allowed sample missingness. User-defined. |
| `--phase2_min_qual` | `30` | Minimum site QUAL for final filtering. |
| `--phase2_max_site_missing` | `0.2` | Maximum site missingness for final filtering. |
| `--phase2_enable_site_mean_dp_filter` | `true` | Filter sites using mean depth limits. |
| `--phase2_site_mean_dp_min_factor` | `0.5` | Lower site mean-DP factor. |
| `--phase2_site_mean_dp_max_factor` | `2.5` | Upper site mean-DP factor. |
| `--phase2_sample_missing_plot_cutoffs_pct` | `5,10,20,40` | Cutoffs used in final sample-missingness plots. |

## Execution Profiles

| Profile | Executor | Notes |
|---|---|---|
| `local` | local | Default local profile. Apptainer is enabled in `nextflow.config`. |
| `slurm` | SLURM | Uses the SLURM resource defaults and `executor.queueSize = 4` from `nextflow.config`. |
| `awsbatch` | AWS Batch | Template profile. Edit queue, region, container, and S3 work directory before use. |

Resource tuning should be done in `nextflow.config`, preferably through profile-level defaults or `withName` process overrides.

## Outputs

```text
results/
├── 00_sra_download/
├── 01_raw_qc/
│   ├── fastqc/
│   └── multiqc/
├── 02_fastp_qc/
│   ├── fastp/
│   ├── fastqc/
│   └── multiqc/
├── 03_markdup/
│   ├── SAMPLE.markdup.cram
│   └── SAMPLE.markdup.cram.crai
├── 04_variant_calling/
│   ├── bcftools/
│   └── gatk/
└── 05_vcf_filtering/
    ├── 01_phase1/
    ├── 02_phase2/
    └── final/
```

Workflow run reports are written to:

```text
pipeline_info/
├── timeline.html
├── report.html
├── trace.txt
└── dag.html
```

## Step Behavior

| `--run_step` | Behavior |
|---|---|
| `Raw_QC` | Runs raw FastQC and MultiQC only, then exits. |
| `Trimming_QC` | Runs trimming, alignment, and calling. VCF filtering runs only if `--enable_vcffilter true`. |
| `Aligning` | Starts from clean FASTQ input or reuses `results/02_fastp_qc/fastp` when available. |
| `Calling` | Reuses `results/03_markdup` when complete, otherwise falls back one step to alignment. |
| `VCF_Filtering` | Runs VCF filtering only from `--input_vcf`. |
| `all` | Runs the full workflow from raw FASTQ or SRA input through variant calling. |

## Important Notes

- Reference sequence names must not use `chr:start-end` style names, because this conflicts with region parsing.
- VCF filtering currently supports biallelic SNP filtering only.
- Indel-only and mixed SNP+indel VCF filtering are outside the current filtering scope.
- For `--caller gatk --use_bqsr true`, the BQSR panel must be built on the same reference genome.
- Use `nextflow run main.nf --help` or `nextflow run main.nf --help_full` for command-line help.

## Repository Layout

```text
.
├── main.nf
├── nextflow.config
├── vcffilter.config
├── subworkflows/
├── modules/
├── vcf_filtering_scripts/
├── Pipeline_figure.png
└── README.md
```
