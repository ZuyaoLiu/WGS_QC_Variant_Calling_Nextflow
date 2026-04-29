# VCF Filtering Scripts

This directory contains helper scripts used by the `vcffilter` subworkflow.

Current scope:

- supports `biallelic SNP` filtering only
- does not support `indel`
- does not support mixed `SNP + indel` workflows

Scripts:

- `sample_dp_mask.py`
  Computes per-sample DP thresholds and masks genotypes outside the allowed range.
  Supports both whole-file and region-parallel execution.

- `plot_sample_dp_distribution.R`
  Plots histogram and density figures from a per-sample DP summary table.

- `plot_sample_missingness.R`
  Plots histogram, density, and cutoff-vs-samples-kept figures from a per-sample missingness table.

Notes:

- These scripts are intended to be called from Nextflow processes.
- New `vcffilter` code should use the helpers in this directory.
