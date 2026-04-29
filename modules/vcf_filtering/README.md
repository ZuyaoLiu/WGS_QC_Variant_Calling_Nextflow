# Step5 VCF Filtering Modules

This directory contains DSL2 modules for the optional `vcffilter` stage.

Implemented modules:

- `phase1_initial_filter_mask.nf`
- `phase1_sample_dp_mask.nf`
- `phase1_reports.nf`
- `phase2_remove_high_missing_samples.nf`
- `phase2_site_mean_dp_limits.nf`
- `phase2_final_filter.nf`
- `phase2_final_reports.nf`

Workflow behavior:

- the workflow entry is `subworkflows/vcf_filtering.nf`
- phase1 always runs first and creates review outputs
- `--vcffilter_stop_after phase1` stops after phase1
- default behavior continues to phase2 and writes the final filtered VCF

Scope limitation:

- current design supports `biallelic SNP` only
- `indel` and mixed `SNP + indel` inputs are out of scope for this module
