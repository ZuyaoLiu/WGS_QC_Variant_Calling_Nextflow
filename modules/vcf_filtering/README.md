# Step5 VCF Filtering Modules

This directory contains DSL2 modules for the optional `vcffilter` stage.

Planned structure:

- `phase1_*`
- `phase2_*`
- `*_stats`
- `*_plots`

Current status:

- workflow entry exists in `subworkflows/vcf_filtering.nf`
- concrete filtering and reporting modules will be added incrementally

Scope limitation:

- current design supports `biallelic SNP` only
- `indel` and mixed `SNP + indel` inputs are out of scope for this module
