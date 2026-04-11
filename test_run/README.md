# test_run scaffold

This directory is prepared for local PE smoke testing with a larger synthetic dataset.

Structure:

- `data/`: synthetic paired-end FASTQ.gz files and a 100 kb reference FASTA
- `run/`: working directory for future test execution

Files in `data/`:

- `sim_ref_100kb.fa`
- `SimA_R1.fastq.gz`
- `SimA_R2.fastq.gz`
- `SimB_R1.fastq.gz`
- `SimB_R2.fastq.gz`
- `SimC_R1.fastq.gz`
- `SimC_R2.fastq.gz`

Dataset notes:

- 3 paired-end samples: `SimA`, `SimB`, `SimC`
- reference length: about `100 kb`
- read length: `150 bp`
- target depth per sample: about `20x`

Run future tests from `test_run/run/`.

Scheduler notes:

- Use `-profile local` for local smoke tests.
- If you want to test with `-profile slurm` or `-profile awsbatch`, edit the scheduler defaults in `nextflow.config` first.
