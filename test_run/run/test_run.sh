nextflow run ../../main.nf \
  -profile local \
  --input_dir ../data \
  --ref ../data/sim_ref_100kb.fa \
  --read_type PE \
  --run_step all \
  --caller gatk \
  --use_bqsr false \
  --threads 4 \
  --sif WGS_Variant_Calling.sif
