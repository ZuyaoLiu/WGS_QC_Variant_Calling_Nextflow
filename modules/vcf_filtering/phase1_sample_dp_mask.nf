process PHASE1_SAMPLE_DP_MASK_PROCESS {
    tag "${callset_id}"
    container "${params.container_image}"
    publishDir 'results/05_vcf_filtering/01_phase1/intermediate', mode: 'copy'

    input:
    tuple val(callset_id), path(input_bcf), path(input_csi)
    path(dp_mask_script)

    output:
    tuple val(callset_id), path("${callset_id}.vcffilter.phase1.bcf"), path("${callset_id}.vcffilter.phase1.bcf.csi"), emit: masked_vcf
    tuple val(callset_id), path("${callset_id}.phase1.sample_mean_dp.tsv"), emit: sample_table

    script:
    """
    python ${dp_mask_script} full \
      --input-vcf ${input_bcf} \
      --output-bcf ${callset_id}.vcffilter.phase1.bcf \
      --sample-table ${callset_id}.phase1.sample_mean_dp.tsv \
      --min-factor ${params.phase1_sample_dp_min_factor} \
      --max-factor ${params.phase1_sample_dp_max_factor} \
      --threads ${task.cpus}

    bcftools index --threads ${task.cpus} -f ${callset_id}.vcffilter.phase1.bcf
    """
}

workflow PHASE1_SAMPLE_DP_MASK {
    take:
    ch_filtered_vcf

    main:
    dp_mask_script = file("${projectDir}/vcf_filtering_scripts/sample_dp_mask.py")
    PHASE1_SAMPLE_DP_MASK_PROCESS(ch_filtered_vcf, dp_mask_script)

    emit:
    masked_vcf = PHASE1_SAMPLE_DP_MASK_PROCESS.out.masked_vcf
    sample_table = PHASE1_SAMPLE_DP_MASK_PROCESS.out.sample_table
}
