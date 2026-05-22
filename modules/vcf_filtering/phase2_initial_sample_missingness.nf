process PHASE2_INITIAL_SAMPLE_MISSINGNESS_PROCESS {
    tag "${callset_id}"
    container "${params.container_image}"
    publishDir 'results/05_vcf_filtering/02_phase2/intermediate', mode: 'copy'

    input:
    tuple val(callset_id), path(input_vcf), path(input_tbi)
    path(compute_script)

    output:
    tuple val(callset_id), path("${callset_id}.initial.sample_missingness.tsv"), emit: sample_missingness

    script:
    """
    python ${compute_script} \
      --input-vcf ${input_vcf} \
      --output-tsv ${callset_id}.initial.sample_missingness.tsv \
      --threads ${task.cpus}
    """
}

workflow PHASE2_INITIAL_SAMPLE_MISSINGNESS {
    take:
    ch_called_vcf

    main:
    compute_script = file("${projectDir}/vcf_filtering_scripts/compute_sample_missingness.py")
    PHASE2_INITIAL_SAMPLE_MISSINGNESS_PROCESS(ch_called_vcf, compute_script)

    emit:
    sample_missingness = PHASE2_INITIAL_SAMPLE_MISSINGNESS_PROCESS.out.sample_missingness
}
