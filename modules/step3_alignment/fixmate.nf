process FIXMATE_PROCESS {
    tag "${sample_id}"
    container "${params.sif}"
    cpus params.threads
    publishDir 'results/03_align', mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.fixmate.bam"), emit: bam

    script:
    """
    samtools fixmate -m -@ ${task.cpus} ${bam} ${sample_id}.fixmate.bam
    """
}

workflow FIXMATE {
    take:
    ch_input

    main:
    FIXMATE_PROCESS(ch_input)

    emit:
    bam_for_sort = FIXMATE_PROCESS.out.bam
}
