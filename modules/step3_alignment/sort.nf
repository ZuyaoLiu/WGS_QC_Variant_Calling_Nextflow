process SORT_BAM_PROCESS {
    tag "${sample_id}"
    container "${params.sif}"
    cpus { (params.threads ?: 1) as Integer }
    publishDir 'results/03_align', mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    """
}

workflow SORT_BAM {
    take:
    ch_input

    main:
    SORT_BAM_PROCESS(ch_input)

    emit:
    sorted_bam = SORT_BAM_PROCESS.out.bam
}
