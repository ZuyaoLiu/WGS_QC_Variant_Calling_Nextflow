process ALIGN_MARKDUP_SE_PROCESS {
    tag "${sample_id}"
    container "${params.sif}"
    publishDir 'results/03_markdup', mode: 'move'

    input:
    tuple val(sample_id), path(fq)
    path(ref_bundle)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), path("${sample_id}.markdup.bam.bai"), emit: markdup

    script:
    """
    bwa mem -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" ${params.bwa_parameters} ${params.ref_base} ${fq} \
      | samtools sort -u -@ ${task.cpus} - \
      | samtools markdup -r -@ ${task.cpus} - ${sample_id}.markdup.bam

    samtools index -@ ${task.cpus} ${sample_id}.markdup.bam
    """
}

workflow ALIGN_MARKDUP_SE {
    take:
    ch_input
    ch_ref_bundle

    main:
    ALIGN_MARKDUP_SE_PROCESS(ch_input, ch_ref_bundle)

    emit:
    markdup_bam = ALIGN_MARKDUP_SE_PROCESS.out.markdup.map { sample_id, bam, bai ->
        tuple(
            sample_id,
            file("results/03_markdup/${bam.name}"),
            file("results/03_markdup/${bai.name}")
        )
    }
}
