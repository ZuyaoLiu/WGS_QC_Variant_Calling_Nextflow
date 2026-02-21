process BWA_SE_PROCESS {
    tag "${sample_id}"
    container "${params.sif}"
    cpus params.bwa_cpus
    publishDir 'results/03_align', mode: 'copy'

    input:
    tuple val(sample_id), path(fq)
    path(ref_bundle)

    output:
    tuple val(sample_id), path("${sample_id}.unsorted.bam"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" ${params.ref_base} ${fq} | samtools view -@ ${task.cpus} -Sb -o ${sample_id}.unsorted.bam -
    """
}

workflow BWA_SE {
    take:
    ch_input
    ch_ref_bundle

    main:
    BWA_SE_PROCESS(ch_input, ch_ref_bundle)

    emit:
    bam_for_sort = BWA_SE_PROCESS.out.bam
}
