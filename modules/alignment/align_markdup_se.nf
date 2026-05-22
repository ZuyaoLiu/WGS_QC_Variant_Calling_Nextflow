process ALIGN_MARKDUP_SE_PROCESS {
    tag "${sample_id}"
    container "${params.container_image}"
    publishDir 'results/03_markdup', mode: 'symlink', enabled: params.publish_cram

    input:
    tuple val(sample_id), path(fq)
    path(ref_bundle)
    path(ref_fa)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.cram"), path("${sample_id}.markdup.cram.crai"), emit: markdup

    script:
    """
    mkdir -p tmp_markdup

    bwa-mem2 mem -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" ${params.bwamem2_parameters} ${ref_fa} ${fq} \
      | samtools sort -@ ${task.cpus} -o ${sample_id}.coordsort.bam -

    sambamba markdup -t ${task.cpus} --tmpdir=tmp_markdup ${sample_id}.coordsort.bam ${sample_id}.markdup.bam
    rm -f ${sample_id}.coordsort.bam

    samtools view -@ ${task.cpus} -C -T ${ref_fa} --output-fmt-option version=3.0,embed_ref=1 -o ${sample_id}.markdup.cram ${sample_id}.markdup.bam
    rm -f ${sample_id}.markdup.bam

    samtools index -@ ${task.cpus} ${sample_id}.markdup.cram
    """
}

workflow ALIGN_MARKDUP_SE {
    take:
    ch_input
    ch_ref_bundle
    ch_ref_fa

    main:
    ALIGN_MARKDUP_SE_PROCESS(ch_input, ch_ref_bundle, ch_ref_fa)

    emit:
    markdup_alignment = ALIGN_MARKDUP_SE_PROCESS.out.markdup
}
