process ALIGN_MARKDUP_PE_PROCESS {
    tag "${sample_id}"
    container "${params.container_image}"
    publishDir 'results/03_markdup', mode: 'move'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path(ref_bundle)
    path(ref_fa)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.cram"), path("${sample_id}.markdup.cram.crai"), emit: markdup

    script:
    """
    mkdir -p tmp_markdup

    bwa-mem2 mem -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" ${params.bwamem2_parameters} ${ref_fa} ${r1} ${r2} \
      | samtools collate -@ ${task.cpus} -u -O -T collate_tmp - \
      | samtools fixmate -m -@ ${task.cpus} -u - - \
      | samtools sort -@ ${task.cpus} -o ${sample_id}.coordsort.bam -

    sambamba markdup -t ${task.cpus} --tmpdir=tmp_markdup ${sample_id}.coordsort.bam ${sample_id}.markdup.bam
    rm -f ${sample_id}.coordsort.bam

    samtools view -@ ${task.cpus} -C -T ${ref_fa} --output-fmt-option version=3.0,embed_ref=1 -o ${sample_id}.markdup.cram ${sample_id}.markdup.bam
    rm -f ${sample_id}.markdup.bam

    samtools index -@ ${task.cpus} ${sample_id}.markdup.cram
    """
}

workflow ALIGN_MARKDUP_PE {
    take:
    ch_input
    ch_ref_bundle
    ch_ref_fa

    main:
    ALIGN_MARKDUP_PE_PROCESS(ch_input, ch_ref_bundle, ch_ref_fa)

    emit:
    markdup_alignment = ALIGN_MARKDUP_PE_PROCESS.out.markdup.map { sample_id, cram, crai ->
        tuple(
            sample_id,
            file("results/03_markdup/${cram.name}"),
            file("results/03_markdup/${crai.name}")
        )
    }
}
