process GATK_BQSR_PROCESS {
    tag "${sample_id}"
    container "${params.container_image}"

    input:
    tuple val(sample_id), path(markdup_alignment), path(markdup_index), path(bqsr_panel), path(ref_fa), path(ref_fai), path(ref_dict)

    output:
    tuple val(sample_id), path("${sample_id}.recal.cram"), path("${sample_id}.recal.cram.crai"), emit: recal_alignment

    script:
    """
    gatk BaseRecalibrator \
      --java-options "-Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=false -Dsamjdk.use_async_io_write_tribble=false -Xmx2g" \
      -R ${ref_fa} \
      -I ${markdup_alignment} \
      --known-sites ${bqsr_panel} \
      ${params.gatk_baserecalibrator_parameters} \
      -O ${sample_id}.recal.table

    gatk ApplyBQSR \
      --java-options "-Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=false -Dsamjdk.use_async_io_write_tribble=false -Xmx2g" \
      -R ${ref_fa} \
      -I ${markdup_alignment} \
      --bqsr-recal-file ${sample_id}.recal.table \
      ${params.gatk_applybqsr_parameters} \
      -O ${sample_id}.recal.cram

    samtools index -@ ${task.cpus} ${sample_id}.recal.cram
    """
}

workflow GATK_BQSR {
    take:
    ch_markdup
    ch_bqsr_panel
    ch_ref

    main:
    // Apply the same external known-sites panel to each sample alignment file for per-sample recalibration.
    ch_joined = ch_markdup.combine(ch_bqsr_panel).map { row ->
        tuple(row[0], row[1], row[2], row[3])
    }
    // ch_ref is a value channel of tuple(ref_fa, ref_fai, ref_dict) from PREPARE_REFERENCE.
    ch_with_ref = ch_joined.combine(ch_ref).map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5], row[6]) }
    GATK_BQSR_PROCESS(ch_with_ref)

    emit:
    recal_alignment = GATK_BQSR_PROCESS.out.recal_alignment
}
