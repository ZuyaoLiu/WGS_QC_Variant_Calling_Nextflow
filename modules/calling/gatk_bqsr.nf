process GATK_BQSR_PROCESS {
    tag "${sample_id}"
    container "${params.container_image}"
    publishDir 'results/04_variant_calling/gatk/bqsr', mode: 'move'

    input:
    tuple val(sample_id), path(markdup_alignment), path(markdup_index), path(bqsr_panel), path(ref_fa)

    output:
    tuple val(sample_id), path("${sample_id}.recal.cram"), path("${sample_id}.recal.cram.crai"), emit: recal_alignment

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export POLARS_MAX_THREADS=${task.cpus}
    export RAYON_NUM_THREADS=${task.cpus}

    if [ ! -f ${ref_fa}.fai ]; then
      samtools faidx ${ref_fa}
    fi

    if [ ! -f ${ref_fa.baseName}.dict ]; then
      gatk CreateSequenceDictionary -R ${ref_fa} -O ${ref_fa.baseName}.dict
    fi

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
    ch_with_ref = ch_joined.combine(ch_ref).map { row -> tuple(row[0], row[1], row[2], row[3], row[4]) }
    GATK_BQSR_PROCESS(ch_with_ref)

    emit:
    recal_alignment = GATK_BQSR_PROCESS.out.recal_alignment.map { sample_id, cram, crai ->
        tuple(
            sample_id,
            file("results/04_variant_calling/gatk/bqsr/${cram.name}"),
            file("results/04_variant_calling/gatk/bqsr/${crai.name}")
        )
    }
}
