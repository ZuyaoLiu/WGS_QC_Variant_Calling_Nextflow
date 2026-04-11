process GATK_BQSR_PROCESS {
    tag "${sample_id}"
    container "${params.sif}"
    publishDir 'results/04_variant_calling/gatk/bqsr', mode: 'move'

    input:
    tuple val(sample_id), path(markdup_bam), path(markdup_bai), path(bqsr_panel), path(ref_fa)

    output:
    tuple val(sample_id), path("${sample_id}.recal.bam"), path("${sample_id}.recal.bam.bai"), emit: recal_bam

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
      -I ${markdup_bam} \
      --known-sites ${bqsr_panel} \
      ${params.gatk_baserecalibrator_parameters} \
      -O ${sample_id}.recal.table

    gatk ApplyBQSR \
      --java-options "-Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=false -Dsamjdk.use_async_io_write_tribble=false -Xmx2g" \
      -R ${ref_fa} \
      -I ${markdup_bam} \
      --bqsr-recal-file ${sample_id}.recal.table \
      ${params.gatk_applybqsr_parameters} \
      -O ${sample_id}.recal.bam

    samtools index -@ ${task.cpus} ${sample_id}.recal.bam
    """
}

workflow GATK_BQSR {
    take:
    ch_markdup
    ch_bqsr_panel
    ch_ref

    main:
    // Apply the same external known-sites panel to each sample BAM for per-sample recalibration.
    ch_joined = ch_markdup.combine(ch_bqsr_panel).map { row ->
        tuple(row[0], row[1], row[2], row[3])
    }
    ch_with_ref = ch_joined.combine(ch_ref).map { row -> tuple(row[0], row[1], row[2], row[3], row[4]) }
    GATK_BQSR_PROCESS(ch_with_ref)

    emit:
    recal_bam = GATK_BQSR_PROCESS.out.recal_bam.map { sample_id, bam, bai ->
        tuple(
            sample_id,
            file("results/04_variant_calling/gatk/bqsr/${bam.name}"),
            file("results/04_variant_calling/gatk/bqsr/${bai.name}")
        )
    }
}
