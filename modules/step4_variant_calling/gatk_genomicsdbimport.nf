process GATK_GENOMICSDBIMPORT_BY_CHR_PROCESS {
    tag "${sample_id}:${chrom}"
    container "${params.sif}"
    cpus params.gatk_cpus
    publishDir 'results/05_variant_calling/gatk/genomicsdb/per_chrom', mode: 'copy'

    input:
    tuple val(sample_id), val(interval_idx), val(chrom), path(gvcf), path(gvcf_index), path(ref_fa)

    output:
    tuple val(sample_id), val(interval_idx), val(chrom), path("genomicsdb_${sample_id}_${interval_idx}"), emit: gendb_per_chr

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

    gatk GenomicsDBImport \
      --java-options "-Xmx2g -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:ActiveProcessorCount=1 -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=false -Dsamjdk.use_async_io_write_tribble=false" \
      --variant ${gvcf} \
      --genomicsdb-workspace-path genomicsdb_${sample_id}_${interval_idx} \
      --intervals ${chrom} \
      --batch-size 1 \
      --reader-threads ${task.cpus}
    """
}

workflow GATK_GENOMICSDBIMPORT {
    take:
    ch_gvcf
    ch_ref

    main:
    ch_input = ch_gvcf.combine(ch_ref).map { row ->
      tuple(row[0], row[1], row[2], row[3], row[4], row[5])
    }

    GATK_GENOMICSDBIMPORT_BY_CHR_PROCESS(ch_input)

    emit:
    gendb = GATK_GENOMICSDBIMPORT_BY_CHR_PROCESS.out.gendb_per_chr
}
