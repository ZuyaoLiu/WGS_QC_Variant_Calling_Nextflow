process GATK_GENOMICSDBIMPORT_BY_CHR_PROCESS {
    tag "cohort:${chrom}"
    container "${params.sif}"
    cpus { (params.gatk_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/05_variant_calling/gatk/genomicsdb/per_chrom', mode: 'copy'

    input:
    tuple val(interval_idx), val(chrom), path(gvcf_files), path(gvcf_indexes), path(ref_fa)

    output:
    tuple val(interval_idx), val(chrom), path("genomicsdb_${interval_idx}"), emit: gendb_per_chr

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

    ls -1 ${gvcf_files} | sort > gvcf.list
    mapfile -t GVCFS < gvcf.list
    VAR_ARGS=""
    for g in "\${GVCFS[@]}"; do
      VAR_ARGS="\${VAR_ARGS} --variant \${g}"
    done

    gatk GenomicsDBImport \
      --java-options "-Xmx2g -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:ActiveProcessorCount=1 -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=false -Dsamjdk.use_async_io_write_tribble=false" \
      \${VAR_ARGS} \
      --genomicsdb-workspace-path genomicsdb_${interval_idx} \
      --intervals ${chrom} \
      --batch-size 1 \
      ${params.gatk_genomicsdbimport_parameters} \
      --reader-threads ${task.cpus}
    """
}

workflow GATK_GENOMICSDBIMPORT {
    take:
    ch_gvcf
    ch_ref

    main:
    ch_grouped = ch_gvcf
      .map { sample_id, interval_idx, chrom, gvcf, gvcf_index -> tuple(interval_idx, chrom, gvcf, gvcf_index) }
      .groupTuple()
      .map { interval_idx, chrom_list, gvcf_list, gvcf_index_list ->
        def uniq_chrom = chrom_list.unique()
        if (uniq_chrom.size() != 1) {
          error "Inconsistent chromosome labels for interval ${interval_idx}: ${uniq_chrom}"
        }
        tuple(interval_idx, uniq_chrom[0], gvcf_list, gvcf_index_list)
      }

    ch_input = ch_grouped.combine(ch_ref).map { row ->
      tuple(row[0], row[1], row[2], row[3], row[4])
    }

    GATK_GENOMICSDBIMPORT_BY_CHR_PROCESS(ch_input)

    emit:
    gendb = GATK_GENOMICSDBIMPORT_BY_CHR_PROCESS.out.gendb_per_chr
}
