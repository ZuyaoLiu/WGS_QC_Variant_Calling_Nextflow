process GATK_GENOMICSDBIMPORT_BY_CHR_PROCESS {
    tag "cohort:${interval_idx}"
    container "${params.container_image}"

    input:
    tuple val(interval_idx), val(interval_spec), path(gvcf_files), path(gvcf_indexes), path(ref_fa), path(ref_fai), path(ref_dict)

    output:
    tuple val(interval_idx), val(interval_spec), path("genomicsdb_${interval_idx}"), emit: gendb_per_chr

    script:
    // Leave ~4 GB of task.memory for native/off-heap allocations; floor at 2 GB.
    def xmx_g = Math.max(2, (task.memory ? task.memory.toGiga() - 4 : 2))
    def consolidate_flag = params.gatk_genomicsdbimport_consolidate ? '--consolidate' : ''
    """
    ls -1 ${gvcf_files} | sort > gvcf.list
    mapfile -t GVCFS < gvcf.list
    VAR_ARGS=""
    for g in "\${GVCFS[@]}"; do
      VAR_ARGS="\${VAR_ARGS} --variant \${g}"
    done

    gatk GenomicsDBImport \
      --java-options "-Xmx${xmx_g}g -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=false -Dsamjdk.use_async_io_write_tribble=false" \
      \${VAR_ARGS} \
      --genomicsdb-workspace-path genomicsdb_${interval_idx} \
      --intervals ${interval_spec} \
      --batch-size ${params.gatk_genomicsdbimport_batch_size} \
      ${consolidate_flag} \
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
      .map { sample_id, interval_idx, interval_spec, gvcf, gvcf_index -> tuple(interval_idx, interval_spec, gvcf, gvcf_index) }
      .groupTuple()
      .map { interval_idx, interval_spec_list, gvcf_list, gvcf_index_list ->
        def uniqIntervals = interval_spec_list.unique()
        if (uniqIntervals.size() != 1) {
          error "Inconsistent interval labels for interval ${interval_idx}: ${uniqIntervals}"
        }
        tuple(interval_idx, uniqIntervals[0], gvcf_list, gvcf_index_list)
      }

    // ch_ref is a value channel of tuple(ref_fa, ref_fai, ref_dict) from PREPARE_REFERENCE.
    ch_input = ch_grouped.combine(ch_ref).map { row ->
      tuple(row[0], row[1], row[2], row[3], row[4], row[5], row[6])
    }

    GATK_GENOMICSDBIMPORT_BY_CHR_PROCESS(ch_input)

    emit:
    gendb = GATK_GENOMICSDBIMPORT_BY_CHR_PROCESS.out.gendb_per_chr
}
