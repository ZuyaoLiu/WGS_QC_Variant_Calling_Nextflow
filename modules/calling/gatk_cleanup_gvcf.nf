process GATK_CLEANUP_GVCF_PROCESS {
    tag "cleanup:${interval_idx}"
    container "${params.container_image}"

    input:
    tuple val(interval_idx), val(interval_spec), val(gvcf_manifest), val(gvcf_index_manifest), val(genomicsdb_dir)

    output:
    tuple val(interval_idx), val(interval_spec), path("${interval_idx}.gvcf_cleanup.done"), emit: cleanup_done

    script:
    """
    if [ ! -d "${genomicsdb_dir}" ]; then
      echo "Expected GenomicsDB directory not found: ${genomicsdb_dir}" >&2
      exit 1
    fi

    while IFS= read -r gvcf; do
      [ -n "\${gvcf}" ] && rm -f "\${gvcf}"
    done <<'EOF_GVCF'
${gvcf_manifest}
EOF_GVCF

    while IFS= read -r gvcf_index; do
      [ -n "\${gvcf_index}" ] && rm -f "\${gvcf_index}"
    done <<'EOF_GVCF_INDEX'
${gvcf_index_manifest}
EOF_GVCF_INDEX

    touch ${interval_idx}.gvcf_cleanup.done
    """
}

workflow GATK_CLEANUP_GVCF {
    take:
    ch_gvcf
    ch_gendb

    main:
    ch_grouped_gvcf = ch_gvcf
      .map { sample_id, interval_idx, interval_spec, gvcf, gvcf_index ->
        tuple(
            interval_idx,
            interval_spec,
            gvcf.toAbsolutePath().toString(),
            gvcf_index.toAbsolutePath().toString()
        )
      }
      .groupTuple()
      .map { interval_idx, interval_spec_list, gvcf_paths, gvcf_index_paths ->
        def uniqIntervals = interval_spec_list.unique()
        if (uniqIntervals.size() != 1) {
          error "Inconsistent interval labels for cleanup ${interval_idx}: ${uniqIntervals}"
        }
        tuple(
            interval_idx,
            uniqIntervals[0],
            gvcf_paths.join('\n'),
            gvcf_index_paths.join('\n')
        )
      }

    ch_ready_for_cleanup = ch_grouped_gvcf
      .join(
        ch_gendb.map { interval_idx, interval_spec, genomicsdb_dir ->
          tuple(interval_idx, interval_spec, genomicsdb_dir.toAbsolutePath().toString())
        }
      )
      .map { interval_idx, gvcf_interval_spec, gvcf_manifest, gvcf_index_manifest, gendb_interval_spec, genomicsdb_dir ->
        if (gvcf_interval_spec != gendb_interval_spec) {
          error "Interval mismatch during gVCF cleanup for ${interval_idx}: ${gvcf_interval_spec} vs ${gendb_interval_spec}"
        }
        tuple(interval_idx, gvcf_interval_spec, gvcf_manifest, gvcf_index_manifest, genomicsdb_dir)
      }

    GATK_CLEANUP_GVCF_PROCESS(ch_ready_for_cleanup)

    emit:
    cleanup_done = GATK_CLEANUP_GVCF_PROCESS.out.cleanup_done
}
