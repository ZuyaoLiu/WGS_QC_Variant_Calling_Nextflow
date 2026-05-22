process PHASE2_SITE_MEAN_DP_LIMITS_PROCESS {
    tag "${callset_id}"
    container "${params.container_image}"
    publishDir 'results/05_vcf_filtering/02_phase2/intermediate', mode: 'copy'

    input:
    tuple val(callset_id), path(input_vcf), path(input_tbi)
    path(dp_mask_script)

    output:
    tuple val(callset_id), path("${callset_id}.phase2.sample_mean_dp.tsv"), path("${callset_id}.phase2.site_mean_dp_limits.tsv"), emit: dp_limits

    script:
    """
    python ${dp_mask_script} compute-stats \
      --input-vcf ${input_vcf} \
      --output-stats ${callset_id}.phase2.sample_dp.stats.tsv \
      --threads ${task.cpus}

    python ${dp_mask_script} aggregate-stats \
      --stats-file ${callset_id}.phase2.sample_dp.stats.tsv \
      --sample-table ${callset_id}.phase2.sample_mean_dp.tsv \
      --min-factor ${params.phase2_site_mean_dp_min_factor} \
      --max-factor ${params.phase2_site_mean_dp_max_factor}

    if [[ "${params.phase2_enable_site_mean_dp_filter}" == "true" ]]; then
      global_mean_dp=\$(awk '
        NR == 1 { next }
        \$2 ~ /^[0-9.]+\$/ { sum += \$2; n += 1 }
        END {
          if (n > 0) {
            printf "%.6f\\n", sum / n
          } else {
            print "NA"
          }
        }
      ' ${callset_id}.phase2.sample_mean_dp.tsv)

      if [[ "\${global_mean_dp}" == "NA" ]]; then
        min_site_mean_dp="NA"
        max_site_mean_dp="NA"
      else
        min_site_mean_dp=\$(awk -v x="\${global_mean_dp}" -v f="${params.phase2_site_mean_dp_min_factor}" 'BEGIN { printf "%.6f\\n", x * f }')
        max_site_mean_dp=\$(awk -v x="\${global_mean_dp}" -v f="${params.phase2_site_mean_dp_max_factor}" 'BEGIN { printf "%.6f\\n", x * f }')
      fi
    else
      global_mean_dp="NA"
      min_site_mean_dp="NA"
      max_site_mean_dp="NA"
    fi

    {
      printf "metric\\tvalue\\n"
      printf "global_mean_dp\\t%s\\n" "\${global_mean_dp}"
      printf "min_site_mean_dp\\t%s\\n" "\${min_site_mean_dp}"
      printf "max_site_mean_dp\\t%s\\n" "\${max_site_mean_dp}"
    } > ${callset_id}.phase2.site_mean_dp_limits.tsv
    """
}

workflow PHASE2_SITE_MEAN_DP_LIMITS {
    take:
    ch_filtered_vcf

    main:
    dp_mask_script = file("${projectDir}/vcf_filtering_scripts/sample_dp_mask.py")
    PHASE2_SITE_MEAN_DP_LIMITS_PROCESS(ch_filtered_vcf, dp_mask_script)

    emit:
    dp_limits = PHASE2_SITE_MEAN_DP_LIMITS_PROCESS.out.dp_limits
}
