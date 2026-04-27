process PHASE2_FINAL_FILTER_PROCESS {
    tag "${callset_id}"
    container "${params.container_image}"
    publishDir 'results/05_vcf_filtering/final', mode: 'copy'

    input:
    tuple val(callset_id), path(input_vcf), path(input_tbi), path(sample_mean_dp_tsv), path(site_dp_limits_tsv)

    output:
    tuple val(callset_id), path("${callset_id}.vcffilter.phase2.vcf.gz"), path("${callset_id}.vcffilter.phase2.vcf.gz.tbi"), emit: final_vcf

    script:
    """
    min_site_mean_dp=\$(awk -F'\\t' '\$1 == "min_site_mean_dp" { print \$2 }' ${site_dp_limits_tsv})
    max_site_mean_dp=\$(awk -F'\\t' '\$1 == "max_site_mean_dp" { print \$2 }' ${site_dp_limits_tsv})

    site_filter_expr="QUAL>=${params.phase2_min_qual} && F_MISSING<=${params.phase2_max_site_missing}"

    if [[ "${params.phase2_enable_site_mean_dp_filter}" == "true" && "\${min_site_mean_dp}" != "NA" && "\${max_site_mean_dp}" != "NA" ]]; then
      site_filter_expr="\${site_filter_expr} && MEAN(FMT/DP)>=\${min_site_mean_dp} && MEAN(FMT/DP)<=\${max_site_mean_dp}"
    fi

    bcftools view \
      --threads ${task.cpus} \
      -m2 -M2 -v snps \
      -Ou \
      ${input_vcf} \
    | bcftools +fill-tags \
      --threads ${task.cpus} \
      -Ou \
      -- -t F_MISSING \
    | bcftools view \
      --threads ${task.cpus} \
      -i "\${site_filter_expr}" \
      -Oz \
      -o ${callset_id}.vcffilter.phase2.vcf.gz

    tabix -f -p vcf ${callset_id}.vcffilter.phase2.vcf.gz
    """
}

workflow PHASE2_FINAL_FILTER {
    take:
    ch_filtered_vcf
    ch_dp_limits

    main:
    ch_input = ch_filtered_vcf.combine(ch_dp_limits).map { row ->
        tuple(row[0], row[1], row[2], row[4], row[5])
    }
    PHASE2_FINAL_FILTER_PROCESS(ch_input)

    emit:
    final_vcf = PHASE2_FINAL_FILTER_PROCESS.out.final_vcf.map { callset_id, vcf, tbi ->
        tuple(
            callset_id,
            file("results/05_vcf_filtering/final/${vcf.name}"),
            file("results/05_vcf_filtering/final/${tbi.name}")
        )
    }
}
