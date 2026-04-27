process PHASE2_REMOVE_HIGH_MISSING_SAMPLES_PROCESS {
    tag "${callset_id}"
    container "${params.container_image}"
    publishDir 'results/05_vcf_filtering/02_phase2/intermediate', mode: 'copy'

    input:
    tuple val(callset_id), path(input_vcf), path(input_tbi), path(sample_missing_tsv)

    output:
    tuple val(callset_id), path("${callset_id}.phase2.sample_filtered.vcf.gz"), path("${callset_id}.phase2.sample_filtered.vcf.gz.tbi"), emit: filtered_vcf
    tuple val(callset_id), path("${callset_id}.samples_to_remove.txt"), path("${callset_id}.samples_to_remove.reasons.tsv"), emit: remove_lists

    script:
    def userList = params.phase2_high_missing_sample_file ? file(params.phase2_high_missing_sample_file) : null
    def userListArg = userList ? userList.toAbsolutePath().toString() : ''
    def maxMissing = params.phase2_max_sample_missing != null ? params.phase2_max_sample_missing.toString() : ''
    """
    : > ${callset_id}.samples_to_remove.txt
    {
      printf "sample\\treason\\n"
    } > ${callset_id}.samples_to_remove.reasons.tsv

    if [[ "${params.phase2_remove_high_missing_samples}" == "true" ]]; then
      if [[ -n "${maxMissing}" ]]; then
        awk -F'\\t' -v max_missing="${maxMissing}" '
          NR == 1 { next }
          (\$4 + 0) > max_missing { print \$1 }
        ' ${sample_missing_tsv} | sort -u >> ${callset_id}.samples_to_remove.txt

        awk -F'\\t' -v max_missing="${maxMissing}" '
          NR == 1 { next }
          (\$4 + 0) > max_missing { printf "%s\\tmissing_rate_gt_%s\\n", \$1, max_missing }
        ' ${sample_missing_tsv} >> ${callset_id}.samples_to_remove.reasons.tsv
      fi

      if [[ -n "${userListArg}" ]]; then
        awk 'NF > 0 { print \$1 }' ${userListArg} | sort -u >> ${callset_id}.samples_to_remove.txt
        awk 'NF > 0 { printf "%s\\tuser_sample_file\\n", \$1 }' ${userListArg} >> ${callset_id}.samples_to_remove.reasons.tsv
      fi
    fi

    sort -u ${callset_id}.samples_to_remove.txt -o ${callset_id}.samples_to_remove.txt
    awk -F'\\t' 'NR == 1 { print; next } !seen[\$1 FS \$2]++ { print }' ${callset_id}.samples_to_remove.reasons.tsv > ${callset_id}.samples_to_remove.reasons.tsv.tmp
    mv ${callset_id}.samples_to_remove.reasons.tsv.tmp ${callset_id}.samples_to_remove.reasons.tsv

    if [[ -s ${callset_id}.samples_to_remove.txt ]]; then
      bcftools view \
        --threads ${task.cpus} \
        -S ^${callset_id}.samples_to_remove.txt \
        --force-samples \
        -Oz \
        -o ${callset_id}.phase2.sample_filtered.vcf.gz \
        ${input_vcf}
    else
      bcftools view \
        --threads ${task.cpus} \
        -Oz \
        -o ${callset_id}.phase2.sample_filtered.vcf.gz \
        ${input_vcf}
    fi

    tabix -f -p vcf ${callset_id}.phase2.sample_filtered.vcf.gz
    """
}

workflow PHASE2_REMOVE_HIGH_MISSING_SAMPLES {
    take:
    ch_called_vcf
    ch_sample_missingness

    main:
    ch_input = ch_called_vcf.combine(ch_sample_missingness).map { row ->
        tuple(row[0], row[1], row[2], row[4])
    }
    PHASE2_REMOVE_HIGH_MISSING_SAMPLES_PROCESS(ch_input)

    emit:
    filtered_vcf = PHASE2_REMOVE_HIGH_MISSING_SAMPLES_PROCESS.out.filtered_vcf.map { callset_id, vcf, tbi ->
        tuple(
            callset_id,
            file("results/05_vcf_filtering/02_phase2/intermediate/${vcf.name}"),
            file("results/05_vcf_filtering/02_phase2/intermediate/${tbi.name}")
        )
    }
    remove_lists = PHASE2_REMOVE_HIGH_MISSING_SAMPLES_PROCESS.out.remove_lists.map { callset_id, remove_list, reason_tsv ->
        tuple(
            callset_id,
            file("results/05_vcf_filtering/02_phase2/intermediate/${remove_list.name}"),
            file("results/05_vcf_filtering/02_phase2/intermediate/${reason_tsv.name}")
        )
    }
}
