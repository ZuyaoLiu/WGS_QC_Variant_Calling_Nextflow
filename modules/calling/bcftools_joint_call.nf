process BCFTOOLS_CALL_BY_CHR_PROCESS {
    tag "${callset_id}:${interval_idx}"
    container "${params.container_image}"

    input:
    tuple val(callset_id), path(markdup_alignments), path(markdup_indexes), path(ref_fa), path(ref_fai), path(ref_dict), val(interval_idx), val(interval_spec)

    output:
    tuple val(callset_id), val(interval_idx), val(interval_spec), path("${callset_id}.${interval_idx}.bcftools.vcf.gz"), path("${callset_id}.${interval_idx}.bcftools.vcf.gz.tbi"), emit: per_chr

    script:
    """
    bcftools mpileup \
      --threads ${task.cpus} \
      -q 20 -a DP,AD \
      -Ou \
      ${params.bcftools_mpileup_parameters} \
      -f ${ref_fa} \
      -r ${interval_spec} \
      ${markdup_alignments} \
      | bcftools call -mv -a GQ,GP -Oz --threads ${task.cpus} ${params.bcftools_call_parameters} -o ${callset_id}.${interval_idx}.bcftools.vcf.gz

    tabix -f -p vcf ${callset_id}.${interval_idx}.bcftools.vcf.gz
    """
}

process BCFTOOLS_MERGE_PROCESS {
    tag "${callset_id}"
    container "${params.container_image}"
    publishDir 'results/04_variant_calling/bcftools', mode: 'symlink'

    input:
    tuple val(callset_id), path(vcf_files), path(vcf_tbis)

    output:
    tuple val(callset_id), path("${callset_id}.bcftools.vcf.gz"), path("${callset_id}.bcftools.vcf.gz.tbi"), emit: merged

    script:
    """
    ls -1 ${vcf_files} | sort > vcf.list
    bcftools concat --threads ${task.cpus} ${params.bcftools_concat_parameters} -f vcf.list -Oz -o ${callset_id}.bcftools.vcf.gz
    tabix -f -p vcf ${callset_id}.bcftools.vcf.gz
    """
}

workflow BCFTOOLS_CALL {
    take:
    ch_markdup
    ch_ref
    ch_intervals

    main:
    ch_callset = ch_markdup
        .map { sample_id, alignment_file, alignment_index -> tuple('cohort', alignment_file, alignment_index) }
        .groupTuple()

    // ch_ref is a value channel of tuple(ref_fa, ref_fai, ref_dict) from PREPARE_REFERENCE.
    ch_input = ch_callset.combine(ch_ref).combine(ch_intervals).map { row ->
        tuple(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7])
    }

    BCFTOOLS_CALL_BY_CHR_PROCESS(ch_input)

    ch_per_chr = BCFTOOLS_CALL_BY_CHR_PROCESS.out.per_chr

    ch_for_merge = ch_per_chr
        .map { callset_id, interval_idx, interval_spec, vcf, tbi -> tuple(callset_id, vcf, tbi) }
        .groupTuple()

    BCFTOOLS_MERGE_PROCESS(ch_for_merge)

    emit:
    bcf_vcf = BCFTOOLS_MERGE_PROCESS.out.merged
    bcf_per_chr = ch_per_chr
}
