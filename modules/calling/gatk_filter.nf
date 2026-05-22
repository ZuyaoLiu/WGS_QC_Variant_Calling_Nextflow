process GATK_FILTER_PROCESS {
    tag "${callset_id}"
    container "${params.container_image}"
    publishDir 'results/04_variant_calling/gatk/filter', mode: 'symlink'

    input:
    tuple val(callset_id), path(raw_vcf), path(raw_tbi), path(ref_fa), path(ref_fai), path(ref_dict)

    output:
    tuple val(callset_id), path("${callset_id}.gatk.filtered.vcf.gz"), path("${callset_id}.gatk.filtered.vcf.gz.tbi"), emit: filtered_vcf

    script:
    """
    gatk SelectVariants -R ${ref_fa} -V ${raw_vcf} --select-type-to-include SNP -O ${callset_id}.snp.raw.vcf.gz
    gatk SelectVariants -R ${ref_fa} -V ${raw_vcf} --select-type-to-include INDEL -O ${callset_id}.indel.raw.vcf.gz

    gatk VariantFiltration \
      -R ${ref_fa} \
      -V ${callset_id}.snp.raw.vcf.gz \
      --filter-name "SNP_HARD_FILTER" \
      --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
      ${params.gatk_variantfiltration_snp_parameters} \
      -O ${callset_id}.snp.filtered.tagged.vcf.gz

    gatk VariantFiltration \
      -R ${ref_fa} \
      -V ${callset_id}.indel.raw.vcf.gz \
      --filter-name "INDEL_HARD_FILTER" \
      --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0" \
      ${params.gatk_variantfiltration_indel_parameters} \
      -O ${callset_id}.indel.filtered.tagged.vcf.gz

    gatk SelectVariants -R ${ref_fa} -V ${callset_id}.snp.filtered.tagged.vcf.gz --exclude-filtered true -O ${callset_id}.snp.pass.vcf.gz
    gatk SelectVariants -R ${ref_fa} -V ${callset_id}.indel.filtered.tagged.vcf.gz --exclude-filtered true -O ${callset_id}.indel.pass.vcf.gz

    gatk MergeVcfs \
      -I ${callset_id}.snp.pass.vcf.gz \
      -I ${callset_id}.indel.pass.vcf.gz \
      -O ${callset_id}.gatk.filtered.vcf.gz

    tabix -f -p vcf ${callset_id}.gatk.filtered.vcf.gz
    """
}

workflow GATK_FILTER {
    take:
    ch_raw_vcf
    ch_ref

    main:
    // ch_ref is a value channel of tuple(ref_fa, ref_fai, ref_dict) from PREPARE_REFERENCE.
    ch_input = ch_raw_vcf.combine(ch_ref).map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5]) }
    GATK_FILTER_PROCESS(ch_input)

    emit:
    filtered_vcf = GATK_FILTER_PROCESS.out.filtered_vcf
}
