process GATK_GENOTYPEGVCFS_BY_CHR_PROCESS {
    tag "cohort:${interval_idx}"
    container "${params.container_image}"

    input:
    tuple val(interval_idx), val(interval_spec), path(genomicsdb_dir), path(ref_fa), path(ref_fai), path(ref_dict)

    output:
    tuple val('cohort'), val(interval_idx), val(interval_spec), path("cohort.${interval_idx}.gatk.raw.vcf.gz"), path("cohort.${interval_idx}.gatk.raw.vcf.gz.tbi"), emit: raw_vcf_per_chr

    script:
    """
    gatk GenotypeGVCFs \
      -R ${ref_fa} \
      -V gendb://\$(basename ${genomicsdb_dir}) \
      ${params.gatk_genotypegvcfs_parameters} \
      -O cohort.${interval_idx}.gatk.raw.vcf.gz

    tabix -f -p vcf cohort.${interval_idx}.gatk.raw.vcf.gz
    """
}

workflow GATK_GENOTYPEGVCFS {
    take:
    ch_gendb
    ch_ref

    main:
    // ch_ref is a value channel of tuple(ref_fa, ref_fai, ref_dict) from PREPARE_REFERENCE.
    ch_input = ch_gendb.combine(ch_ref).map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5]) }
    GATK_GENOTYPEGVCFS_BY_CHR_PROCESS(ch_input)

    emit:
    raw_vcf = GATK_GENOTYPEGVCFS_BY_CHR_PROCESS.out.raw_vcf_per_chr
}
