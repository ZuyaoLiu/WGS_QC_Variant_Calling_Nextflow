process GATK_HAPLOTYPECALLER_BY_CHR_PROCESS {
    tag "${sample_id}:${interval_idx}"
    container "${params.container_image}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bai), path(ref_fa), path(ref_fai), path(ref_dict), val(interval_idx), val(interval_spec)

    output:
    tuple val(sample_id), val(interval_idx), val(interval_spec), path("${sample_id}.${interval_idx}.g.vcf.gz"), path("${sample_id}.${interval_idx}.g.vcf.gz.*"), emit: gvcf_per_chr

    script:
    """
    gatk HaplotypeCaller \
      -R ${ref_fa} \
      -I ${input_bam} \
      -O ${sample_id}.${interval_idx}.g.vcf.gz \
      -L ${interval_spec} \
      -ERC GVCF \
      ${params.gatk_haplotypecaller_parameters} \
      --sample-name ${sample_id} \
      --native-pair-hmm-threads ${task.cpus}
    """
}

workflow GATK_HAPLOTYPECALLER {
    take:
    ch_bam
    ch_ref
    ch_intervals

    main:
    // ch_ref is a value channel of tuple(ref_fa, ref_fai, ref_dict) from PREPARE_REFERENCE.
    ch_input = ch_bam.combine(ch_ref).combine(ch_intervals).map { row ->
      tuple(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7])
    }
    GATK_HAPLOTYPECALLER_BY_CHR_PROCESS(ch_input)

    emit:
    gvcf = GATK_HAPLOTYPECALLER_BY_CHR_PROCESS.out.gvcf_per_chr
}
