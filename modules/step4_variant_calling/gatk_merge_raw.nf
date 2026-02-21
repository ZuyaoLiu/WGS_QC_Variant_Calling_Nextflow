process GATK_MERGE_RAW_PROCESS {
    tag "${sample_id}"
    container "${params.sif}"
    cpus params.gatk_cpus
    publishDir 'results/05_variant_calling/gatk/genotypegvcfs', mode: 'copy'

    input:
    tuple val(sample_id), path(raw_vcf_files), path(raw_vcf_tbis)

    output:
    tuple val(sample_id), path("${sample_id}.gatk.raw.vcf.gz"), path("${sample_id}.gatk.raw.vcf.gz.tbi"), emit: raw_merged

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export POLARS_MAX_THREADS=${task.cpus}
    export RAYON_NUM_THREADS=${task.cpus}

    ls -1 ${raw_vcf_files} | sort > raw_vcf.list
    bcftools concat --threads ${task.cpus} -a -f raw_vcf.list -Oz -o ${sample_id}.gatk.raw.vcf.gz
    tabix -f -p vcf ${sample_id}.gatk.raw.vcf.gz
    """
}

workflow GATK_MERGE_RAW {
    take:
    ch_raw_vcf_per_chr

    main:
    ch_for_merge = ch_raw_vcf_per_chr
      .map { sample_id, interval_idx, chrom, raw_vcf, raw_tbi -> tuple(sample_id, raw_vcf, raw_tbi) }
      .groupTuple()

    GATK_MERGE_RAW_PROCESS(ch_for_merge)

    emit:
    raw_merged = GATK_MERGE_RAW_PROCESS.out.raw_merged
}
