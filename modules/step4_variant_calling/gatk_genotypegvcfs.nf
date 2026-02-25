process GATK_GENOTYPEGVCFS_BY_CHR_PROCESS {
    tag "cohort:${chrom}"
    container "${params.sif}"
    cpus { (params.gatk_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/05_variant_calling/gatk/genotypegvcfs/per_chrom', mode: 'copy'

    input:
    tuple val(interval_idx), val(chrom), path(genomicsdb_dir), path(ref_fa)

    output:
    tuple val('cohort'), val(interval_idx), val(chrom), path("cohort.${interval_idx}.gatk.raw.vcf.gz"), path("cohort.${interval_idx}.gatk.raw.vcf.gz.tbi"), emit: raw_vcf_per_chr

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export POLARS_MAX_THREADS=${task.cpus}
    export RAYON_NUM_THREADS=${task.cpus}

    if [ ! -f ${ref_fa}.fai ]; then
      samtools faidx ${ref_fa}
    fi

    if [ ! -f ${ref_fa.baseName}.dict ]; then
      gatk CreateSequenceDictionary -R ${ref_fa} -O ${ref_fa.baseName}.dict
    fi

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
    ch_input = ch_gendb.combine(ch_ref).map { row -> tuple(row[0], row[1], row[2], row[3]) }
    GATK_GENOTYPEGVCFS_BY_CHR_PROCESS(ch_input)

    emit:
    raw_vcf = GATK_GENOTYPEGVCFS_BY_CHR_PROCESS.out.raw_vcf_per_chr
}
