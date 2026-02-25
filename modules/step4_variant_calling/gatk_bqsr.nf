process GATK_BQSR_PROCESS {
    tag "${sample_id}"
    container "${params.sif}"
    cpus { (params.gatk_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/05_variant_calling/gatk/bqsr', mode: 'copy'

    input:
    tuple val(sample_id), path(markdup_bam), path(markdup_bai), path(seed_vcf), path(seed_tbi), path(ref_fa)

    output:
    tuple val(sample_id), path("${sample_id}.recal.bam"), path("${sample_id}.recal.bam.bai"), emit: recal_bam
    tuple val(sample_id), path("${sample_id}.highconf.snp.vcf.gz"), path("${sample_id}.highconf.snp.vcf.gz.tbi"), emit: highconf_snp

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

    bcftools view -v snps -i 'QUAL>=40 && DP>=10 && DP<=200' -Oz -o ${sample_id}.highconf.snp.vcf.gz ${seed_vcf}
    tabix -f -p vcf ${sample_id}.highconf.snp.vcf.gz

    gatk BaseRecalibrator \
      --java-options "-Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=false -Dsamjdk.use_async_io_write_tribble=false -Xmx2g" \
      -R ${ref_fa} \
      -I ${markdup_bam} \
      --known-sites ${sample_id}.highconf.snp.vcf.gz \
      ${params.gatk_baserecalibrator_parameters} \
      -O ${sample_id}.recal.table

    gatk ApplyBQSR \
      --java-options "-Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=false -Dsamjdk.use_async_io_write_tribble=false -Xmx2g" \
      -R ${ref_fa} \
      -I ${markdup_bam} \
      --bqsr-recal-file ${sample_id}.recal.table \
      ${params.gatk_applybqsr_parameters} \
      -O ${sample_id}.recal.bam

    samtools index -@ ${task.cpus} ${sample_id}.recal.bam
    """
}

workflow GATK_BQSR {
    take:
    ch_markdup
    ch_seed_vcf
    ch_ref

    main:
    // Apply the same cohort seed VCF to each sample BAM for per-sample recalibration.
    ch_joined = ch_markdup.combine(ch_seed_vcf).map { row ->
        tuple(row[0], row[1], row[2], row[4], row[5])
    }
    ch_with_ref = ch_joined.combine(ch_ref).map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5]) }
    GATK_BQSR_PROCESS(ch_with_ref)

    emit:
    recal_bam = GATK_BQSR_PROCESS.out.recal_bam
    highconf_snp = GATK_BQSR_PROCESS.out.highconf_snp
}
