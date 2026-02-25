process BCFTOOLS_CALL_BY_CHR_PROCESS {
    tag "${callset_id}:${chrom}"
    container "${params.sif}"
    cpus { (params.bcftools_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/05_variant_calling/bcftools/per_chrom', mode: 'copy'

    input:
    tuple val(callset_id), path(markdup_bams), path(markdup_bais), path(ref_fa), val(interval_idx), val(chrom)

    output:
    tuple val(callset_id), val(interval_idx), val(chrom), path("${callset_id}.${interval_idx}.bcftools.vcf.gz"), path("${callset_id}.${interval_idx}.bcftools.vcf.gz.tbi"), emit: per_chr

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

    bcftools mpileup \
      --threads ${task.cpus} \
      -Ou \
      ${params.bcftools_mpileup_parameters} \
      -f ${ref_fa} \
      -r ${chrom} \
      ${markdup_bams} \
      | bcftools call -mv -Oz --threads ${task.cpus} ${params.bcftools_call_parameters} -o ${callset_id}.${interval_idx}.bcftools.vcf.gz

    tabix -f -p vcf ${callset_id}.${interval_idx}.bcftools.vcf.gz
    """
}

process BCFTOOLS_MERGE_PROCESS {
    tag "${callset_id}"
    container "${params.sif}"
    cpus { (params.bcftools_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/05_variant_calling/bcftools', mode: 'copy'

    input:
    tuple val(callset_id), path(vcf_files), path(vcf_tbis)

    output:
    tuple val(callset_id), path("${callset_id}.bcftools.vcf.gz"), path("${callset_id}.bcftools.vcf.gz.tbi"), emit: merged

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export POLARS_MAX_THREADS=${task.cpus}
    export RAYON_NUM_THREADS=${task.cpus}

    ls -1 ${vcf_files} | sort > vcf.list
    bcftools concat --threads ${task.cpus} -a ${params.bcftools_concat_parameters} -f vcf.list -Oz -o ${callset_id}.bcftools.vcf.gz
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
        .map { sample_id, bam, bai -> tuple('cohort', bam, bai) }
        .groupTuple()

    ch_input = ch_callset.combine(ch_ref).combine(ch_intervals).map { row ->
        tuple(row[0], row[1], row[2], row[3], row[4], row[5])
    }

    BCFTOOLS_CALL_BY_CHR_PROCESS(ch_input)

    ch_for_merge = BCFTOOLS_CALL_BY_CHR_PROCESS.out.per_chr
        .map { callset_id, interval_idx, chrom, vcf, tbi -> tuple(callset_id, vcf, tbi) }
        .groupTuple()

    BCFTOOLS_MERGE_PROCESS(ch_for_merge)

    emit:
    bcf_vcf = BCFTOOLS_MERGE_PROCESS.out.merged
    bcf_per_chr = BCFTOOLS_CALL_BY_CHR_PROCESS.out.per_chr
}
