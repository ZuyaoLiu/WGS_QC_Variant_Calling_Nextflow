process BCFTOOLS_CALL_BY_CHR_PROCESS {
    tag "${sample_id}:${chrom}"
    container "${params.sif}"
    cpus params.bcftools_cpus
    publishDir 'results/05_variant_calling/bcftools/per_chrom', mode: 'copy'

    input:
    tuple val(sample_id), path(markdup_bam), path(markdup_bai), path(ref_fa), val(interval_idx), val(chrom)

    output:
    tuple val(sample_id), val(interval_idx), val(chrom), path("${sample_id}.${interval_idx}.bcftools.vcf.gz"), path("${sample_id}.${interval_idx}.bcftools.vcf.gz.tbi"), emit: per_chr

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
      -f ${ref_fa} \
      -r ${chrom} \
      ${markdup_bam} \
      | bcftools call -mv -Oz --threads ${task.cpus} -o ${sample_id}.${interval_idx}.bcftools.vcf.gz

    tabix -f -p vcf ${sample_id}.${interval_idx}.bcftools.vcf.gz
    """
}

process BCFTOOLS_MERGE_PROCESS {
    tag "${sample_id}"
    container "${params.sif}"
    cpus params.bcftools_cpus
    publishDir 'results/05_variant_calling/bcftools', mode: 'copy'

    input:
    tuple val(sample_id), path(vcf_files), path(vcf_tbis)

    output:
    tuple val(sample_id), path("${sample_id}.bcftools.vcf.gz"), path("${sample_id}.bcftools.vcf.gz.tbi"), emit: merged

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export POLARS_MAX_THREADS=${task.cpus}
    export RAYON_NUM_THREADS=${task.cpus}

    ls -1 ${vcf_files} | sort > vcf.list
    bcftools concat --threads ${task.cpus} -a -f vcf.list -Oz -o ${sample_id}.bcftools.vcf.gz
    tabix -f -p vcf ${sample_id}.bcftools.vcf.gz
    """
}

workflow BCFTOOLS_CALL {
    take:
    ch_markdup
    ch_ref
    ch_intervals

    main:
    ch_input = ch_markdup.combine(ch_ref).combine(ch_intervals).map { row ->
        tuple(row[0], row[1], row[2], row[3], row[4], row[5])
    }

    BCFTOOLS_CALL_BY_CHR_PROCESS(ch_input)

    ch_for_merge = BCFTOOLS_CALL_BY_CHR_PROCESS.out.per_chr
        .map { sample_id, interval_idx, chrom, vcf, tbi -> tuple(sample_id, vcf, tbi) }
        .groupTuple()

    BCFTOOLS_MERGE_PROCESS(ch_for_merge)

    emit:
    bcf_vcf = BCFTOOLS_MERGE_PROCESS.out.merged
    bcf_per_chr = BCFTOOLS_CALL_BY_CHR_PROCESS.out.per_chr
}
