process GATK_HAPLOTYPECALLER_BY_CHR_PROCESS {
    tag "${sample_id}:${chrom}"
    container "${params.sif}"
    cpus params.gatk_cpus
    publishDir 'results/05_variant_calling/gatk/haplotypecaller/per_chrom', mode: 'copy'

    input:
    tuple val(sample_id), path(input_bam), path(input_bai), path(ref_fa), val(interval_idx), val(chrom)

    output:
    tuple val(sample_id), val(interval_idx), val(chrom), path("${sample_id}.${interval_idx}.g.vcf.gz"), path("${sample_id}.${interval_idx}.g.vcf.gz.*"), emit: gvcf_per_chr

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

    gatk HaplotypeCaller \
      -R ${ref_fa} \
      -I ${input_bam} \
      -O ${sample_id}.${interval_idx}.g.vcf.gz \
      -L ${chrom} \
      -ERC GVCF \
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
    ch_input = ch_bam.combine(ch_ref).combine(ch_intervals).map { row ->
      tuple(row[0], row[1], row[2], row[3], row[4], row[5])
    }
    GATK_HAPLOTYPECALLER_BY_CHR_PROCESS(ch_input)

    emit:
    gvcf = GATK_HAPLOTYPECALLER_BY_CHR_PROCESS.out.gvcf_per_chr
}
