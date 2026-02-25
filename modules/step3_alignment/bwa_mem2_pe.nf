process BWAMEM2_PE_PROCESS {
    tag "${sample_id}"
    container "${params.sif}"
    cpus { (params.bwamem2_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/03_align', mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path(ref_bundle)

    output:
    tuple val(sample_id), path("${sample_id}.unsorted.bam"), emit: bam

    script:
    """
    bwa-mem2 mem -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" ${params.bwamem2_parameters} ${params.ref_base} ${r1} ${r2} | samtools view -@ ${task.cpus} -Sb -o ${sample_id}.unsorted.bam -
    """
}

workflow BWAMEM2_PE {
    take:
    ch_input
    ch_ref_bundle

    main:
    BWAMEM2_PE_PROCESS(ch_input, ch_ref_bundle)

    emit:
    bam_for_fixmate = BWAMEM2_PE_PROCESS.out.bam
}
