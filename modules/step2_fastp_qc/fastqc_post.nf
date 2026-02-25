process FASTQC_POST_PROCESS {
    tag "${sample_id}:${read_label}"
    container "${params.sif}"
    cpus { (params.fastqc_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/02_fastp_qc/fastqc', mode: 'copy'

    input:
    tuple val(sample_id), val(read_label), path(fq)

    output:
    tuple val(sample_id), val(read_label), path('*_fastqc.zip'), path('*_fastqc.html'), emit: qc

    script:
    """
    fastqc -t ${task.cpus} ${fq}
    """
}

workflow FASTQC_POST {
    take:
    ch_fastq

    main:
    FASTQC_POST_PROCESS(ch_fastq)

    emit:
    fastqc_reports = FASTQC_POST_PROCESS.out.qc
}
