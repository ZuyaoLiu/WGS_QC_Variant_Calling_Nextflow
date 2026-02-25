process FASTQC_RAW_PROCESS {
    tag "${sample_id}:${read_label}"
    container "${params.sif}"
    cpus { (params.fastqc_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/01_raw_qc/fastqc', mode: 'copy'

    input:
    tuple val(sample_id), val(read_label), path(fq)

    output:
    tuple val(sample_id), val(read_label), path("*_fastqc.zip"), path("*_fastqc.html"), emit: qc

    script:
    """
    fastqc -t ${task.cpus} ${fq}
    """
}

workflow FASTQC_RAW {
    take:
    ch_fastq

    main:
    FASTQC_RAW_PROCESS(ch_fastq)

    emit:
    fastqc_reports = FASTQC_RAW_PROCESS.out.qc
}
