process FASTP_SE_PROCESS {
    tag "${sample_id}:SE"
    container "${params.sif}"
    cpus { (params.fastp_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/02_fastp_qc/fastp', mode: 'copy'

    input:
    tuple val(sample_id), path(fq)

    output:
    tuple val(sample_id), path("${sample_id}.clean.fastq.gz"), emit: cleaned_se
    tuple val(sample_id), val('SE'), path("${sample_id}.clean.fastq.gz"), emit: cleaned_for_qc_se
    path("${sample_id}.fastp.html"), emit: html
    path("${sample_id}.fastp.json"), emit: json

    script:
    """
    fastp \
      -i ${fq} \
      -o ${sample_id}.clean.fastq.gz \
      -h ${sample_id}.fastp.html \
      -j ${sample_id}.fastp.json \
      -w ${task.cpus} \
      ${params.fastp_parameters}
    """
}

process FASTP_PE_PROCESS {
    tag "${sample_id}:PE"
    container "${params.sif}"
    cpus { (params.fastp_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/02_fastp_qc/fastp', mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.clean.R1.fastq.gz"), path("${sample_id}.clean.R2.fastq.gz"), emit: cleaned_pe
    tuple val(sample_id), val('R1'), path("${sample_id}.clean.R1.fastq.gz"), emit: cleaned_for_qc_r1
    tuple val(sample_id), val('R2'), path("${sample_id}.clean.R2.fastq.gz"), emit: cleaned_for_qc_r2
    path("${sample_id}.fastp.html"), emit: html
    path("${sample_id}.fastp.json"), emit: json

    script:
    """
    fastp \
      -i ${r1} \
      -I ${r2} \
      -o ${sample_id}.clean.R1.fastq.gz \
      -O ${sample_id}.clean.R2.fastq.gz \
      -h ${sample_id}.fastp.html \
      -j ${sample_id}.fastp.json \
      -w ${task.cpus} \
      ${params.fastp_parameters}
    """
}

workflow FASTP {
    take:
    ch_input

    main:
    if (params.read_type == 'SE') {
        FASTP_SE_PROCESS(ch_input)
    } else {
        FASTP_PE_PROCESS(ch_input)
    }

    emit:
    cleaned_reads = params.read_type == 'SE' ? FASTP_SE_PROCESS.out.cleaned_se : FASTP_PE_PROCESS.out.cleaned_pe
    cleaned_for_qc = params.read_type == 'SE' ? FASTP_SE_PROCESS.out.cleaned_for_qc_se : FASTP_PE_PROCESS.out.cleaned_for_qc_r1.mix(FASTP_PE_PROCESS.out.cleaned_for_qc_r2)
    fastp_html = params.read_type == 'SE' ? FASTP_SE_PROCESS.out.html : FASTP_PE_PROCESS.out.html
    fastp_json = params.read_type == 'SE' ? FASTP_SE_PROCESS.out.json : FASTP_PE_PROCESS.out.json
}
