process FASTP_SE_PROCESS {
    tag "${sample_id}:SE"
    container "${params.sif}"
    cpus { (params.fastp_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/02_fastp_qc/fastp', mode: 'move'

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
    publishDir 'results/02_fastp_qc/fastp', mode: 'move'

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
    cleaned_reads = params.read_type == 'SE'
        ? FASTP_SE_PROCESS.out.cleaned_se.map { sample_id, fq ->
            tuple(sample_id, file("results/02_fastp_qc/fastp/${fq.name}"))
        }
        : FASTP_PE_PROCESS.out.cleaned_pe.map { sample_id, r1, r2 ->
            tuple(
                sample_id,
                file("results/02_fastp_qc/fastp/${r1.name}"),
                file("results/02_fastp_qc/fastp/${r2.name}")
            )
        }

    cleaned_for_qc = params.read_type == 'SE'
        ? FASTP_SE_PROCESS.out.cleaned_for_qc_se.map { sample_id, read_label, fq ->
            tuple(sample_id, read_label, file("results/02_fastp_qc/fastp/${fq.name}"))
        }
        : FASTP_PE_PROCESS.out.cleaned_for_qc_r1.map { sample_id, read_label, fq ->
            tuple(sample_id, read_label, file("results/02_fastp_qc/fastp/${fq.name}"))
          }.mix(
            FASTP_PE_PROCESS.out.cleaned_for_qc_r2.map { sample_id, read_label, fq ->
                tuple(sample_id, read_label, file("results/02_fastp_qc/fastp/${fq.name}"))
            }
          )

    fastp_html = (params.read_type == 'SE' ? FASTP_SE_PROCESS.out.html : FASTP_PE_PROCESS.out.html)
        .map { html_file -> file("results/02_fastp_qc/fastp/${html_file.name}") }
    fastp_json = (params.read_type == 'SE' ? FASTP_SE_PROCESS.out.json : FASTP_PE_PROCESS.out.json)
        .map { json_file -> file("results/02_fastp_qc/fastp/${json_file.name}") }
}
