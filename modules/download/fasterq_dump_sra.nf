process FASTERQ_DUMP_SE_PROCESS {
    tag "${sample_id}:${accession}"
    container "${params.container_image}"
    publishDir 'results/00_sra_download', mode: 'move'

    input:
    tuple val(sample_id), val(accession)

    output:
    tuple val(sample_id), val('SE'), path("${sample_id}.fastq.gz"), emit: raw_fastq_se
    tuple val(sample_id), path("${sample_id}.fastq.gz"), emit: fastp_input_se

    script:
    """
    mkdir -p tmp_fasterq

    fasterq-dump \
      -O . \
      -t tmp_fasterq \
      -e ${task.cpus} \
      -f \
      ${accession}

    gzip -f ${accession}.fastq
    mv ${accession}.fastq.gz ${sample_id}.fastq.gz
    """
}

process FASTERQ_DUMP_PE_PROCESS {
    tag "${sample_id}:${accession}"
    container "${params.container_image}"
    publishDir 'results/00_sra_download', mode: 'move'

    input:
    tuple val(sample_id), val(accession)

    output:
    tuple val(sample_id), val('R1'), path("${sample_id}_R1.fastq.gz"), emit: raw_fastq_r1
    tuple val(sample_id), val('R2'), path("${sample_id}_R2.fastq.gz"), emit: raw_fastq_r2
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: fastp_input_pe

    script:
    """
    mkdir -p tmp_fasterq

    fasterq-dump \
      --split-files \
      -O . \
      -t tmp_fasterq \
      -e ${task.cpus} \
      -f \
      ${accession}

    gzip -f ${accession}_1.fastq
    gzip -f ${accession}_2.fastq

    mv ${accession}_1.fastq.gz ${sample_id}_R1.fastq.gz
    mv ${accession}_2.fastq.gz ${sample_id}_R2.fastq.gz
    """
}

workflow FASTERQ_DUMP_RAW {
    take:
    ch_sra_samples

    main:
    if (params.read_type == 'SE') {
        FASTERQ_DUMP_SE_PROCESS(ch_sra_samples)
    } else {
        FASTERQ_DUMP_PE_PROCESS(ch_sra_samples)
    }

    emit:
    raw_fastq = params.read_type == 'SE'
        ? FASTERQ_DUMP_SE_PROCESS.out.raw_fastq_se.map { sample_id, read_label, fq ->
            tuple(sample_id, read_label, file("results/00_sra_download/${fq.name}"))
        }
        : FASTERQ_DUMP_PE_PROCESS.out.raw_fastq_r1.map { sample_id, read_label, fq ->
            tuple(sample_id, read_label, file("results/00_sra_download/${fq.name}"))
          }.mix(
            FASTERQ_DUMP_PE_PROCESS.out.raw_fastq_r2.map { sample_id, read_label, fq ->
                tuple(sample_id, read_label, file("results/00_sra_download/${fq.name}"))
            }
          )

    fastp_input = params.read_type == 'SE'
        ? FASTERQ_DUMP_SE_PROCESS.out.fastp_input_se.map { sample_id, fq ->
            tuple(sample_id, file("results/00_sra_download/${fq.name}"))
        }
        : FASTERQ_DUMP_PE_PROCESS.out.fastp_input_pe.map { sample_id, r1, r2 ->
            tuple(
                sample_id,
                file("results/00_sra_download/${r1.name}"),
                file("results/00_sra_download/${r2.name}")
            )
        }
}
