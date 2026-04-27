include { FASTQC_RAW } from '../modules/qc/raw_fastqc'
include { MULTIQC_RAW } from '../modules/qc/raw_multiqc'
include { FASTERQ_DUMP_RAW } from '../modules/download/fasterq_dump_sra'

// Raw QC is split into two entry helpers:
// 1. optional SRA download
// 2. FastQC + MultiQC over raw FASTQ inputs

workflow DOWNLOAD_SRA {
    take:
    ch_sra_samples

    main:
    FASTERQ_DUMP_RAW(ch_sra_samples)

    emit:
    raw_fastq = FASTERQ_DUMP_RAW.out.raw_fastq
    fastp_input = FASTERQ_DUMP_RAW.out.fastp_input
}

workflow RAW_QC {
    take:
    ch_raw_fastq

    main:
    raw_qc = FASTQC_RAW(ch_raw_fastq)
    ch_fastqc_zip = raw_qc.fastqc_reports.map { sample, read_label, zip_file, html_file -> zip_file }
    raw_mqc = MULTIQC_RAW(ch_fastqc_zip)

    emit:
    fastqc_reports = raw_qc.fastqc_reports
    multiqc_report = raw_mqc.multiqc_report
}
