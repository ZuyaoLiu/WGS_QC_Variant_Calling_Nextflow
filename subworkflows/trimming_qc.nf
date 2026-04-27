include { FASTP } from '../modules/qc/trim_fastp'
include { FASTQC_POST } from '../modules/qc/trimmed_fastqc'
include { MULTIQC_POST } from '../modules/qc/trimmed_multiqc'

// Trimming QC wraps the fastp run plus post-trim FastQC/MultiQC reporting.

workflow TRIMMING_QC {
    take:
    ch_fastp_input

    main:
    trim = FASTP(ch_fastp_input)
    post_qc = FASTQC_POST(trim.cleaned_for_qc)
    ch_post_fastqc_zip = post_qc.fastqc_reports.map { sample, read_label, zip_file, html_file -> zip_file }
    post_mqc = MULTIQC_POST(ch_post_fastqc_zip)

    emit:
    cleaned_reads = trim.cleaned_reads
    cleaned_for_qc = trim.cleaned_for_qc
    multiqc_report = post_mqc.multiqc_report
}
