process MULTIQC_POST_PROCESS {
    tag 'post_fastqc'
    container "${params.container_image}"
    publishDir 'results/02_fastp_qc/multiqc', mode: 'symlink'

    input:
    path fastqc_zip_files

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data', emit: data_dir

    script:
    """
    multiqc -o . ${fastqc_zip_files}
    """
}

workflow MULTIQC_POST {
    take:
    ch_fastqc_zip

    main:
    MULTIQC_POST_PROCESS(ch_fastqc_zip.collect())

    emit:
    multiqc_report = MULTIQC_POST_PROCESS.out.report
    multiqc_data   = MULTIQC_POST_PROCESS.out.data_dir
}
