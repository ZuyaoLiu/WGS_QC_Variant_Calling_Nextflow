process MULTIQC_RAW_PROCESS {
    tag 'raw_fastqc'
    container "${params.container_image}"
    publishDir 'results/01_raw_qc/multiqc', mode: 'symlink'

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

workflow MULTIQC_RAW {
    take:
    ch_fastqc_zip

    main:
    MULTIQC_RAW_PROCESS(ch_fastqc_zip.collect())

    emit:
    multiqc_report = MULTIQC_RAW_PROCESS.out.report
    multiqc_data   = MULTIQC_RAW_PROCESS.out.data_dir
}
