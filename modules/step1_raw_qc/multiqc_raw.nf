process MULTIQC_RAW_PROCESS {
    tag 'raw_fastqc'
    container "${params.sif}"
    cpus { (params.threads ?: 1) as Integer }
    publishDir 'results/01_raw_qc/multiqc', mode: 'copy'

    input:
    path fastqc_zip_files

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data', emit: data_dir

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export POLARS_MAX_THREADS=${task.cpus}
    export RAYON_NUM_THREADS=${task.cpus}
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
