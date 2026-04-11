process BWA_INDEX_PROCESS {
    tag "bwa_index:${ref_fa.baseName}"
    container "${params.sif}"
    cpus { (params.bwa_cpus ?: params.threads ?: 1) as Integer }

    input:
    path ref_fa

    output:
    path("${ref_fa.name}*"), emit: ref_bundle

    script:
    """
    bwa index ${ref_fa}
    """
}

workflow BWA_INDEX {
    take:
    ch_ref

    main:
    BWA_INDEX_PROCESS(ch_ref)

    emit:
    ref_bundle = BWA_INDEX_PROCESS.out.ref_bundle.collect()
}
