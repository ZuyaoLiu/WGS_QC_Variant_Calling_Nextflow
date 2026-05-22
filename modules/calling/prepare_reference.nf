process PREPARE_REFERENCE_PROCESS {
    tag "${ref_fa.baseName}"
    container "${params.container_image}"

    input:
    path(ref_fa)

    output:
    tuple path(ref_fa, includeInputs: true), path("${ref_fa}.fai"), path("${ref_fa.baseName}.dict"), emit: ref_indexed

    script:
    """
    if [ ! -f ${ref_fa}.fai ]; then
      samtools faidx ${ref_fa}
    fi

    if [ ! -f ${ref_fa.baseName}.dict ]; then
      gatk CreateSequenceDictionary -R ${ref_fa} -O ${ref_fa.baseName}.dict
    fi
    """
}

workflow PREPARE_REFERENCE {
    take:
    ch_ref

    main:
    PREPARE_REFERENCE_PROCESS(ch_ref)

    emit:
    ref_indexed = PREPARE_REFERENCE_PROCESS.out.ref_indexed
}
