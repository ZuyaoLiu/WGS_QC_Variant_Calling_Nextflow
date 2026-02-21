process BWAMEM2_INDEX_PROCESS {
    tag "bwamem2_index:${ref_fa.baseName}"
    container "${params.sif}"
    cpus params.bwamem2_cpus
    publishDir 'results/03_align/index', mode: 'copy'

    input:
    path ref_fa

    output:
    path("${ref_fa.name}*"), emit: ref_bundle

    script:
    """
    bwa-mem2 index ${ref_fa}
    """
}

workflow BWAMEM2_INDEX {
    take:
    ch_ref

    main:
    BWAMEM2_INDEX_PROCESS(ch_ref)

    emit:
    ref_bundle = BWAMEM2_INDEX_PROCESS.out.ref_bundle.collect()
}
