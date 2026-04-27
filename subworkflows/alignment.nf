include { ALIGN_MARKDUP_SE } from '../modules/alignment/align_markdup_se'
include { ALIGN_MARKDUP_PE } from '../modules/alignment/align_markdup_pe'
include { BWAMEM2_INDEX } from '../modules/alignment/bwamem2_index'

// Alignment generates a shared bwa-mem2 index bundle once, then routes reads
// to the SE or PE align-and-markdup implementation.

workflow ALIGNMENT {
    take:
    ch_clean_reads
    ch_ref

    main:
    ref_index = BWAMEM2_INDEX(ch_ref)

    if (params.read_type == 'SE') {
        alignment = ALIGN_MARKDUP_SE(ch_clean_reads, ref_index.ref_bundle, ch_ref)
    } else {
        alignment = ALIGN_MARKDUP_PE(ch_clean_reads, ref_index.ref_bundle, ch_ref)
    }

    emit:
    markdup_alignment = alignment.markdup_alignment
}
