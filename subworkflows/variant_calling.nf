include { BCFTOOLS_CALL } from '../modules/calling/bcftools_joint_call'
include { GATK_BQSR } from '../modules/calling/gatk_bqsr'
include { GATK_CLEANUP_GVCF } from '../modules/calling/gatk_cleanup_gvcf'
include { GATK_HAPLOTYPECALLER } from '../modules/calling/gatk_haplotypecaller'
include { GATK_GENOMICSDBIMPORT } from '../modules/calling/gatk_genomicsdbimport'
include { GATK_GENOTYPEGVCFS } from '../modules/calling/gatk_genotypegvcfs'
include { GATK_MERGE_RAW } from '../modules/calling/gatk_merge_raw'
include { GATK_FILTER } from '../modules/calling/gatk_filter'
include { PREPARE_REFERENCE } from '../modules/calling/prepare_reference'

// Variant calling groups the complete bcftools and GATK calling branches behind
// one stable subworkflow interface.

workflow VARIANT_CALLING {
    take:
    ch_markdup_alignment
    ch_ref
    ch_intervals

    main:
    ch_ref_indexed = PREPARE_REFERENCE(ch_ref).ref_indexed

    if (params.caller == 'bcftools') {
        bcftools_out = BCFTOOLS_CALL(ch_markdup_alignment, ch_ref_indexed, ch_intervals)
        final_vcf = bcftools_out.bcf_vcf
    } else {
        hc_input = ch_markdup_alignment

        if (params.use_bqsr) {
            ch_bqsr_panel = Channel.value(file(params.bqsr_panel_abs))
            recalibrated_alignment = GATK_BQSR(ch_markdup_alignment, ch_bqsr_panel, ch_ref_indexed)
            hc_input = recalibrated_alignment.recal_alignment
        }

        hc = GATK_HAPLOTYPECALLER(hc_input, ch_ref_indexed, ch_intervals)
        gdb = GATK_GENOMICSDBIMPORT(hc.gvcf, ch_ref_indexed)

        if (params.cleanup_GVCFS) {
            GATK_CLEANUP_GVCF(hc.gvcf, gdb.gendb)
        }

        gt = GATK_GENOTYPEGVCFS(gdb.gendb, ch_ref_indexed)
        merged_raw = GATK_MERGE_RAW(gt.raw_vcf)
        gatk_filtered = GATK_FILTER(merged_raw.raw_merged, ch_ref_indexed)
        final_vcf = gatk_filtered.filtered_vcf
    }

    emit:
    final_vcf = final_vcf
}
