include { PHASE1_INITIAL_FILTER_MASK } from '../modules/vcf_filtering/phase1_initial_filter_mask'
include { PHASE1_SAMPLE_DP_MASK } from '../modules/vcf_filtering/phase1_sample_dp_mask'
include { PHASE1_REPORTS } from '../modules/vcf_filtering/phase1_reports'
include { PHASE2_REMOVE_HIGH_MISSING_SAMPLES } from '../modules/vcf_filtering/phase2_remove_high_missing_samples'
include { PHASE2_SITE_MEAN_DP_LIMITS } from '../modules/vcf_filtering/phase2_site_mean_dp_limits'
include { PHASE2_FINAL_FILTER } from '../modules/vcf_filtering/phase2_final_filter'
include { PHASE2_FINAL_REPORTS } from '../modules/vcf_filtering/phase2_final_reports'

workflow VCF_FILTERING {
    take:
    ch_called_vcf

    main:
    log.info "Step5 VCF_Filtering active: current implementation provides phase1 and phase2"
    log.info "Step5 VCF_Filtering scope: supports biallelic SNP only"

    phase1_initial = PHASE1_INITIAL_FILTER_MASK(ch_called_vcf)
    phase1_masked = PHASE1_SAMPLE_DP_MASK(phase1_initial.filtered_vcf)
    phase1_reports = PHASE1_REPORTS(phase1_masked.masked_vcf, phase1_masked.sample_table)
    def final_vcf_channel

    if (params.vcffilter_stop_after == 'phase1') {
        final_vcf_channel = phase1_masked.masked_vcf
    } else {
        sample_filtered = PHASE2_REMOVE_HIGH_MISSING_SAMPLES(phase1_masked.masked_vcf, phase1_reports.sample_missingness)
        dp_limits = PHASE2_SITE_MEAN_DP_LIMITS(sample_filtered.filtered_vcf)
        final_filter = PHASE2_FINAL_FILTER(sample_filtered.filtered_vcf, dp_limits.dp_limits)
        PHASE2_FINAL_REPORTS(final_filter.final_vcf)
        final_vcf_channel = final_filter.final_vcf
    }

    emit:
    final_vcf = final_vcf_channel
}
