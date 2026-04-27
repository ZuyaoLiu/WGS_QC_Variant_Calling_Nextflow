process PHASE1_REPORTS_PROCESS {
    tag "${callset_id}"
    container "${params.container_image}"
    publishDir 'results/05_vcf_filtering/01_phase1/review', mode: 'copy'

    input:
    tuple val(callset_id), path(masked_bcf), path(masked_csi), path(sample_mean_dp_tsv)
    path(compute_script)
    path(plot_missing_script)
    path(plot_dp_script)

    output:
    tuple val(callset_id), path("${callset_id}.phase1.sample_missingness.tsv"), emit: sample_missingness
    tuple val(callset_id), path("${callset_id}.phase1.sample_missingness.tsv"), path("${callset_id}.phase1.sample_missingness.histogram.pdf"), path("${callset_id}.phase1.sample_missingness.density.pdf"), path("${callset_id}.phase1.sample_missingness.cutoff_vs_samples_kept.pdf"), path("${callset_id}.phase1.sample_mean_dp.tsv"), path("${callset_id}.phase1.sample_mean_dp.histogram.pdf"), path("${callset_id}.phase1.sample_mean_dp.density.pdf"), path("${callset_id}.phase1.bcftools.stats.txt"), emit: reports

    script:
    """

    python ${compute_script} \
      --input-vcf ${masked_bcf} \
      --output-tsv ${callset_id}.phase1.sample_missingness.tsv \
      --threads ${task.cpus}

    Rscript ${plot_missing_script} \
      --input-tsv ${callset_id}.phase1.sample_missingness.tsv \
      --out-prefix ${callset_id}.phase1.sample_missingness \
      --cutoffs-pct ${params.phase1_sample_missing_plot_cutoffs_pct} \
      --title-prefix "Phase1 sample missingness"

    Rscript ${plot_dp_script} \
      --input-tsv ${callset_id}.phase1.sample_mean_dp.tsv \
      --out-prefix ${callset_id}.phase1.sample_mean_dp \
      --dp-column mean_dp \
      --title-prefix "Phase1 sample mean DP"

    bcftools stats --threads ${task.cpus} ${masked_bcf} > ${callset_id}.phase1.bcftools.stats.txt
    """
}

workflow PHASE1_REPORTS {
    take:
    ch_masked_vcf
    ch_sample_table

    main:
    compute_script = file("${projectDir}/vcf_filtering_scripts/compute_sample_missingness.py")
    plot_missing_script = file("${projectDir}/vcf_filtering_scripts/plot_sample_missingness.R")
    plot_dp_script = file("${projectDir}/vcf_filtering_scripts/plot_sample_dp_distribution.R")
    ch_input = ch_masked_vcf.combine(ch_sample_table).map { row ->
        tuple(row[0], row[1], row[2], row[4])
    }
    PHASE1_REPORTS_PROCESS(ch_input, compute_script, plot_missing_script, plot_dp_script)

    emit:
    sample_missingness = PHASE1_REPORTS_PROCESS.out.sample_missingness.map { callset_id, missing_tsv ->
        tuple(
            callset_id,
            file("results/05_vcf_filtering/01_phase1/review/${missing_tsv.name}")
        )
    }
    reports = PHASE1_REPORTS_PROCESS.out.reports.map { callset_id, missing_tsv, hist_pdf, density_pdf, cutoff_pdf, sample_dp_tsv, dp_hist_pdf, dp_density_pdf, stats_txt ->
        tuple(
            callset_id,
            file("results/05_vcf_filtering/01_phase1/review/${missing_tsv.name}"),
            file("results/05_vcf_filtering/01_phase1/review/${hist_pdf.name}"),
            file("results/05_vcf_filtering/01_phase1/review/${density_pdf.name}"),
            file("results/05_vcf_filtering/01_phase1/review/${cutoff_pdf.name}"),
            file("results/05_vcf_filtering/01_phase1/review/${sample_dp_tsv.name}"),
            file("results/05_vcf_filtering/01_phase1/review/${dp_hist_pdf.name}"),
            file("results/05_vcf_filtering/01_phase1/review/${dp_density_pdf.name}"),
            file("results/05_vcf_filtering/01_phase1/review/${stats_txt.name}")
        )
    }
}
