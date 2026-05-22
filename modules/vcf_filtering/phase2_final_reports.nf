process PHASE2_FINAL_REPORTS_PROCESS {
    tag "${callset_id}"
    container "${params.container_image}"
    publishDir 'results/05_vcf_filtering/final', mode: 'copy'

    input:
    tuple val(callset_id), path(final_vcf), path(final_tbi)
    path(compute_script)
    path(plot_script)

    output:
    tuple val(callset_id), path("${callset_id}.final.sample_missingness.tsv"), path("${callset_id}.final.sample_missingness.histogram.pdf"), path("${callset_id}.final.sample_missingness.density.pdf"), path("${callset_id}.final.sample_missingness.cutoff_vs_samples_kept.pdf"), path("${callset_id}.final.bcftools.stats.txt"), emit: reports

    script:
    """
    python ${compute_script} \
      --input-vcf ${final_vcf} \
      --output-tsv ${callset_id}.final.sample_missingness.tsv \
      --threads ${task.cpus}

    Rscript ${plot_script} \
      --input-tsv ${callset_id}.final.sample_missingness.tsv \
      --out-prefix ${callset_id}.final.sample_missingness \
      --cutoffs-pct ${params.phase2_sample_missing_plot_cutoffs_pct} \
      --title-prefix "Final sample missingness"

    bcftools stats --threads ${task.cpus} ${final_vcf} > ${callset_id}.final.bcftools.stats.txt
    """
}

workflow PHASE2_FINAL_REPORTS {
    take:
    ch_final_vcf

    main:
    compute_script = file("${projectDir}/vcf_filtering_scripts/compute_sample_missingness.py")
    plot_script = file("${projectDir}/vcf_filtering_scripts/plot_sample_missingness.R")
    PHASE2_FINAL_REPORTS_PROCESS(ch_final_vcf, compute_script, plot_script)

    emit:
    reports = PHASE2_FINAL_REPORTS_PROCESS.out.reports
}
