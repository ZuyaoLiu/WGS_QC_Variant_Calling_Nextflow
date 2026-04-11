nextflow.enable.dsl = 2

params.help             = params.help ?: false
params.input_dir        = params.input_dir ?: null
params.ref              = params.ref ?: null
params.read_type        = params.read_type ?: 'PE'
params.caller           = params.caller ?: 'bcftools'
params.use_bqsr         = params.use_bqsr ?: false
params.bqsr_panel       = params.containsKey('bqsr_panel') ? params.get('bqsr_panel') : null
params.run_step         = params.run_step ?: 'all'
params.sif              = params.sif ?: 'WGS_Variant_Calling.sif'
params.fastp_parameters = params.fastp_parameters ?: ''
params.bwa_parameters = params.bwa_parameters ?: ''
params.bwamem2_parameters = params.bwamem2_parameters ?: ''
params.bcftools_mpileup_parameters = params.bcftools_mpileup_parameters ?: ''
params.bcftools_call_parameters = params.bcftools_call_parameters ?: ''
params.bcftools_concat_parameters = params.bcftools_concat_parameters ?: ''
params.gatk_haplotypecaller_parameters = params.gatk_haplotypecaller_parameters ?: ''
params.gatk_genomicsdbimport_parameters = params.gatk_genomicsdbimport_parameters ?: ''
params.gatk_genotypegvcfs_parameters = params.gatk_genotypegvcfs_parameters ?: ''
params.gatk_baserecalibrator_parameters = params.gatk_baserecalibrator_parameters ?: ''
params.gatk_applybqsr_parameters = params.gatk_applybqsr_parameters ?: ''
params.gatk_variantfiltration_snp_parameters = params.gatk_variantfiltration_snp_parameters ?: ''
params.gatk_variantfiltration_indel_parameters = params.gatk_variantfiltration_indel_parameters ?: ''
params.ref_abs          = params.ref ? file(params.ref).toAbsolutePath().toString() : null
params.ref_base         = params.ref_abs ? file(params.ref_abs).getName() : null
params.bqsr_panel_abs   = params.bqsr_panel ? file(params.bqsr_panel).toAbsolutePath().toString() : null


include { FASTQC_RAW }  from './modules/step1_raw_qc/fastqc_raw'
include { MULTIQC_RAW } from './modules/step1_raw_qc/multiqc_raw'

include { FASTP }        from './modules/step2_fastp_qc/fastp'
include { FASTQC_POST }  from './modules/step2_fastp_qc/fastqc_post'
include { MULTIQC_POST } from './modules/step2_fastp_qc/multiqc_post'

include { ALIGN_MARKDUP_SE }  from './modules/step3_alignment/align_markdup_se'
include { ALIGN_MARKDUP_PE }  from './modules/step3_alignment/align_markdup_pe'
include { BWA_INDEX }         from './modules/step3_alignment/bwa_index'
include { BWAMEM2_INDEX }     from './modules/step3_alignment/bwa_mem2_index'

include { BCFTOOLS_CALL }         from './modules/step4_variant_calling/bcftools_call'
include { GATK_BQSR }             from './modules/step4_variant_calling/gatk_bqsr'
include { GATK_HAPLOTYPECALLER }  from './modules/step4_variant_calling/gatk_haplotypecaller'
include { GATK_GENOMICSDBIMPORT } from './modules/step4_variant_calling/gatk_genomicsdbimport'
include { GATK_GENOTYPEGVCFS }    from './modules/step4_variant_calling/gatk_genotypegvcfs'
include { GATK_MERGE_RAW }        from './modules/step4_variant_calling/gatk_merge_raw'
include { GATK_FILTER }           from './modules/step4_variant_calling/gatk_filter'

def print_help() {
    log.info """
NGS Variant Calling Nextflow Pipeline (DSL2, Nextflow 25.10.4)

Execution modes:
  1) Full pipeline (default): --run_step all -> Raw_QC -> Trimming_QC -> Aligning -> Calling
  2) Start from a specific step and continue to the end:
     --run_step Raw_QC       -> Raw_QC -> Trimming_QC -> Aligning -> Calling
     --run_step Trimming_QC  -> Trimming_QC -> Aligning -> Calling
     --run_step Aligning     -> Aligning -> Calling
     --run_step Calling      -> Calling

Required inputs:
  --input_dir   Input FASTQ directory                          [required]
  --ref         Reference FASTA                                [required]

Input example:
  PE input_dir:
    SimA_R1.fastq.gz
    SimA_R2.fastq.gz
    SimB_R1.fastq.gz
    SimB_R2.fastq.gz
    SimC_R1.fastq.gz
    SimC_R2.fastq.gz
  SE input_dir:
    SimSingleA.fastq.gz
    SimSingleB.fq.gz

Core controls:
  --read_type         SE|PE                                    [default: PE]
  --caller            bcftools|gatk                            [default: bcftools]
  --use_bqsr          true|false                               [default: false]
  --bqsr_panel        External known-sites VCF/VCF.GZ for BQSR [default: null]
  --run_step          Raw_QC|Trimming_QC|Aligning|Calling|all [default: all]
  --help              Show this help and exit                  [default: false]
  --fastp_parameters  Extra fastp options string               [default: '']
  --bwa_parameters    Extra bwa mem options string             [default: '']
  --bwamem2_parameters Extra bwa-mem2 mem options string       [default: '']
  --bcftools_mpileup_parameters Extra bcftools mpileup options [default: '']
  --bcftools_call_parameters Extra bcftools call options        [default: '']
  --bcftools_concat_parameters Extra bcftools concat options    [default: '']
  --gatk_haplotypecaller_parameters Extra GATK HaplotypeCaller options [default: '']
  --gatk_genomicsdbimport_parameters Extra GATK GenomicsDBImport options [default: '']
  --gatk_genotypegvcfs_parameters Extra GATK GenotypeGVCFs options [default: '']
  --gatk_baserecalibrator_parameters Extra GATK BaseRecalibrator options [default: '']
  --gatk_applybqsr_parameters Extra GATK ApplyBQSR options      [default: '']
  --gatk_variantfiltration_snp_parameters Extra GATK SNP VariantFiltration options [default: '']
  --gatk_variantfiltration_indel_parameters Extra GATK INDEL VariantFiltration options [default: '']
  --sif               Container image path                     [default: ./WGS_Variant_Calling.sif]

Profile switching:
  -profile local      Local executor (default profile behavior)
  -profile slurm      SLURM scheduler
  -profile awsbatch   AWS Batch scheduler

Scheduler configuration:
  SLURM and AWS Batch resource defaults are defined in nextflow.config.
  Edit nextflow.config to change global queue/cpu/memory/time or per-process withName overrides.

Workflow steps:
  Raw_QC:
    FastQC (raw FASTQ) -> MultiQC
  Trimming_QC:
    fastp -> FastQC (post-fastp) -> MultiQC
  Aligning:
    SE: BWA index -> bwa mem | samtools sort | samtools markdup -r
    PE: bwa-mem2 index -> bwa-mem2 mem | samtools sort -n | samtools fixmate | samtools sort | samtools markdup -r
  Calling (parallel by chromosome/scaffold from reference headers):
    bcftools: per-chrom call -> merge
    gatk(use_bqsr=false): per-chrom HC/GDB/Genotype -> merge raw VCF -> SNP/INDEL hard filter -> merge
    gatk(use_bqsr=true): external known-sites panel -> FULL-BAM BQSR -> per-chrom HC/GDB/Genotype -> merge raw VCF -> SNP/INDEL hard filter -> merge

Step-start behavior:
  --run_step Raw_QC       starts at Raw_QC and runs to the end
  --run_step Trimming_QC  starts at Trimming_QC and runs to the end
  --run_step Aligning     starts at Aligning and runs to the end; if Step2 outputs are incomplete, it auto-runs Trimming_QC once
  --run_step Calling      starts at Calling; if Step3 outputs are incomplete, it auto-runs Aligning once
                          (no multi-level fallback beyond one step)

Output root:
  results/

Examples:
  From test_run/run:
  nextflow run ../../main.nf -profile local --input_dir ../data --ref ../data/sim_ref_100kb.fa --run_step all
  nextflow run ../../main.nf -profile slurm --input_dir ../data --ref ../data/sim_ref_100kb.fa --run_step Calling --caller gatk
  nextflow run ../../main.nf -profile awsbatch --input_dir ../data --ref ../data/sim_ref_100kb.fa --run_step Calling --caller bcftools
  nextflow run ../../main.nf --input_dir ../data --ref ../data/sim_ref_100kb.fa --read_type PE --run_step Raw_QC
  nextflow run ../../main.nf --input_dir ../data --ref ../data/sim_ref_100kb.fa --read_type PE --run_step Trimming_QC --fastp_parameters '--qualified_quality_phred 20 --length_required 50'
  nextflow run ../../main.nf --input_dir ../data --ref ../data/sim_ref_100kb.fa --read_type PE --run_step Aligning
  nextflow run ../../main.nf --input_dir ../data --ref ../data/sim_ref_100kb.fa --read_type PE --run_step Calling --caller gatk --use_bqsr true --bqsr_panel /path/to/known_sites.vcf.gz
""".stripIndent()
}

// Validate user-facing parameters before channel creation.
def validate_params() {
    if (!(params.read_type in ['SE', 'PE'])) {
        error "Invalid --read_type: ${params.read_type}. Allowed: SE, PE"
    }
    if (!(params.caller in ['bcftools', 'gatk'])) {
        error "Invalid --caller: ${params.caller}. Allowed: bcftools, gatk"
    }
    if (params.use_bqsr && params.caller != 'gatk') {
        error "--use_bqsr is only supported when --caller gatk"
    }
    if (!(params.run_step in ['Raw_QC', 'Trimming_QC', 'Aligning', 'Calling', 'all'])) {
        error "Invalid --run_step: ${params.run_step}. Allowed: Raw_QC,Trimming_QC,Aligning,Calling,all"
    }
    if (!params.input_dir && !params.help) {
        error "Missing required parameter: --input_dir"
    }
    if (!params.ref && !params.help) {
        error "Missing required parameter: --ref"
    }
    if (!params.ref_abs && !params.help) {
        error "Failed to resolve absolute reference path from --ref: ${params.ref}"
    }
    if (!params.ref_base && !params.help) {
        error "Failed to resolve reference basename from --ref: ${params.ref}"
    }
    if (params.use_bqsr && !params.bqsr_panel) {
        error "Missing required parameter: --bqsr_panel (required when --caller gatk --use_bqsr true)"
    }
    if (params.use_bqsr && !params.bqsr_panel_abs) {
        error "Failed to resolve absolute BQSR panel path from --bqsr_panel: ${params.bqsr_panel}"
    }
    if (params.use_bqsr && !file(params.bqsr_panel_abs).exists()) {
        error "BQSR panel file not found: ${params.bqsr_panel_abs}"
    }
}

// Build raw FASTQ channel for Step1: tuple(sample, read_label, fastq).
def build_raw_fastq_channel() {
    if (params.read_type == 'SE') {
        return Channel
            .fromPath("${params.input_dir}/*.fastq")
            .mix(Channel.fromPath("${params.input_dir}/*.fastq.gz"))
            .mix(Channel.fromPath("${params.input_dir}/*.fq"))
            .mix(Channel.fromPath("${params.input_dir}/*.fq.gz"))
            .ifEmpty { error "No SE FASTQ files found under: ${params.input_dir}" }
            .map { fq ->
                def sample = fq.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
                tuple(sample, 'SE', fq)
            }
    }

    def chR1 = Channel
        .fromPath("${params.input_dir}/*_R1.fastq")
        .mix(Channel.fromPath("${params.input_dir}/*_R1.fastq.gz"))
        .mix(Channel.fromPath("${params.input_dir}/*_R1.fq"))
        .mix(Channel.fromPath("${params.input_dir}/*_R1.fq.gz"))
        .ifEmpty { error "No PE R1 FASTQ files found under: ${params.input_dir}" }

    return chR1.flatMap { r1 ->
        def matcher = (r1.name =~ /^(.+)_R1\.(fastq|fq)(\.gz)?$/)
        if (!matcher.matches()) {
            error "Unsupported PE R1 filename pattern: ${r1.name}"
        }
        def sample = matcher[0][1]
        def ext = "${matcher[0][2]}${matcher[0][3] ?: ''}"
        def r2 = file("${r1.parent}/${sample}_R2.${ext}")
        if (!r2.exists()) {
            error "Missing R2 for sample ${sample}: expected ${r2}"
        }
        [tuple(sample, 'R1', r1), tuple(sample, 'R2', r2)]
    }
}

// Build fastp input channel for Step2/3/4 minimal dependency runs.
def build_fastp_input_channel() {
    if (params.read_type == 'SE') {
        return Channel
            .fromPath("${params.input_dir}/*.fastq")
            .mix(Channel.fromPath("${params.input_dir}/*.fastq.gz"))
            .mix(Channel.fromPath("${params.input_dir}/*.fq"))
            .mix(Channel.fromPath("${params.input_dir}/*.fq.gz"))
            .ifEmpty { error "No SE FASTQ files found under: ${params.input_dir}" }
            .map { fq ->
                def sample = fq.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
                tuple(sample, fq)
            }
    }

    def chR1 = Channel
        .fromPath("${params.input_dir}/*_R1.fastq")
        .mix(Channel.fromPath("${params.input_dir}/*_R1.fastq.gz"))
        .mix(Channel.fromPath("${params.input_dir}/*_R1.fq"))
        .mix(Channel.fromPath("${params.input_dir}/*_R1.fq.gz"))
        .ifEmpty { error "No PE R1 FASTQ files found under: ${params.input_dir}" }

    return chR1.map { r1 ->
        def matcher = (r1.name =~ /^(.+)_R1\.(fastq|fq)(\.gz)?$/)
        if (!matcher.matches()) {
            error "Unsupported PE R1 filename pattern: ${r1.name}"
        }
        def sample = matcher[0][1]
        def ext = "${matcher[0][2]}${matcher[0][3] ?: ''}"
        def r2 = file("${r1.parent}/${sample}_R2.${ext}")
        if (!r2.exists()) {
            error "Missing R2 for sample ${sample}: expected ${r2}"
        }
        tuple(sample, r1, r2)
    }
}

// Collect sample IDs from input FASTQ files using the same naming rules.
def collect_sample_ids_from_input() {
    def samples = [] as Set
    if (params.read_type == 'SE') {
        file(params.input_dir).listFiles()?.each { f ->
            if (f.name ==~ /.+\.(fastq|fq)(\.gz)?$/) {
                samples << f.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
            }
        }
        return samples as List
    }

    file(params.input_dir).listFiles()?.each { f ->
        if (f.name ==~ /.+_R1\.(fastq|fq)(\.gz)?$/) {
            samples << f.name.replaceFirst(/_R1\.(fastq|fq)(\.gz)?$/, '')
        }
    }
    return samples as List
}

// Check whether Step2 trimmed FASTQ outputs exist for all input samples.
def has_complete_step2_outputs() {
    def samples = collect_sample_ids_from_input()
    if (!samples) return false
    def step2Dir = file('results/02_fastp_qc/fastp')
    if (!step2Dir.exists()) return false

    return samples.every { s ->
        if (params.read_type == 'SE') {
            file("${step2Dir}/${s}.clean.fastq.gz").exists()
        } else {
            file("${step2Dir}/${s}.clean.R1.fastq.gz").exists() &&
            file("${step2Dir}/${s}.clean.R2.fastq.gz").exists()
        }
    }
}

// Build cleaned reads channel directly from Step2 published outputs.
def build_clean_reads_from_step2_results() {
    def samples = collect_sample_ids_from_input()
    if (!samples) {
        error "No samples detected from input_dir for Step2 reuse: ${params.input_dir}"
    }
    def step2Dir = file('results/02_fastp_qc/fastp')

    def tuples = samples.collect { s ->
        if (params.read_type == 'SE') {
            def fq = file("${step2Dir}/${s}.clean.fastq.gz")
            if (!fq.exists()) error "Missing Step2 output for sample ${s}: ${fq}"
            tuple(s, fq)
        } else {
            def r1 = file("${step2Dir}/${s}.clean.R1.fastq.gz")
            def r2 = file("${step2Dir}/${s}.clean.R2.fastq.gz")
            if (!r1.exists() || !r2.exists()) {
                error "Missing Step2 output for sample ${s}: expected ${r1} and ${r2}"
            }
            tuple(s, r1, r2)
        }
    }
    return Channel.from(tuples)
}

// Check whether Step3 markdup BAM outputs exist for all input samples.
def has_complete_step3_outputs() {
    def samples = collect_sample_ids_from_input()
    if (!samples) return false
    def step3Dir = file('results/03_markdup')
    if (!step3Dir.exists()) return false

    return samples.every { s ->
        file("${step3Dir}/${s}.markdup.bam").exists() &&
        file("${step3Dir}/${s}.markdup.bam.bai").exists()
    }
}

// Build markdup BAM channel directly from Step3 published outputs.
def build_markdup_from_step3_results() {
    def samples = collect_sample_ids_from_input()
    if (!samples) {
        error "No samples detected from input_dir for Step3 reuse: ${params.input_dir}"
    }
    def step3Dir = file('results/03_markdup')

    def tuples = samples.collect { s ->
        def bam = file("${step3Dir}/${s}.markdup.bam")
        def bai = file("${step3Dir}/${s}.markdup.bam.bai")
        if (!bam.exists() || !bai.exists()) {
            error "Missing Step3 output for sample ${s}: expected ${bam} and ${bai}"
        }
        tuple(s, bam, bai)
    }
    return Channel.from(tuples)
}

// Build interval channel (idx, chrom) from reference FASTA headers for per-chrom parallel calling.
def build_interval_channel() {
    def intervals = []
    int idx = 1
    file(params.ref_abs).eachLine { line ->
        if (line.startsWith('>')) {
            def chrom = line.substring(1).tokenize(' \t')[0]
            def interval_idx = String.format('%04d', idx)
            intervals << tuple(interval_idx, chrom)
            idx++
        }
    }

    if (!intervals) {
        error "No chromosome/scaffold headers found in reference: ${params.ref_abs}"
    }

    return Channel.from(intervals)
}

// Step3 dispatcher: SE and PE branches converge on markdup BAM.
def run_step3(ch_clean_reads) {
    def ch_ref_fa = Channel.value(file(params.ref_abs))

    if (params.read_type == 'SE') {
        ref_index = BWA_INDEX(ch_ref_fa)
        return ALIGN_MARKDUP_SE(ch_clean_reads, ref_index.ref_bundle)
    }

    ref_index = BWAMEM2_INDEX(ch_ref_fa)
    return ALIGN_MARKDUP_PE(ch_clean_reads, ref_index.ref_bundle)
}

// Step4 dispatcher: all calling branches run per chromosome/scaffold and then merge.
def run_step4(ch_markdup) {
    def ch_ref_fa = Channel.value(file(params.ref_abs))
    def ch_intervals = build_interval_channel()

    if (params.caller == 'bcftools') {
        return BCFTOOLS_CALL(ch_markdup, ch_ref_fa, ch_intervals)
    }

    def hc_input = ch_markdup

    if (params.use_bqsr) {
        def ch_bqsr_panel = Channel.value(file(params.bqsr_panel_abs))
        full_bam_bqsr = GATK_BQSR(ch_markdup, ch_bqsr_panel, ch_ref_fa)
        hc_input = full_bam_bqsr.recal_bam
    }

    hc = GATK_HAPLOTYPECALLER(hc_input, ch_ref_fa, ch_intervals)
    gdb = GATK_GENOMICSDBIMPORT(hc.gvcf, ch_ref_fa)
    gt = GATK_GENOTYPEGVCFS(gdb.gendb, ch_ref_fa)
    merged_raw = GATK_MERGE_RAW(gt.raw_vcf)
    return GATK_FILTER(merged_raw.raw_merged, ch_ref_fa)
}

def run_trimming_stage(step1_done_signal = null) {
    ch_fastp_input = build_fastp_input_channel()
    if (step1_done_signal != null) {
        ch_fastp_input = ch_fastp_input.combine(step1_done_signal).map { row ->
            def data = row.take(row.size() - 1)
            tuple(*data)
        }
    }

    trim = FASTP(ch_fastp_input)
    post_qc = FASTQC_POST(trim.cleaned_for_qc)
    ch_post_fastqc_zip = post_qc.fastqc_reports.map { sample, read_label, zip_file, html_file -> zip_file }
    MULTIQC_POST(ch_post_fastqc_zip)
    return trim
}

workflow {
    if (params.help) {
        print_help()
        exit 0
    }

    validate_params()

    def step1_done_signal = null
    def trim = null
    def markdup = null

    if (params.run_step in ['Raw_QC', 'all']) {
        ch_raw_fastq = build_raw_fastq_channel()
        raw_qc = FASTQC_RAW(ch_raw_fastq)
        ch_fastqc_zip = raw_qc.fastqc_reports.map { sample, read_label, zip_file, html_file -> zip_file }
        raw_mqc = MULTIQC_RAW(ch_fastqc_zip)
        step1_done_signal = raw_mqc.multiqc_report
    }

    if (params.run_step == 'all') {
        trim = run_trimming_stage(step1_done_signal)
        markdup = run_step3(trim.cleaned_reads)
        run_step4(markdup.markdup_bam)
    }

    if (params.run_step == 'Raw_QC') {
        trim = run_trimming_stage(step1_done_signal)
        markdup = run_step3(trim.cleaned_reads)
        run_step4(markdup.markdup_bam)
    }

    if (params.run_step == 'Trimming_QC') {
        trim = run_trimming_stage()
        markdup = run_step3(trim.cleaned_reads)
        run_step4(markdup.markdup_bam)
    }

    if (params.run_step == 'Aligning') {
        def clean_reads_for_align
        if (has_complete_step2_outputs()) {
            log.info "run_step=Aligning: reusing Step2 outputs from results/02_fastp_qc/fastp"
            clean_reads_for_align = build_clean_reads_from_step2_results()
        } else {
            log.info "run_step=Aligning: Step2 outputs incomplete, auto-running Trimming_QC once"
            trim = run_trimming_stage()
            clean_reads_for_align = trim.cleaned_reads
        }

        markdup = run_step3(clean_reads_for_align)
        run_step4(markdup.markdup_bam)
    }

    if (params.run_step == 'Calling') {
        def markdup_for_calling
        if (has_complete_step3_outputs()) {
            log.info "run_step=Calling: reusing Step3 outputs from results/03_markdup"
            markdup_for_calling = build_markdup_from_step3_results()
        } else {
            log.info "run_step=Calling: Step3 outputs incomplete, auto-running Aligning once"
            if (!has_complete_step2_outputs()) {
                error "run_step=Calling can auto-fallback only one step to Aligning, but Step2 outputs are incomplete under results/02_fastp_qc/fastp"
            }
            markdup = run_step3(build_clean_reads_from_step2_results())
            markdup_for_calling = markdup.markdup_bam
        }
        run_step4(markdup_for_calling)
    }
}
