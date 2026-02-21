nextflow.enable.dsl = 2

params.help             = params.help ?: false
params.input_dir        = params.input_dir ?: null
params.ref              = params.ref ?: null
params.read_type        = params.read_type ?: 'PE'
params.caller           = params.caller ?: 'bcftools'
params.use_bqsr         = params.use_bqsr ?: false
params.run_step         = params.run_step ?: 'all'
params.sif              = params.sif ?: 'WGS_Variant_Calling.sif'
params.fastp_parameters = params.fastp_parameters ?: ''
params.ref_abs          = params.ref ? file(params.ref).toAbsolutePath().toString() : null
params.ref_base         = params.ref_abs ? file(params.ref_abs).getName() : null

params.threads        = (params.threads ?: 4) as Integer
params.fastqc_cpus    = (params.fastqc_cpus ?: params.threads) as Integer
params.fastp_cpus     = (params.fastp_cpus ?: params.threads) as Integer
params.bwa_cpus       = (params.bwa_cpus ?: params.threads) as Integer
params.bwamem2_cpus   = (params.bwamem2_cpus ?: params.threads) as Integer
params.sambamba_cpus  = (params.sambamba_cpus ?: params.threads) as Integer
params.bcftools_cpus  = (params.bcftools_cpus ?: params.threads) as Integer
params.gatk_cpus      = (params.gatk_cpus ?: params.threads) as Integer

include { FASTQC_RAW }  from './modules/step1_raw_qc/fastqc_raw'
include { MULTIQC_RAW } from './modules/step1_raw_qc/multiqc_raw'

include { FASTP }        from './modules/step2_fastp_qc/fastp'
include { FASTQC_POST }  from './modules/step2_fastp_qc/fastqc_post'
include { MULTIQC_POST } from './modules/step2_fastp_qc/multiqc_post'

include { BWA_SE }            from './modules/step3_alignment/bwa_se'
include { BWAMEM2_PE }        from './modules/step3_alignment/bwa_mem2_pe'
include { BWA_INDEX }         from './modules/step3_alignment/bwa_index'
include { BWAMEM2_INDEX }     from './modules/step3_alignment/bwa_mem2_index'
include { FIXMATE }           from './modules/step3_alignment/fixmate'
include { SORT_BAM }          from './modules/step3_alignment/sort'
include { SAMBAMBA_MARKDUP }  from './modules/step3_alignment/sambamba_markdup'

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
  1) Full pipeline (default): Raw_QC -> Trimming_QC -> Aligning -> Calling
  2) Single-step: --run_step Raw_QC|Trimming_QC|Aligning|Calling

Required inputs:
  --input_dir   Input FASTQ directory
  --ref         Reference FASTA

Core controls:
  --read_type         SE|PE                                    [default: PE]
  --caller            bcftools|gatk                            [default: bcftools]
  --use_bqsr          true|false                               [default: false]
  --run_step          Raw_QC|Trimming_QC|Aligning|Calling|all [default: all]
  --help              Show this help and exit                  [default: false]
  --fastp_parameters  Extra fastp options string               [default: '']
  --sif               Container image path                     [default: WGS_Variant_Calling.sif]

Thread controls:
  --threads        Global fallback threads             [default: 4]
  --fastqc_cpus    FastQC threads (overrides --threads)
  --fastp_cpus     fastp threads (overrides --threads)
  --bwa_cpus       BWA threads (overrides --threads)
  --bwamem2_cpus   bwa-mem2 threads (overrides --threads)
  --sambamba_cpus  sambamba threads (overrides --threads)
  --bcftools_cpus  bcftools threads (overrides --threads)
  --gatk_cpus      GATK threads (overrides --threads)

Workflow steps:
  Raw_QC:
    FastQC (raw FASTQ) -> MultiQC
  Trimming_QC:
    fastp -> FastQC (post-fastp) -> MultiQC
  Aligning:
    SE: BWA index -> BWA -> sort -> sambamba markdup
    PE: bwa-mem2 index -> bwa-mem2 -> fixmate -> sort -> sambamba markdup
  Calling (parallel by chromosome/scaffold from reference headers):
    bcftools: per-chrom call -> merge
    gatk(use_bqsr=false): per-chrom HC/GDB/Genotype -> merge raw VCF -> SNP/INDEL hard filter -> merge
    gatk(use_bqsr=true): per-chrom bcftools seed->merge -> high-conf SNP -> FULL-BAM BQSR -> per-chrom HC/GDB/Genotype -> merge raw VCF -> SNP/INDEL hard filter -> merge

Output root:
  results/

Examples:
  nextflow run main.nf --input_dir data --ref ref.fa --read_type PE --run_step Raw_QC
  nextflow run main.nf --input_dir data --ref ref.fa --read_type PE --run_step Trimming_QC --fastp_parameters '--qualified_quality_phred 20 --length_required 50'
  nextflow run main.nf --input_dir data --ref ref.fa --read_type PE --run_step Aligning
  nextflow run main.nf --input_dir data --ref ref.fa --read_type PE --run_step Calling --caller gatk --use_bqsr true
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
}

// Build raw FASTQ channel for Step1: tuple(sample, read_label, fastq).
def build_raw_fastq_channel() {
    if (params.read_type == 'SE') {
        return Channel
            .fromPath("${params.input_dir}/*.fastq.gz")
            .mix(Channel.fromPath("${params.input_dir}/*.fq.gz"))
            .ifEmpty { error "No SE FASTQ files found under: ${params.input_dir}" }
            .map { fq ->
                def sample = fq.baseName.replaceFirst(/\.(fastq|fq)$/, '')
                tuple(sample, 'SE', fq)
            }
    }

    def chR1 = Channel
        .fromPath("${params.input_dir}/*_R1.fastq.gz")
        .mix(Channel.fromPath("${params.input_dir}/*_R1.fq.gz"))
        .ifEmpty { error "No PE R1 FASTQ files found under: ${params.input_dir}" }

    return chR1.flatMap { r1 ->
        def sample = r1.name.replaceFirst(/_R1\.(fastq|fq)\.gz$/, '')
        def ext = r1.name.endsWith('.fastq.gz') ? 'fastq.gz' : 'fq.gz'
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
            .fromPath("${params.input_dir}/*.fastq.gz")
            .mix(Channel.fromPath("${params.input_dir}/*.fq.gz"))
            .ifEmpty { error "No SE FASTQ files found under: ${params.input_dir}" }
            .map { fq ->
                def sample = fq.baseName.replaceFirst(/\.(fastq|fq)$/, '')
                tuple(sample, fq)
            }
    }

    def chR1 = Channel
        .fromPath("${params.input_dir}/*_R1.fastq.gz")
        .mix(Channel.fromPath("${params.input_dir}/*_R1.fq.gz"))
        .ifEmpty { error "No PE R1 FASTQ files found under: ${params.input_dir}" }

    return chR1.map { r1 ->
        def sample = r1.name.replaceFirst(/_R1\.(fastq|fq)\.gz$/, '')
        def ext = r1.name.endsWith('.fastq.gz') ? 'fastq.gz' : 'fq.gz'
        def r2 = file("${r1.parent}/${sample}_R2.${ext}")
        if (!r2.exists()) {
            error "Missing R2 for sample ${sample}: expected ${r2}"
        }
        tuple(sample, r1, r2)
    }
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
        se_aln = BWA_SE(ch_clean_reads, ref_index.ref_bundle)
        sorted = SORT_BAM(se_aln.bam_for_sort)
        return SAMBAMBA_MARKDUP(sorted.sorted_bam)
    }

    ref_index = BWAMEM2_INDEX(ch_ref_fa)
    pe_aln = BWAMEM2_PE(ch_clean_reads, ref_index.ref_bundle)
    fixed = FIXMATE(pe_aln.bam_for_fixmate)
    sorted = SORT_BAM(fixed.bam_for_sort)
    return SAMBAMBA_MARKDUP(sorted.sorted_bam)
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
        // 1) Seed calling is parallelized by chromosome/scaffold then merged into one genome-wide seed VCF.
        seed_call = BCFTOOLS_CALL(ch_markdup, ch_ref_fa, ch_intervals)
        // 2) High-confidence SNP extraction + BQSR are applied once on the full sample BAM.
        full_bam_bqsr = GATK_BQSR(ch_markdup, seed_call.bcf_vcf, ch_ref_fa)
        hc_input = full_bam_bqsr.recal_bam
    }

    hc = GATK_HAPLOTYPECALLER(hc_input, ch_ref_fa, ch_intervals)
    gdb = GATK_GENOMICSDBIMPORT(hc.gvcf, ch_ref_fa)
    gt = GATK_GENOTYPEGVCFS(gdb.gendb, ch_ref_fa)
    merged_raw = GATK_MERGE_RAW(gt.raw_vcf)
    return GATK_FILTER(merged_raw.raw_merged, ch_ref_fa)
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

    if (params.run_step in ['Trimming_QC', 'all']) {
        ch_fastp_input = build_fastp_input_channel()
        if (params.run_step == 'all') {
            ch_fastp_input = ch_fastp_input.combine(step1_done_signal).map { reads, done -> reads }
        }

        trim = FASTP(ch_fastp_input)
        post_qc = FASTQC_POST(trim.cleaned_for_qc)
        ch_post_fastqc_zip = post_qc.fastqc_reports.map { sample, read_label, zip_file, html_file -> zip_file }
        MULTIQC_POST(ch_post_fastqc_zip)
    }

    if (params.run_step == 'Aligning') {
        trim = FASTP(build_fastp_input_channel())
        markdup = run_step3(trim.cleaned_reads)
    }

    if (params.run_step == 'Calling') {
        trim = FASTP(build_fastp_input_channel())
        markdup = run_step3(trim.cleaned_reads)
        run_step4(markdup.markdup_bam)
    }

    if (params.run_step == 'all') {
        markdup = run_step3(trim.cleaned_reads)
        run_step4(markdup.markdup_bam)
    }
}
