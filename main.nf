nextflow.enable.dsl = 2

params.help             = params.containsKey('help') ? params.get('help') : false
params.help_full        = params.containsKey('help_full') ? params.get('help_full') : false
params.input_dir        = params.input_dir ?: null
params.input_vcf        = params.input_vcf ?: null
params.SRA_list         = params.containsKey('SRA_list') ? params.get('SRA_list') : null
params.trimmed_input_dir = params.trimmed_input_dir ?: null
params.ref              = params.ref ?: null
params.read_type        = params.read_type ?: 'PE'
params.caller           = params.caller ?: 'bcftools'
params.use_bqsr         = params.use_bqsr ?: false
params.cleanup_GVCFS    = params.cleanup_GVCFS ?: false
params.enable_vcffilter = params.enable_vcffilter ?: false
params.vcffilter_stop_after = params.vcffilter_stop_after ?: 'final'
params.interval_size    = params.interval_size ?: null
params.bqsr_panel       = params.containsKey('bqsr_panel') ? params.get('bqsr_panel') : null
params.run_step         = params.run_step ?: 'all'
params.image            = params.image ?: 'leonardliu0910/wgs_variant_calling:latest'
params.sif              = params.sif ?: null
params.container_image  = params.sif ?: params.image
params.fastp_parameters = params.fastp_parameters ?: ''
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


include { DOWNLOAD_SRA; RAW_QC }  from './subworkflows/raw_qc'
include { TRIMMING_QC }           from './subworkflows/trimming_qc'
include { ALIGNMENT }             from './subworkflows/alignment'
include { VARIANT_CALLING }       from './subworkflows/variant_calling'
include { VCF_FILTERING }         from './subworkflows/vcf_filtering'

def print_help(boolean full = false) {
    def trimmingAdvanced = full ? """
  Advanced:
    --fastp_parameters  Extra fastp options string             [default: '']
""" : ""

    def aligningAdvanced = full ? """
  Advanced:
    --bwamem2_parameters Extra bwa-mem2 mem options string     [default: '']
""" : ""

    def callingAdvanced = full ? """
  Advanced:
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
""" : ""

    def vcfFilteringAdvanced = full ? """
  Config-driven parameters:
    VCF filtering phase parameters are read from vcffilter.config
    Current scope: biallelic SNP only
    Current implementation status: phase1 and phase2 available
""" : ""

    log.info """
NGS Variant Calling Nextflow Pipeline (DSL2, Nextflow 25.10.4)

Execution modes:
  1) Full pipeline (default): --run_step all -> Raw_QC -> Trimming_QC -> Aligning -> Calling
  2) Start from a specific step and continue to the end:
     --run_step Raw_QC       -> Raw_QC only
     --run_step Trimming_QC  -> Trimming_QC -> Aligning -> Calling [-> VCF_Filtering if enabled]
     --run_step Aligning     -> Aligning -> Calling [-> VCF_Filtering if enabled]
     --run_step Calling      -> Calling [-> VCF_Filtering if enabled]
     --run_step VCF_Filtering -> VCF_Filtering only

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

Parameters by scope:

Global:
  --run_step          Raw_QC|Trimming_QC|Aligning|Calling|VCF_Filtering|all [default: all]
  --read_type         SE|PE                                    [default: PE]
  --ref               Reference FASTA                          [required for Raw_QC/Trimming_QC/Aligning/Calling/all; not required for standalone VCF_Filtering]
  --help              Show this help and exit                  [default: false]
  --help_full         Show full help with advanced options     [default: false]
  --image             Docker/registry image                    [default: leonardliu0910/wgs_variant_calling:latest]
  --sif               Local Singularity/Apptainer image path   [default: none; overrides --image for those engines]


Step1 Raw_QC:
  Inputs:
    --input_dir       Raw FASTQ directory                      [raw input mode 1]
    --SRA_list        Two-column SRA list file                 [raw input mode 2]
  Requirement:
    For --run_step Raw_QC or --run_step all, provide exactly one of --input_dir or --SRA_list


Step2 Trimming_QC:
  Inputs:
    --input_dir       Raw FASTQ directory                      [required when starting at Trimming_QC]
  Requirement:
    For --run_step Trimming_QC, --input_dir is required
${trimmingAdvanced}\


Step3 Aligning:
  Inputs:
    --trimmed_input_dir Trimmed FASTQ directory                [optional direct input]
  Behavior:
    If --trimmed_input_dir is not provided, the workflow reuses Step2 outputs from results/02_fastp_qc/fastp
    If Step2 outputs are missing and --run_step Aligning is used, the workflow auto-runs Trimming_QC once
${aligningAdvanced}\


Step4 Calling:
  Inputs:
    Step3 CRAM outputs under results/03_markdup               [preferred when available]
    --trimmed_input_dir Trimmed FASTQ directory                [optional one-step fallback to Aligning]
    --ref              Reference FASTA                         [required]
    --caller            bcftools|gatk                            [default: bcftools]
    --interval_size     Interval block size for calling          [required e.g：5K or 5M or 5G or 5000000]
    --use_bqsr          true|false                               [default: false]
    --bqsr_panel        External known-sites VCF/VCF.GZ for BQSR [required when BQSR is enabled]
    --cleanup_GVCFS     Delete interval gVCFs after GenomicsDBImport succeeds [default: false]
  Requirement:
    For any path that executes Calling, --interval_size is required
${callingAdvanced}\


Step5 VCF_Filtering:
  Inputs:
    Step4 cohort VCF                                          [used when --enable_vcffilter true]
    --input_vcf       Input VCF.GZ                             [required for standalone VCF_Filtering]
    --enable_vcffilter  true|false                               [default: false]
    --vcffilter_stop_after phase1|final                          [default: final]
  Notes:
    supports biallelic SNP only
    does not support indel or mixed SNP+indel VCF filtering
    standalone mode expects a tabix index at INPUT_VCF.tbi
  Requirement:
    For --run_step VCF_Filtering, --input_vcf is required
${vcfFilteringAdvanced}\


Profile switching:
  -profile local      Local executor (default profile behavior; Apptainer enabled unless overridden)
  -profile slurm      SLURM scheduler
  -profile awsbatch   AWS Batch scheduler


Workflow steps:
  Raw_QC:
    FastQC (raw FASTQ) -> MultiQC
  Trimming_QC:
    fastp -> FastQC (post-fastp) -> MultiQC
  Aligning:
    SE: bwa-mem2 index -> bwa-mem2 mem | samtools sort -> sambamba markdup -> reference-embedded CRAM
    PE: bwa-mem2 index -> bwa-mem2 mem | samtools collate | samtools fixmate | samtools sort -> sambamba markdup -> reference-embedded CRAM
  Calling (parallel by interval_size blocks derived from reference sequence lengths):
    bcftools: per-interval call -> merge
    gatk(use_bqsr=false): per-interval HC/GDB/Genotype -> merge raw VCF -> SNP/INDEL hard filter -> merge
    gatk(use_bqsr=true): external known-sites panel -> full-alignment BQSR -> per-interval HC/GDB/Genotype -> merge raw VCF -> SNP/INDEL hard filter -> merge
  VCF_Filtering:
    optional post-calling workflow
    phase1: biallelic SNP filter -> optional mislabeled-sample removal -> GQ/DP masking -> optional AB masking -> sample-specific DP masking -> review reports
    phase2: remove high-missing samples -> final site filter -> final reports

Step-start behavior:
  --run_step Raw_QC       runs Raw_QC only and exits
  --run_step Trimming_QC  starts at Trimming_QC and runs to the end; Step5 runs only if --enable_vcffilter true
  --run_step Aligning     starts at Aligning and runs to the end; if Step2 outputs are incomplete, it auto-runs Trimming_QC once
                          Step5 runs only if --enable_vcffilter true
  --run_step Calling      starts at Calling; if Step3 outputs are incomplete, it auto-runs Aligning once
                          Step5 runs only if --enable_vcffilter true
                          (no multi-level fallback beyond one step)
  --run_step VCF_Filtering runs Step5 only from --input_vcf

Output root:
  results/
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
    if (!(params.run_step in ['Raw_QC', 'Trimming_QC', 'Aligning', 'Calling', 'VCF_Filtering', 'all'])) {
        error "Invalid --run_step: ${params.run_step}. Allowed: Raw_QC,Trimming_QC,Aligning,Calling,VCF_Filtering,all"
    }
    if (!(params.vcffilter_stop_after in ['phase1', 'final'])) {
        error "Invalid --vcffilter_stop_after: ${params.vcffilter_stop_after}. Allowed: phase1, final"
    }
    if (params.input_vcf && !file(params.input_vcf).exists()) {
        error "input_vcf not found: ${params.input_vcf}"
    }
    if (params.run_step == 'VCF_Filtering' && params.input_vcf && !file("${params.input_vcf}.tbi").exists()) {
        error "Standalone VCF filtering requires a tabix index next to the VCF: ${params.input_vcf}.tbi"
    }
    if (!params.ref && params.run_step != 'VCF_Filtering' && !params.help && !params.help_full) {
        error "Missing required parameter: --ref"
    }
    if (!params.ref_abs && params.run_step != 'VCF_Filtering' && !params.help && !params.help_full) {
        error "Failed to resolve absolute reference path from --ref: ${params.ref}"
    }
    if (!params.ref_base && params.run_step != 'VCF_Filtering' && !params.help && !params.help_full) {
        error "Failed to resolve reference basename from --ref: ${params.ref}"
    }
    if (params.input_dir && !file(params.input_dir).exists()) {
        error "input_dir not found: ${params.input_dir}"
    }
    if (params.SRA_list && !file(params.SRA_list).exists()) {
        error "SRA_list not found: ${params.SRA_list}"
    }
    if (params.trimmed_input_dir && !file(params.trimmed_input_dir).exists()) {
        error "trimmed_input_dir not found: ${params.trimmed_input_dir}"
    }
    def needsRawInputDir = params.run_step in ['Raw_QC', 'Trimming_QC', 'all']
    def needsAnyFastqInput = params.run_step in ['Aligning', 'Calling']
    def needsCallingIntervals = params.run_step in ['Trimming_QC', 'Aligning', 'Calling', 'all']
    def usesRawEntry = params.run_step in ['Raw_QC', 'all']
    def usesStandaloneVcfFiltering = params.run_step == 'VCF_Filtering'
    if (usesRawEntry && params.input_dir && params.SRA_list) {
        error "Provide only one raw input mode: --input_dir or --SRA_list"
    }
    if (usesRawEntry && !params.input_dir && !params.SRA_list && !params.help && !params.help_full) {
        error "Missing raw input source: provide --input_dir or --SRA_list"
    }
    if (params.run_step == 'Trimming_QC' && !params.input_dir && !params.help && !params.help_full) {
        error "Missing required parameter: --input_dir"
    }
    if (needsAnyFastqInput && !params.input_dir && !params.trimmed_input_dir && !params.help && !params.help_full) {
        error "Missing input FASTQ source: provide --input_dir or --trimmed_input_dir when --run_step ${params.run_step}"
    }
    if (usesStandaloneVcfFiltering && !params.input_vcf && !params.help && !params.help_full) {
        error "Missing required parameter: --input_vcf"
    }
    if (needsCallingIntervals && !params.interval_size && !params.help && !params.help_full) {
        error "Missing required parameter: --interval_size"
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

def build_fastq_channel_from_dir(baseDir, stageLabel = 'input_dir') {
    if (!baseDir) {
        error "No directory provided for ${stageLabel}"
    }

    if (params.read_type == 'SE') {
        return Channel
            .fromPath("${baseDir}/*.fastq")
            .mix(Channel.fromPath("${baseDir}/*.fastq.gz"))
            .mix(Channel.fromPath("${baseDir}/*.fq"))
            .mix(Channel.fromPath("${baseDir}/*.fq.gz"))
            .ifEmpty { error "No SE FASTQ files found under ${stageLabel}: ${baseDir}" }
            .map { fq ->
                def sample = fq.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
                tuple(sample, fq)
            }
    }

    def chR1 = Channel
        .fromPath("${baseDir}/*_R1.fastq")
        .mix(Channel.fromPath("${baseDir}/*_R1.fastq.gz"))
        .mix(Channel.fromPath("${baseDir}/*_R1.fq"))
        .mix(Channel.fromPath("${baseDir}/*_R1.fq.gz"))
        .ifEmpty { error "No PE R1 FASTQ files found under ${stageLabel}: ${baseDir}" }

    return chR1.map { r1 ->
        def matcher = (r1.name =~ /^(.+)_R1\.(fastq|fq)(\.gz)?$/)
        if (!matcher.matches()) {
            error "Unsupported PE R1 filename pattern in ${stageLabel}: ${r1.name}"
        }
        def sample = matcher[0][1]
        def ext = "${matcher[0][2]}${matcher[0][3] ?: ''}"
        def r2 = file("${r1.parent}/${sample}_R2.${ext}")
        if (!r2.exists()) {
            error "Missing R2 for sample ${sample} in ${stageLabel}: expected ${r2}"
        }
        tuple(sample, r1, r2)
    }
}

def collect_sample_ids_from_dir(baseDir, stageLabel = 'input_dir') {
    if (!baseDir) return []

    def samples = [] as Set
    if (params.read_type == 'SE') {
        file(baseDir).listFiles()?.each { f ->
            if (f.name ==~ /.+\.(fastq|fq)(\.gz)?$/) {
                samples << f.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
            }
        }
        return samples as List
    }

    file(baseDir).listFiles()?.each { f ->
        if (f.name ==~ /.+_R1\.(fastq|fq)(\.gz)?$/) {
            samples << f.name.replaceFirst(/_R1\.(fastq|fq)(\.gz)?$/, '')
        }
    }
    return samples as List
}

def resolve_alignment_input_dir() {
    return params.trimmed_input_dir ?: params.input_dir
}

def build_sra_sample_channel() {
    if (!params.SRA_list) {
        error "SRA_list is required to build SRA sample channel"
    }

    def rows = []
    file(params.SRA_list).eachLine { rawLine ->
        def line = rawLine.trim()
        if (!line || line.startsWith('#')) return
        def cols = line.split(/\t/)
        if (cols.size() < 2) {
            error "Invalid SRA_list row: expected at least 2 tab-separated columns, got: ${rawLine}"
        }
        def accession = cols[0].trim()
        def sampleId = cols[1].trim()
        if (!accession || !sampleId) {
            error "Invalid SRA_list row with empty accession or sample name: ${rawLine}"
        }
        rows << tuple(sampleId, accession)
    }

    if (!rows) {
        error "No usable rows found in SRA_list: ${params.SRA_list}"
    }

    return Channel.from(rows)
}

def parse_interval_size_to_bases(intervalSizeRaw) {
    if (!intervalSizeRaw) {
        error "interval_size is required for calling"
    }

    def normalized = intervalSizeRaw.toString().trim().toUpperCase()
    def matcher = (normalized =~ /^(\d+)([KMG])?$/)
    if (!matcher.matches()) {
        error "Invalid --interval_size: ${intervalSizeRaw}. Use forms like 500K, 5M, or 1000000"
    }

    long value = matcher[0][1] as long
    def suffix = matcher[0][2]
    long multiplier = 1L
    if (suffix == 'K') multiplier = 1_000L
    if (suffix == 'M') multiplier = 1_000_000L
    if (suffix == 'G') multiplier = 1_000_000_000L

    long bases = value * multiplier
    if (bases <= 0L) {
        error "interval_size must be > 0: ${intervalSizeRaw}"
    }
    return bases
}

def load_reference_lengths() {
    def refFai = file("${params.ref_abs}.fai")
    if (refFai.exists()) {
        def lengths = []
        refFai.eachLine { line ->
            def cols = line.split('\t')
            if (cols.size() >= 2) {
                lengths << [name: cols[0], length: cols[1] as long]
            }
        }
        if (lengths) {
            return lengths
        }
    }

    def lengths = []
    String currentName = null
    long currentLength = 0L
    file(params.ref_abs).eachLine { line ->
        if (line.startsWith('>')) {
            if (currentName != null) {
                lengths << [name: currentName, length: currentLength]
            }
            currentName = line.substring(1).tokenize(' \t')[0]
            currentLength = 0L
        } else {
            currentLength += line.trim().size()
        }
    }
    if (currentName != null) {
        lengths << [name: currentName, length: currentLength]
    }

    if (!lengths) {
        error "No reference sequence entries found in reference: ${params.ref_abs}"
    }
    return lengths
}

// Build fastp input channel for Step2/3/4 minimal dependency runs.
def build_fastp_input_channel() {
    return build_fastq_channel_from_dir(params.input_dir, 'input_dir')
}

// Collect sample IDs from input FASTQ files using the same naming rules.
def collect_sample_ids_from_input() {
    return collect_sample_ids_from_dir(params.input_dir, 'input_dir')
}

def collect_sample_ids_from_alignment_source() {
    def baseDir = resolve_alignment_input_dir()
    def label = params.trimmed_input_dir ? 'trimmed_input_dir' : 'input_dir'
    return collect_sample_ids_from_dir(baseDir, label)
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

def build_clean_reads_from_trimmed_input_dir() {
    return build_fastq_channel_from_dir(params.trimmed_input_dir, 'trimmed_input_dir')
}

// Check whether Step3 markdup CRAM outputs exist for all input samples.
def has_complete_step3_outputs() {
    def samples = collect_sample_ids_from_alignment_source()
    if (!samples) return false
    def step3Dir = file('results/03_markdup')
    if (!step3Dir.exists()) return false

    return samples.every { s ->
        file("${step3Dir}/${s}.markdup.cram").exists() &&
        file("${step3Dir}/${s}.markdup.cram.crai").exists()
    }
}

// Build markdup CRAM channel directly from Step3 published outputs.
def build_markdup_from_step3_results() {
    def samples = collect_sample_ids_from_alignment_source()
    if (!samples) {
        error "No samples detected from alignment input source for Step3 reuse"
    }
    def step3Dir = file('results/03_markdup')

    def tuples = samples.collect { s ->
        def cram = file("${step3Dir}/${s}.markdup.cram")
        def crai = file("${step3Dir}/${s}.markdup.cram.crai")
        if (!cram.exists() || !crai.exists()) {
            error "Missing Step3 output for sample ${s}: expected ${cram} and ${crai}"
        }
        tuple(s, cram, crai)
    }
    return Channel.from(tuples)
}

// Build interval channel (idx, interval_spec) from reference lengths using the requested interval size.
def build_interval_channel() {
    def intervals = []
    long blockSize = parse_interval_size_to_bases(params.interval_size)
    int idx = 1
    load_reference_lengths().each { entry ->
        def contig = entry.name
        long contigLength = entry.length
        long start = 1L
        while (start <= contigLength) {
            long end = Math.min(start + blockSize - 1L, contigLength)
            def intervalIdx = String.format('%04d', idx)
            def intervalSpec = "${contig}:${start}-${end}"
            intervals << tuple(intervalIdx, intervalSpec)
            start = end + 1L
            idx++
        }
    }

    if (!intervals) {
        error "No calling intervals could be generated from reference: ${params.ref_abs}"
    }

    return Channel.from(intervals)
}

// Step3 dispatches clean reads into the alignment subworkflow.
def run_step3(ch_clean_reads) {
    def ch_ref_fa = Channel.value(file(params.ref_abs))
    return ALIGNMENT(ch_clean_reads, ch_ref_fa)
}

// Step4 dispatches markdup alignments into the variant-calling subworkflow.
def run_step4(ch_markdup) {
    def ch_ref_fa = Channel.value(file(params.ref_abs))
    def ch_intervals = build_interval_channel()
    return VARIANT_CALLING(ch_markdup, ch_ref_fa, ch_intervals).final_vcf
}

def run_step5_vcffilter(ch_called_vcf) {
    if (!params.enable_vcffilter) {
        return ch_called_vcf
    }
    return VCF_FILTERING(ch_called_vcf).final_vcf
}

def build_called_vcf_from_input_vcf() {
    def inputVcf = file(params.input_vcf)
    def inputTbi = file("${params.input_vcf}.tbi")
    def callsetId = inputVcf.getName()
        .replaceFirst(/\.vcf\.gz$/, '')
        .replaceFirst(/\.vcf\.bgz$/, '')
        .replaceFirst(/\.vcf$/, '')

    return Channel.value(tuple(callsetId, inputVcf, inputTbi))
}

def run_trimming_stage(step1_done_signal = null, ch_fastp_input_override = null) {
    ch_fastp_input = ch_fastp_input_override ?: build_fastp_input_channel()
    if (step1_done_signal != null) {
        ch_fastp_input = ch_fastp_input.combine(step1_done_signal).map { row ->
            def data = row.take(row.size() - 1)
            tuple(*data)
        }
    }
    return TRIMMING_QC(ch_fastp_input)
}

workflow {
    if (params.help_full) {
        print_help(true)
        exit 0
    }

    if (params.help) {
        print_help(false)
        exit 0
    }

    validate_params()

    def step1_done_signal = null
    def trim = null
    def markdup = null
    def raw_fastq_source = null
    def fastp_input_source = null

    if (params.run_step in ['Raw_QC', 'all']) {
        if (params.SRA_list) {
            sra_download = DOWNLOAD_SRA(build_sra_sample_channel())
            raw_fastq_source = sra_download.raw_fastq
            fastp_input_source = sra_download.fastp_input
        } else {
            raw_fastq_source = build_raw_fastq_channel()
        }

        raw_qc = RAW_QC(raw_fastq_source)
        step1_done_signal = raw_qc.multiqc_report
    }

    if (params.run_step == 'all') {
        trim = run_trimming_stage(step1_done_signal, fastp_input_source)
        markdup = run_step3(trim.cleaned_reads)
        run_step5_vcffilter(run_step4(markdup.markdup_alignment))
    } else if (params.run_step == 'Raw_QC') {
        log.info "run_step=Raw_QC: finished raw QC only; downstream steps are skipped by design"
    } else if (params.run_step == 'Trimming_QC') {
        trim = run_trimming_stage()
        markdup = run_step3(trim.cleaned_reads)
        run_step5_vcffilter(run_step4(markdup.markdup_alignment))
    } else if (params.run_step == 'Aligning') {
        def clean_reads_for_align
        if (params.trimmed_input_dir) {
            log.info "run_step=Aligning: using external trimmed reads from --trimmed_input_dir ${params.trimmed_input_dir}"
            clean_reads_for_align = build_clean_reads_from_trimmed_input_dir()
        } else if (has_complete_step2_outputs()) {
            log.info "run_step=Aligning: reusing Step2 outputs from results/02_fastp_qc/fastp"
            clean_reads_for_align = build_clean_reads_from_step2_results()
        } else {
            log.info "run_step=Aligning: Step2 outputs incomplete, auto-running Trimming_QC once"
            trim = run_trimming_stage()
            clean_reads_for_align = trim.cleaned_reads
        }

        markdup = run_step3(clean_reads_for_align)
        run_step5_vcffilter(run_step4(markdup.markdup_alignment))
    } else if (params.run_step == 'Calling') {
        def markdup_for_calling
        if (has_complete_step3_outputs()) {
            log.info "run_step=Calling: reusing Step3 outputs from results/03_markdup"
            markdup_for_calling = build_markdup_from_step3_results()
        } else {
            log.info "run_step=Calling: Step3 outputs incomplete, auto-running Aligning once"
            if (params.trimmed_input_dir) {
                log.info "run_step=Calling: using external trimmed reads from --trimmed_input_dir ${params.trimmed_input_dir} for one-step fallback to Aligning"
                markdup = run_step3(build_clean_reads_from_trimmed_input_dir())
            } else if (!has_complete_step2_outputs()) {
                error "run_step=Calling can auto-fallback only one step to Aligning, but Step2 outputs are incomplete under results/02_fastp_qc/fastp"
            } else {
                markdup = run_step3(build_clean_reads_from_step2_results())
            }
            markdup_for_calling = markdup.markdup_alignment
        }
        run_step5_vcffilter(run_step4(markdup_for_calling))
    } else if (params.run_step == 'VCF_Filtering') {
        log.info "run_step=VCF_Filtering: using external VCF from --input_vcf ${params.input_vcf}"
        VCF_FILTERING(build_called_vcf_from_input_vcf())
    }
}
