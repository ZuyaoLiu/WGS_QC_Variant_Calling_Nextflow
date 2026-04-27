process PHASE1_INITIAL_FILTER_MASK_PROCESS {
    tag "${callset_id}"
    container "${params.container_image}"
    publishDir 'results/05_vcf_filtering/01_phase1/intermediate', mode: 'copy'

    input:
    tuple val(callset_id), path(input_vcf), path(input_index)

    output:
    tuple val(callset_id), path("${callset_id}.phase1.pre_dp_mask.bcf"), path("${callset_id}.phase1.pre_dp_mask.bcf.csi"), emit: filtered_vcf

    script:
    def mislabeledFile = params.phase1_mislabeled_samples ? file(params.phase1_mislabeled_samples) : null
    def mislabeledArg = mislabeledFile ? mislabeledFile.toAbsolutePath().toString() : ''
    """
    bcftools view \
      --threads ${task.cpus} \
      -m2 -M2 -v snps \
      -Ou \
      ${input_vcf} \
    | {
        if [[ "${params.phase1_remove_mislabeled}" == "true" ]]; then
          if [[ -z "${mislabeledArg}" ]]; then
            echo "phase1_remove_mislabeled=true but phase1_mislabeled_samples is not set" >&2
            exit 1
          fi
          bcftools view \
            --threads ${task.cpus} \
            -S ^${mislabeledArg} \
            --force-samples \
            -Ou
        else
          bcftools view \
            --threads ${task.cpus} \
            -Ou
        fi
      } \
    | bcftools +setGT \
      --threads ${task.cpus} \
      -Ou \
      -- \
      -t q -n . \
      -i "FMT/GQ<${params.phase1_min_gq} | FMT/DP<${params.phase1_min_dp}" \
    | {
        if [[ "${params.phase1_enable_ab_mask}" == "true" ]]; then
          bcftools +setGT \
            --threads ${task.cpus} \
            -Ou \
            -- \
            -t q -n . \
            -i 'GT="het" & (FMT/AD[:0]=0 | FMT/AD[:1]=0)' \
          | bcftools +setGT \
            --threads ${task.cpus} \
            -Ou \
            -- \
            -t "b:AD<0.005" -n X \
          | bcftools +setGT \
            --threads ${task.cpus} \
            -Ob \
            -o ${callset_id}.phase1.pre_dp_mask.bcf \
            -- \
            -t "b:AD<0.01" -n .
        else
          bcftools view \
            --threads ${task.cpus} \
            -Ob \
            -o ${callset_id}.phase1.pre_dp_mask.bcf
        fi
      }

    bcftools index --threads ${task.cpus} -f ${callset_id}.phase1.pre_dp_mask.bcf
    """
}

workflow PHASE1_INITIAL_FILTER_MASK {
    take:
    ch_called_vcf

    main:
    PHASE1_INITIAL_FILTER_MASK_PROCESS(ch_called_vcf)

    emit:
    filtered_vcf = PHASE1_INITIAL_FILTER_MASK_PROCESS.out.filtered_vcf.map { callset_id, vcf, csi ->
        tuple(
            callset_id,
            file("results/05_vcf_filtering/01_phase1/intermediate/${vcf.name}"),
            file("results/05_vcf_filtering/01_phase1/intermediate/${csi.name}")
        )
    }
}
