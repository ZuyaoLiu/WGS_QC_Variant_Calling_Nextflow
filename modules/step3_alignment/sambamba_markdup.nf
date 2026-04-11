process SAMBAMBA_MARKDUP_PROCESS {
    tag "${sample_id}"
    container "${params.sif}"
    cpus { (params.sambamba_cpus ?: params.threads ?: 1) as Integer }
    publishDir 'results/03_markdup', mode: 'move'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), path("${sample_id}.markdup.bam.bai"), emit: markdup

    script:
    """
    sambamba markdup -t ${task.cpus} --tmpdir . ${sorted_bam} ${sample_id}.markdup.bam
    sambamba index -t ${task.cpus} ${sample_id}.markdup.bam
    """
}

workflow SAMBAMBA_MARKDUP {
    take:
    ch_input

    main:
    SAMBAMBA_MARKDUP_PROCESS(ch_input)

    emit:
    markdup_bam = SAMBAMBA_MARKDUP_PROCESS.out.markdup.map { sample_id, bam, bai ->
        tuple(
            sample_id,
            file("results/03_markdup/${bam.name}"),
            file("results/03_markdup/${bai.name}")
        )
    }
}
