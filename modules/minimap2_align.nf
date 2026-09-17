process MINIMAP2_ALIGN {
    tag "${sample}:${cluster_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/minimap2:2.28--he4a0461_3'
    publishDir(path: { "${params.outdir}/${sample}/05_racon/${cluster_id}" }, mode: 'copy')

    input:
    tuple val(sample), val(cluster_id), path(draft), path(cluster_fastq)

    output:
    tuple val(sample), val(cluster_id), path(draft), path(cluster_fastq), path("${cluster_id}.sam"), emit: aligned

    script:
    """
    minimap2 -ax map-ont ${draft} ${cluster_fastq} > ${cluster_id}.sam
    """
}
