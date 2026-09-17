process SPOA_CONSENSUS {
    tag "${sample}:${cluster_id}"
    label 'process_low'
    container 'quay.io/biocontainers/spoa:4.1.4--h077b44d_3'
    publishDir(path: { "${params.outdir}/${sample}/04_draft/${cluster_id}" }, mode: 'copy')

    input:
    tuple val(sample), val(cluster_id), path(cluster_fastq)

    output:
    tuple val(sample), val(cluster_id), path("${cluster_id}.draft.fasta"), emit: draft

    script:
    // the spoa container ships only spoa itself, no seqtk -- convert fastq to
    // fasta with a plain awk one-liner instead of pulling in another tool/container
    """
    awk 'NR % 4 == 1 { print ">" substr(\$0, 2) } NR % 4 == 2 { print }' ${cluster_fastq} > ${cluster_id}.reads.fasta
    spoa ${cluster_id}.reads.fasta -r 0 > ${cluster_id}.draft.fasta
    """
}
