process MINIMAP2_ALIGN {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/minimap2:2.28--he4a0461_3"

    input:
    tuple val(sample), path(drafts), path(cluster_fastqs)

    output:
    tuple val(sample), path(drafts), path(cluster_fastqs), path("*.sam"), emit: aligned

    script:
    // per sample, one single-threaded minimap2 per cluster, side by side.
    // The SAMs are only an intermediate for RACON, so they aren't published.
    """
    . par_each.sh

    align_one() {
        cid=\$(basename "\$1" .fastq)
        minimap2 -ax map-ont -t 1 "\$cid.draft.fasta" "\$1" > "\$cid.sam"
    }

    par_each ${task.cpus} align_one *.fastq
    """
}
