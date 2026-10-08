process SPOA_CONSENSUS {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/spoa:4.1.4--h077b44d_3"
    publishDir(path: { "${params.outdir}/${sample}/04_draft" }, mode: 'copy')

    input:
    tuple val(sample), path(cluster_fastqs)

    output:
    tuple val(sample), path("*.draft.fasta"), emit: draft

    script:
    // one task per sample, looping over its clusters (spoa is single-threaded,
    // so clusters run side by side, one per cpu) -- a task per cluster spends
    // far longer on VM start-up and staging than on the consensus itself
    //
    // the spoa container ships only spoa itself, no seqtk -- convert fastq to
    // fasta with a plain awk one-liner instead of pulling in another tool/container
    """
    . par_each.sh

    draft_one() {
        cid=\$(basename "\$1" .fastq)
        awk 'NR % 4 == 1 { print ">" substr(\$0, 2) } NR % 4 == 2 { print }' "\$1" > "\$cid.reads.fasta"
        spoa "\$cid.reads.fasta" -r 0 > "\$cid.draft.fasta"
        rm "\$cid.reads.fasta"
    }

    par_each ${task.cpus} draft_one *.fastq
    """
}
