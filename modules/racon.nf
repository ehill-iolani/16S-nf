process RACON {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/racon:1.5.0--h21ec9f0_2"
    publishDir(path: { "${params.outdir}/${sample}/05_racon" }, mode: 'copy')

    input:
    tuple val(sample), path(drafts), path(cluster_fastqs), path(sams)

    output:
    tuple val(sample), path("*.racon.fasta"), path(cluster_fastqs), emit: polished

    script:
    // filenames are sample-prefixed (not just cluster_id) because cluster ids
    // are only unique within a sample, and these files later get collected
    // across all samples into one BUILD_REPORT call -- a bare "0.racon.fasta"
    // would collide with every other sample's cluster 0
    //
    // the racon biocontainer ships only racon itself, no minimap2 -- alignment
    // happens upstream in MINIMAP2_ALIGN. One task per sample, one single-
    // threaded racon per cluster, side by side.
    """
    . par_each.sh

    polish_one() {
        cid=\$(basename "\$1" .fastq)
        racon -t 1 "\$1" "\$cid.sam" "\$cid.draft.fasta" > "\$cid.racon.raw"

        # racon always names the record "Consensus" and regenerates its own
        # description (LN:/RC:/XC: tags), discarding anything set upstream -- so
        # the sample/cluster-size header has to be stamped on here, after racon runs
        n_reads=\$(( \$(wc -l < "\$1") / 4 ))
        awk -v s="${sample}" -v c="\$cid" -v n="\$n_reads" \\
            'NR==1 { print ">" s "_" c " sample=" s " cluster_size=" n; next } { print }' \\
            "\$cid.racon.raw" > "${sample}.\$cid.racon.fasta"
        rm "\$cid.racon.raw"
    }

    par_each ${task.cpus} polish_one *.fastq
    """
}
