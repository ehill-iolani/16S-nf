process RACON {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/racon:1.5.0--h21ec9f0_2"
    publishDir(path: { "${params.outdir}/${sample}/05_racon" }, mode: 'copy')

    input:
    tuple val(sample), path(drafts), path(cluster_fastqs), path(sams)

    output:
    // optional: a sample whose every cluster was skipped (below) has no
    // polished fasta, and drops out of the steps after this one
    tuple val(sample), path("*.racon.fasta"), path(cluster_fastqs), optional: true, emit: polished
    path("${sample}.racon_skipped.tsv"), optional: true, emit: skipped

    script:
    // filenames are sample-prefixed (not just cluster_id) because cluster ids
    // are only unique within a sample, and these files later get collected
    // across all samples into one BUILD_REPORT call -- a bare "0.racon.fasta"
    // would collide with every other sample's cluster 0
    //
    // the racon biocontainer ships only racon itself, no minimap2 -- alignment
    // happens upstream in MINIMAP2_ALIGN. One task per sample, one single-
    // threaded racon per cluster, side by side.
    //
    // A cluster none of whose reads aligned to its own draft is skipped, not
    // polished: racon exits 1 on it ("empty overlap set!"), which used to end
    // the whole run (found in validation, 2026-11-12: a cluster of 33 bp reads
    // in a 2 kb amplicon run, too short for minimap2 to place). A draft its own
    // reads don't support isn't a consensus worth identifying, so it isn't
    // passed on unpolished either. The cluster is recorded in
    // <sample>.racon_skipped.tsv (sample, cluster, reads, draft_length,
    // reason), which BUILD_REPORT gathers into final_report/skipped_clusters.tsv.
    """
    . par_each.sh

    polish_one() {
        local cid n_reads draft_len
        cid=\$(basename "\$1" .fastq)
        n_reads=\$(( \$(wc -l < "\$1") / 4 ))

        # RNAME is "*" for an unmapped read; any other value is an overlap
        if ! awk '!/^@/ && \$3 != "*" { found = 1; exit } END { exit !found }' "\$cid.sam"; then
            draft_len=\$(awk '!/^>/ { n += length(\$0) } END { print n + 0 }' "\$cid.draft.fasta")
            echo "RACON: no read of ${sample} cluster \$cid (\$n_reads reads, draft \$draft_len bp) aligned to its draft; skipping the cluster" >&2
            printf '%s\\t%s\\t%s\\t%s\\tno_reads_aligned\\n' "${sample}" "\$cid" "\$n_reads" "\$draft_len" > "\$cid.skipped.row"
            return 0
        fi

        racon -t 1 "\$1" "\$cid.sam" "\$cid.draft.fasta" > "\$cid.racon.raw"

        # racon always names the record "Consensus" and regenerates its own
        # description (LN:/RC:/XC: tags), discarding anything set upstream -- so
        # the sample/cluster-size header has to be stamped on here, after racon runs
        awk -v s="${sample}" -v c="\$cid" -v n="\$n_reads" \\
            'NR==1 { print ">" s "_" c " sample=" s " cluster_size=" n; next } { print }' \\
            "\$cid.racon.raw" > "${sample}.\$cid.racon.fasta"
        rm "\$cid.racon.raw"
    }

    par_each ${task.cpus} polish_one *.fastq

    # one row per skipped cluster: each worker writes its own file (they run
    # side by side, so they don't append to one), joined here
    if ls *.skipped.row >/dev/null 2>&1; then
        cat *.skipped.row > ${sample}.racon_skipped.tsv
    fi
    """
}
