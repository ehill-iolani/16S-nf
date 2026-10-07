process MERGE_CONSENSUS {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/vsearch:2.30.6--h0bb26bb_0"
    publishDir(path: { "${params.outdir}/${sample}/06_merged" }, mode: 'copy')

    input:
    tuple val(sample), val(cluster_ids), path(consensus_fastas)

    output:
    tuple val(sample), path("merged/*.merged.fasta"), emit: merged
    tuple val(sample), path("merge_map.tsv"), emit: merge_map

    script:
    /*
     * Second vsearch pass, on the per-cluster consensus sequences rather than
     * raw reads. The read-level pass (VSEARCH_CLUSTER) has to run loose enough
     * to cope with raw ONT error, and a loose-but-noisy centroid still splits
     * one species into several clusters. Consensus sequences are far more
     * accurate than reads, so a tight identity here (merge_id) cleanly folds
     * those fragments back together without merging distinct species.
     *
     * Each group keeps its largest cluster's consensus as the representative
     * and sums the members' cluster_size into it, so downstream steps
     * (BLAST_TAX / SORT_CONSENSUS / BUILD_REPORT) see one fasta per group,
     * stamped exactly like an unmerged consensus.
     */
    """
    set -o pipefail

    # vsearch weights centroids by a ;size=N label annotation: turn the
    # stamped cluster_size=N into one (largest cluster becomes the centroid).
    # awk reads the files itself (not cat) so a file missing its trailing
    # newline can't run its last line into the next file's header
    awk '
        /^>/ {
            n = 1
            for (i = 2; i <= NF; i++) if (\$i ~ /^cluster_size=/) n = substr(\$i, 14)
            print ">" substr(\$1, 2) ";size=" n
            next
        }
        { print }
    ' ${consensus_fastas} > consensus.fasta

    # --strand both: clusters were oriented to their own centroid, so two
    # fragments of the same species can still sit on opposite strands
    vsearch --cluster_size consensus.fasta \\
        --id ${params.merge_id} \\
        --strand both \\
        --sizein --sizeout \\
        --fasta_width 0 \\
        --threads ${task.cpus} \\
        --centroids centroids.fasta \\
        --uc merge.uc

    # one fasta per merged group, restamped in the header format the rest of
    # the pipeline reads (>{sample}_{cluster_id} sample={sample} cluster_size={n}),
    # with cluster_size now the summed read count. sample-prefixed filenames
    # because these get collected across samples in BUILD_REPORT
    mkdir -p merged
    awk -v s="${sample}" '
        /^>/ {
            if (out != "") close(out)
            split(substr(\$0, 2), a, ";")
            id = a[1]
            size = a[2]
            sub(/size=/, "", size)
            out = "merged/" s "." substr(id, length(s) + 2) ".merged.fasta"
            print ">" id " sample=" s " cluster_size=" size > out
            next
        }
        { print > out }
    ' centroids.fasta

    # which clusters were folded into which representative (S = centroid itself)
    printf 'member\\trepresentative\\n' > merge_map.tsv
    awk -F'\\t' '\$1 == "S" || \$1 == "H" {
        q = \$9; sub(/;.*/, "", q)
        t = (\$1 == "S") ? q : \$10; sub(/;.*/, "", t)
        print q "\\t" t
    }' merge.uc >> merge_map.tsv

    echo "merged \$(grep -c '^>' consensus.fasta) cluster consensus sequences into \$(ls merged/*.merged.fasta | wc -l) at --merge_id ${params.merge_id}" >&2
    """
}
