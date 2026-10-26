process SPOA_CONSENSUS {
    tag "$sample"
    // process_medium, not process_low: more cpus for samples with thousands of
    // clusters, and a time limit that only catches a hung task -- a timeout
    // restarts the task from zero, so a too-tight limit just burns compute
    label 'process_medium'
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
    //
    // spoa only builds the draft, so it gets at most --spoa_max_reads reads,
    // spread evenly through the cluster file (which holds the plus-strand reads,
    // then the flipped minus-strand ones). Its run time grows much faster than
    // the read count: a 143k-read cluster took over 3 h, where a few hundred
    // reads already give a draft as good as racon needs. RACON still polishes
    // with every read, and cluster_size still counts every read.
    """
    . par_each.sh

    draft_one() {
        cid=\$(basename "\$1" .fastq)
        n=\$(( \$(wc -l < "\$1") / 4 ))
        awk -v n="\$n" -v cap=${params.spoa_max_reads} '
            NR % 4 == 1 {
                i = (NR - 1) / 4
                keep = (cap <= 0 || n <= cap || int((i + 1) * cap / n) > int(i * cap / n))
                if (keep) print ">" substr(\$0, 2)
            }
            NR % 4 == 2 && keep { print }
        ' "\$1" > "\$cid.reads.fasta"
        spoa "\$cid.reads.fasta" -r 0 > "\$cid.draft.fasta"
        rm "\$cid.reads.fasta"
    }

    par_each ${task.cpus} draft_one *.fastq
    """
}
