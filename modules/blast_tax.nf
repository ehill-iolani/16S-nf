process BLAST_TAX {
    tag "$sample $chunk"
    // process_medium, not process_low: blastn scales with threads, and a
    // several-thousand-query soil sample runs for hours -- the time limit is
    // only there to catch a hung task, since a timeout restarts from zero
    label 'process_medium'
    container "${params.container_registry}/biocontainers/blast:2.15.0--pl5321h6f7f691_1"

    input:
    tuple val(sample), val(chunk), path(consensus_fastas)
    path db_files
    val db_name

    output:
    tuple val(sample), path("${sample}.chunk${chunk}.hits.tsv"), emit: hits

    script:
    // one blastn over a chunk of a sample's consensus sequences (the workflow
    // splits them by --blast_chunk_size), not one task per cluster -- each
    // task stages the whole BLAST db (hundreds of MB for SILVA), so a single
    // sequence per task would spend far longer staging than searching. The
    // file is sample- and chunk-prefixed since these all get collected
    // together in BUILD_REPORT, which joins each sample's chunks back into
    // one published <sample>.hits.tsv.
    //
    // `awk 1` rather than cat so a fasta missing its trailing newline can't
    // run its last line into the next file's header
    """
    awk 1 ${consensus_fastas} > queries.fasta

    blastn -query queries.fasta -db ${db_name} \\
        -num_threads ${task.cpus} \\
        -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \\
        -max_target_seqs 5 -evalue 1e-10 \\
        -out blast.tsv

    # unmatched clusters (no hit) are flagged, not dropped -- a community
    # sample can easily contain strains that are simply absent from whatever
    # reference set was supplied, and that's a real (if uncertain) result,
    # not a pipeline failure
    #
    # qseqid here must match the id blastn used (the consensus fasta's header
    # token, "sample_clusterid") -- otherwise this row can't be joined back to
    # its sample in BUILD_REPORT. The FILENAME test (not NR == FNR) keeps this
    # right when blast.tsv is empty, i.e. when nothing at all hit.
    awk -F'\\t' '
        FILENAME == ARGV[1] { seen[\$1] = 1; next }
        /^>/ {
            id = substr(\$0, 2)
            sub(/[ \t].*/, "", id)
            if (!(id in seen)) printf "%s\\tNO_HIT\\tNA\\tNA\\tNA\\tNA\\tflag_for_manual_review\\n", id
        }
    ' blast.tsv queries.fasta > no_hit.tsv

    cat blast.tsv no_hit.tsv > ${sample}.chunk${chunk}.hits.tsv
    """
}
