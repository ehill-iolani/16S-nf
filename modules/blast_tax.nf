process BLAST_TAX {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/blast:2.15.0--pl5321h6f7f691_1"
    publishDir(path: { "${params.outdir}/${sample}/07_taxonomy" }, mode: 'copy')

    input:
    tuple val(sample), path(consensus_fastas)
    path db_files
    val db_name

    output:
    tuple val(sample), path("${sample}.hits.tsv"), emit: hits

    script:
    // one blastn over all of a sample's consensus sequences, not one task per
    // cluster -- each task would otherwise stage the whole BLAST db (hundreds
    // of MB for SILVA) just to search a single sequence. File is sample-
    // prefixed since these get collected together in BUILD_REPORT.
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

    cat blast.tsv no_hit.tsv > ${sample}.hits.tsv
    """
}
