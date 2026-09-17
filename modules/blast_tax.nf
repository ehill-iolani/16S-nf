process BLAST_TAX {
    tag "${sample}:${cluster_id}"
    label 'process_low'
    container 'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'
    publishDir(path: { "${params.outdir}/${sample}/07_taxonomy/${cluster_id}" }, mode: 'copy')

    input:
    tuple val(sample), val(cluster_id), path(consensus)
    path db_files
    val db_name

    output:
    tuple val(sample), val(cluster_id), path("${sample}.${cluster_id}.hits.tsv"), emit: hits

    script:
    // sample-prefixed for the same reason as RACON's output -- cluster ids
    // repeat across samples and these get collected together in BUILD_REPORT
    """
    blastn -query ${consensus} -db ${db_name} \\
        -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \\
        -max_target_seqs 5 -evalue 1e-10 \\
        -out ${sample}.${cluster_id}.hits.tsv

    # unmatched clusters (no hit) are flagged, not dropped -- a community
    # sample can easily contain strains that are simply absent from whatever
    # reference set was supplied, and that's a real (if uncertain) result,
    # not a pipeline failure
    #
    # qseqid here must match the id blastn would have used (the consensus
    # fasta's header token, "sample_clusterid"), not the bare cluster_id --
    # otherwise this row can't be joined back to its sample in BUILD_REPORT
    if [ ! -s ${sample}.${cluster_id}.hits.tsv ]; then
        echo -e "${sample}_${cluster_id}\\tNO_HIT\\tNA\\tNA\\tNA\\tNA\\tflag_for_manual_review" > ${sample}.${cluster_id}.hits.tsv
    fi
    """
}
