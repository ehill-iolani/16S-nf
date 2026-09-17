process SORT_CONSENSUS {
    tag "${sample}:${cluster_id}"
    label 'process_low'
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir "${params.outdir}/consensus_by_confidence", mode: 'copy'

    input:
    tuple val(sample), val(cluster_id), path(consensus), path(hits)

    output:
    path("*/${sample}.${cluster_id}.consensus.fasta")

    script:
    // blastn's outfmt 6 output is sorted best-hit-first per query by default,
    // so the top line of ${hits} is the same "best hit" build_report.py picks
    """
    subject=\$(awk -F'\\t' 'NR==1{print \$2}' ${hits})
    pident=\$(awk -F'\\t' 'NR==1{print \$3}' ${hits})

    if [ "\$subject" = "NO_HIT" ]; then
        category="no_hit"
    elif awk -v p="\$pident" 'BEGIN{exit !(p < ${params.min_pident})}'; then
        category="low_confidence"
    else
        category="confident"
    fi

    mkdir -p "\$category"
    cp ${consensus} "\$category/${sample}.${cluster_id}.consensus.fasta"
    """
}
