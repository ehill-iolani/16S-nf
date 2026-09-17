process RACON {
    tag "${sample}:${cluster_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/racon:1.5.0--h21ec9f0_2'
    publishDir(path: { "${params.outdir}/${sample}/05_racon/${cluster_id}" }, mode: 'copy')

    input:
    tuple val(sample), val(cluster_id), path(draft), path(cluster_fastq), path(sam)

    output:
    tuple val(sample), val(cluster_id), path("${sample}.${cluster_id}.racon.fasta"), path(cluster_fastq), emit: polished

    script:
    // filenames are sample-prefixed (not just cluster_id) because cluster ids
    // are only unique within a sample, and these files later get collected
    // across all samples into one BUILD_REPORT call -- a bare "0.racon.fasta"
    // would collide with every other sample's cluster 0
    //
    // the racon biocontainer ships only racon itself, no minimap2 -- alignment
    // happens upstream in MINIMAP2_ALIGN
    """
    racon ${cluster_fastq} ${sam} ${draft} > raw.racon.fasta

    # racon always names the record "Consensus" and regenerates its own
    # description (LN:/RC:/XC: tags), discarding anything set upstream -- so
    # the sample/cluster-size header has to be stamped on here, after racon runs
    n_reads=\$(( \$(wc -l < ${cluster_fastq}) / 4 ))
    awk -v s="${sample}" -v c="${cluster_id}" -v n="\$n_reads" \\
        'NR==1 { print ">" s "_" c " sample=" s " cluster_size=" n; next } { print }' \\
        raw.racon.fasta > ${sample}.${cluster_id}.racon.fasta
    """
}
