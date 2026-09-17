process MEDAKA {
    tag "${sample}:${cluster_id}"
    label 'process_medium'
    // ONT's own multi-arch image, not biocontainers: the biocontainers medaka
    // build is amd64-only and its TensorFlow backend SIGILLs under Docker's
    // amd64 emulation on Apple Silicon; this one has a native arm64 build too
    container 'ontresearch/medaka:v1.11.3'
    publishDir(path: { "${params.outdir}/${sample}/06_consensus" }, mode: 'copy')

    input:
    tuple val(sample), val(cluster_id), path(racon_fasta), path(cluster_fastq)

    output:
    tuple val(sample), val(cluster_id), path("${sample}.${cluster_id}.medaka.consensus.fasta"), emit: consensus

    script:
    // sample-prefixed for the same reason as RACON's output -- cluster ids
    // repeat across samples and these get collected together in BUILD_REPORT
    """
    medaka_consensus -i ${cluster_fastq} -d ${racon_fasta} -o medaka_out -t ${task.cpus}

    # stamp the sample/cluster-size header on, same as RACON -- medaka renames
    # the record based on the racon draft's header rather than keeping it as-is
    n_reads=\$(( \$(wc -l < ${cluster_fastq}) / 4 ))
    awk -v s="${sample}" -v c="${cluster_id}" -v n="\$n_reads" \\
        'NR==1 { print ">" s "_" c " sample=" s " cluster_size=" n; next } { print }' \\
        medaka_out/consensus.fasta > ${sample}.${cluster_id}.medaka.consensus.fasta
    """
}
