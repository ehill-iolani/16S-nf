process MEDAKA {
    tag "$sample"
    label 'process_medium'
    // ONT's own multi-arch image, not biocontainers: the biocontainers medaka
    // build is amd64-only and its TensorFlow backend SIGILLs under Docker's
    // amd64 emulation on Apple Silicon; this one has a native arm64 build too
    container "${params.dockerhub_registry}/ontresearch/medaka:v1.11.3"
    publishDir(path: { "${params.outdir}/${sample}/06_consensus" }, mode: 'copy')

    input:
    tuple val(sample), path(racon_fastas), path(cluster_fastqs)

    output:
    tuple val(sample), path("*.medaka.consensus.fasta"), emit: consensus

    script:
    // one task per sample; clusters go one after another with all the task's
    // cpus on each medaka call, since medaka already multi-threads internally
    // (and its model memory makes running several at once a poor trade)
    //
    // sample-prefixed for the same reason as RACON's output -- cluster ids
    // repeat across samples and these get collected together in BUILD_REPORT
    """
    for fq in *.fastq; do
        cid=\$(basename "\$fq" .fastq)
        medaka_consensus -i "\$fq" -d "${sample}.\$cid.racon.fasta" -o "medaka_\$cid" -t ${task.cpus}

        # stamp the sample/cluster-size header on, same as RACON -- medaka renames
        # the record based on the racon draft's header rather than keeping it as-is
        n_reads=\$(( \$(wc -l < "\$fq") / 4 ))
        awk -v s="${sample}" -v c="\$cid" -v n="\$n_reads" \\
            'NR==1 { print ">" s "_" c " sample=" s " cluster_size=" n; next } { print }' \\
            "medaka_\$cid/consensus.fasta" > "${sample}.\$cid.medaka.consensus.fasta"
        rm -r "medaka_\$cid"
    done
    """
}
