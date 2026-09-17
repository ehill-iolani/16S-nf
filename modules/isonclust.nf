process ISONCLUST {
    tag "$sample"
    label 'process_medium'
    container 'quay.io/biocontainers/isonclust:0.0.6.1--py_0'
    publishDir(path: { "${params.outdir}/${sample}/03_clusters" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("clusters/*.fastq"), emit: clusters
    tuple val(sample), path("isonclust_out/final_clusters.tsv"), emit: cluster_tsv

    script:
    /*
     * NOTE: cluster_id here maps to isONclust's internal similarity/quality
     * thresholds via --min_shared / --mapped_threshold. This is the single
     * most sensitivity-critical parameter in the whole pipeline: too loose
     * merges distinct strains/species into one cluster, too strict splits
     * one organism's reads (or its multiple rRNA operon copies, which can
     * differ slightly within a single genome) into several clusters --
     * re-tune per primer set / amplicon region (full-length 16S vs. a single
     * V-region behave differently here).
     */
    """
    # isONclust reads its input as plain text and chokes on gzip -- decompress first.
    # Modern basecaller (Dorado) headers carry tab-separated tags after the read id
    # (qs:f:.. mx:i:.. etc); isONclust's clustering step and its write_fastq step
    # parse that header differently, which produces a stale accession -> KeyError
    # inside write_fastq. Strip everything but the bare read id to sidestep it.
    zcat -f ${fastq} | awk 'NR % 4 == 1 { sub(/[ \\t].*/, "", \$0) } { print }' > input.fastq

    isONclust \\
        --fastq input.fastq \\
        --outfolder isonclust_out \\
        --mapped_threshold ${params.cluster_id}

    mkdir -p clusters
    isONclust write_fastq \\
        --clusters isonclust_out/final_clusters.tsv \\
        --fastq input.fastq \\
        --outfolder clusters \\
        --N ${params.min_cluster}
    """
}
