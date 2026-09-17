process FETCH_TAXDB_REFSEQ {
    label 'process_low'
    // needs /bin/bash to exist (the Nextflow docker executor execs it
    // directly, regardless of the `shell` directive) -- ubuntu has it plus
    // gunzip/awk already installed, only curl needs adding
    container 'ubuntu:22.04'
    storeDir params.tax_db_cache

    output:
    path "refseq_16s.fasta", emit: fasta

    script:
    // NCBI's Targeted Loci 16S export is split by domain -- pull both
    // Bacteria and Archaea so this matches SILVA's domain coverage above;
    // already plain DNA fasta, no alphabet conversion needed
    """
    apt-get update -qq && apt-get install -y -qq --no-install-recommends curl ca-certificates

    curl -fsSL --retry 3 -o bacteria.fna.gz \\
        "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
    curl -fsSL --retry 3 -o archaea.fna.gz \\
        "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz"

    gunzip -c bacteria.fna.gz archaea.fna.gz > refseq_16s.fasta
    rm -f bacteria.fna.gz archaea.fna.gz

    if [ ! -s refseq_16s.fasta ]; then
        echo "ERROR: downloaded RefSeq 16S reference is empty" >&2
        exit 1
    fi
    """
}
