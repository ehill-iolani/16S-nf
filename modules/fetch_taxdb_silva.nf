process FETCH_TAXDB_SILVA {
    label 'process_low'
    // needs /bin/bash to exist (the Nextflow docker executor execs it
    // directly, regardless of the `shell` directive) -- ubuntu has it plus
    // gunzip/awk already installed, only curl needs adding
    container "${params.dockerhub_registry}/library/ubuntu:22.04"
    storeDir params.tax_db_cache

    output:
    path "silva_138_16s.fasta", emit: fasta
    // SILVA's rank for every taxonomy path, which BUILD_REPORT needs to name
    // each lineage level (the fasta's lineages skip ranks now and then). A
    // cache holding only the fasta, from before this was an output, is
    // re-fetched once, since storeDir only skips a task when all its outputs exist
    path "silva_138_taxmap.txt", emit: taxmap

    script:
    // SILVA's SSURef export spans all three domains (small-subunit rRNA in
    // general, not 16S-specific) and stores sequences as RNA (U); this is a
    // bacterial/archaeal 16S pipeline, so restrict to Bacteria/Archaea and
    // convert U->T, otherwise blastn -- which scores in the DNA alphabet --
    // would count every U as a mismatch against the (DNA) consensus queries
    """
    apt-get update -qq && apt-get install -y -qq --no-install-recommends curl ca-certificates

    curl -fsSL --retry 3 -o silva.fasta.gz \\
        "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"

    gunzip -c silva.fasta.gz \\
        | awk '
            /^>/ { keep = (\$0 ~ /Bacteria;/ || \$0 ~ /Archaea;/); if (keep) print; next }
            keep { gsub(/U/, "T"); gsub(/u/, "t"); print }
          ' \\
        > silva_138_16s.fasta
    rm -f silva.fasta.gz

    if [ ! -s silva_138_16s.fasta ]; then
        echo "ERROR: filtered SILVA reference is empty -- check the Bacteria/Archaea filter against the current SILVA release's header format" >&2
        exit 1
    fi

    curl -fsSL --retry 3 -o taxmap.txt.gz \\
        "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz"
    gunzip -c taxmap.txt.gz > silva_138_taxmap.txt
    rm -f taxmap.txt.gz

    if [ ! -s silva_138_taxmap.txt ]; then
        echo "ERROR: downloaded SILVA taxonomy rank map is empty" >&2
        exit 1
    fi
    """
}
