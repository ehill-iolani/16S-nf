process PCOA {
    label 'process_low'
    // Same image BUILD_REPORT already uses -- proven to work on this exact
    // Batch setup, and pandas' own dependency on numpy is all pcoa.py needs
    // (classical PCoA is implemented directly with numpy, not scipy/
    // scikit-bio, to avoid pulling in a second, unverified container).
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir "${params.outdir}/final_report", mode: 'copy'

    input:
    path abundance_table
    path metadata_csv

    output:
    path "pcoa_coordinates.tsv", emit: coordinates

    script:
    // metadata_csv is the NO_METADATA placeholder (see main.nf) when
    // --metadata wasn't set -- detected by name, not by content, since an
    // input path is always staged for a fixed-arity process input
    def metadata_arg = metadata_csv.name != 'NO_METADATA' ? "--metadata ${metadata_csv}" : ""
    """
    pcoa.py \\
        --abundance ${abundance_table} \\
        ${metadata_arg} \\
        --out pcoa_coordinates.tsv
    """
}
