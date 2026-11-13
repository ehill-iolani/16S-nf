process BUILD_REPORT {
    label 'process_low'
    container "${params.container_registry}/biocontainers/pandas:2.2.1"
    publishDir "${params.outdir}/final_report", mode: 'copy', pattern: '{abundance_table.tsv,run_qc_summary.html,skipped_samples.tsv,skipped_clusters.tsv}'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'consensus_by_confidence'
    publishDir "${params.outdir}", mode: 'copy', pattern: '*/07_taxonomy/*.hits.tsv'

    input:
    path hit_files
    path consensus_files
    path silva_taxmap   // SILVA's rank map, or the NO_TAXMAP placeholder
    path skipped    // VSEARCH_CLUSTER's <sample>.skipped.tsv files, or [] if none
    path skipped_clusters  // RACON's <sample>.racon_skipped.tsv files, or [] if none

    output:
    path "abundance_table.tsv", emit: report
    path "run_qc_summary.html"
    path "skipped_samples.tsv", optional: true
    path "skipped_clusters.tsv", optional: true
    path "consensus_by_confidence"
    path "*/07_taxonomy/*.hits.tsv"

    script:
    // detected by name, as PCOA does with NO_METADATA
    def taxmap_arg = silva_taxmap.name != 'NO_TAXMAP' ? "--silva-taxmap ${silva_taxmap}" : ""
    def skipped_arg = skipped ? "--skipped ${skipped}" : ''
    def skipped_clusters_arg = skipped_clusters ? "--skipped-clusters ${skipped_clusters}" : ''

    // BLAST_TAX runs per chunk of a sample's sequences; each sample's chunk
    // files are joined back into one <sample>/07_taxonomy/<sample>.hits.tsv
    // here, the one step that sees them all. Chunk numbers are zero-padded,
    // so the glob keeps them in order.
    """
    for f in *.chunk*.hits.tsv; do
        s=\${f%.chunk*.hits.tsv}
        mkdir -p "\$s/07_taxonomy"
        cat "\$f" >> "\$s/07_taxonomy/\$s.hits.tsv"
    done

    build_report.py \\
        --hits ${hit_files} \\
        --consensus ${consensus_files} \\
        --out-table abundance_table.tsv \\
        --out-html run_qc_summary.html \\
        --out-sorted consensus_by_confidence \\
        --min-pident ${params.min_pident} \\
        --min-abundance ${params.min_abundance} \\
        --min-rel-abundance ${params.min_rel_abundance} \\
        ${taxmap_arg} \\
        ${skipped_arg} \\
        ${skipped_clusters_arg}
    """
}
