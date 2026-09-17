process BUILD_REPORT {
    label 'process_low'
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir "${params.outdir}/final_report", mode: 'copy'

    input:
    path hit_files
    path consensus_files

    output:
    path "abundance_table.tsv", emit: report
    path "run_qc_summary.html"

    script:
    """
    build_report.py \\
        --hits ${hit_files} \\
        --consensus ${consensus_files} \\
        --out-table abundance_table.tsv \\
        --out-html run_qc_summary.html \\
        --min-pident ${params.min_pident}
    """
}
