// Read-level QC (length + mean Q) before vs. after filtering. Per-sample
// intermediates stay in the work dir; READ_STATS_REPORT combines them and is
// the only thing published (into final_report/, which downstream UIs already
// list, so no new top-level output dir).
process READ_STATS {
    tag "$sample"
    label 'process_low'
    container 'quay.io/biocontainers/pandas:2.2.1'

    input:
    tuple val(sample), path(raw), path(filtered)

    output:
    tuple val(sample), path("${sample}.read_stats.tsv"), path("${sample}.read_hist.tsv"), emit: stats

    script:
    // raw/filtered can share a basename across samples but never within one
    // task, so the staged names are safe to pass straight through
    """
    read_stats.py summarize \\
        --sample ${sample} \\
        --stage raw=${raw} \\
        --stage filtered=${filtered} \\
        --out-prefix ${sample}
    """
}

process READ_STATS_REPORT {
    label 'process_low'
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir "${params.outdir}/final_report", mode: 'copy'

    input:
    path stats_files
    path hist_files

    output:
    path "read_stats.tsv",        emit: table
    path "read_length_qscore_hist.tsv"
    path "read_qc_summary.html"

    script:
    """
    read_stats.py report \\
        --stats ${stats_files} \\
        --hist ${hist_files} \\
        --out-table read_stats.tsv \\
        --out-hist read_length_qscore_hist.tsv \\
        --out-html read_qc_summary.html
    """
}
