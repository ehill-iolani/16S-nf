process CHOPPER {
    tag "$sample"
    label 'process_low'
    container 'quay.io/biocontainers/chopper:0.7.0--hdcf5f25_0'
    publishDir(path: { "${params.outdir}/${sample}/01_filtered" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.filtered.fastq.gz"), emit: filtered

    script:
    """
    set -o pipefail

    zcat -f ${fastq} \\
      | chopper -q ${params.min_qual} -l ${params.min_len} --maxlength ${params.max_len} \\
      | gzip > ${sample}.filtered.fastq.gz
    """
}
