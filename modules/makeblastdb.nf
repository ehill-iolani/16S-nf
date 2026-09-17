process MAKEBLASTDB {
    label 'process_low'
    container 'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'
    publishDir "${params.outdir}/blastdb", mode: 'copy'

    input:
    path ref_fasta

    output:
    path "blastdb.*", emit: db_files
    val "blastdb", emit: db_name

    script:
    """
    makeblastdb -in ${ref_fasta} -dbtype nucl -parse_seqids -out blastdb
    """
}
