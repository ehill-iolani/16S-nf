process CUTADAPT {
    tag "$sample"
    label 'process_low'
    container 'quay.io/biocontainers/cutadapt:4.9--py310h1fe012e_3'
    publishDir(path: { "${params.outdir}/${sample}/02_trimmed" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.trimmed.fastq.gz"), emit: trimmed

    script:
    // --times 2: without it, cutadapt removes only the single best-matching
    // adapter per read, so the 5' primer (-g) or the 3' primer (-a) -- not
    // both -- would silently survive on most reads
    def primer_args = (params.fwd_primer && params.rev_primer)
        ? "-g ${params.fwd_primer} -a ${params.rev_primer} --times 2 --discard-untrimmed"
        : ""
    """
    if [ -n "${primer_args}" ]; then
        cutadapt ${primer_args} -o ${sample}.trimmed.fastq.gz ${fastq}
    else
        cp ${fastq} ${sample}.trimmed.fastq.gz
    fi
    """
}
