process MERGE_FASTQ {
    tag "$sample"
    label 'process_low'
    publishDir(path: { "${params.outdir}/${sample}/00_merged" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.merged.fastq.gz"), emit: merged

    script:
    // ONT barcode folders are usually split into many part-files
    // (fastq_pass/barcodeNN/*.fastq.gz); cat them into one file per sample
    // so downstream steps don't need to care how many parts came in.
    def infiles = fastq instanceof List ? fastq : [fastq]
        """
        cat ${infiles.join(' ')} > ${sample}.merged.fastq.gz
        """
}
