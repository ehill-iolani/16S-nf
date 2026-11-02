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
    //
    // Each gzipped part is checked first (this image has GNU gzip). A damaged
    // part can't be caught later: the busybox zcat in CHOPPER's container
    // prints "crc error", stops reading and still exits 0, so the sample
    // silently loses every read after the bad part.
    //
    // Plain FASTQ parts (no gzip magic) are gzipped on the way in, so the
    // merged file is always real gzip: that busybox zcat ignores -f and fails
    // on plain text ("no gzip/bzip2/xz magic"), and a plain part cat'ed after
    // gzipped ones would be dropped as trailing garbage. -1 because CHOPPER
    // reads it straight back.
    def infiles = fastq instanceof List ? fastq : [fastq]
        """
        for f in ${infiles.join(' ')}; do
            if [ "\$(head -c 2 "\$f" | od -An -tx1 | tr -d ' \\n')" = "1f8b" ]; then
                if ! msg=\$(gzip -t "\$f" 2>&1); then
                    echo "ERROR: sample ${sample}: input file \$f is damaged and can't be read (gzip -t: \${msg##*: })." >&2
                    echo "Copy it again from the sequencer, or leave it out of the samplesheet." >&2
                    exit 1
                fi
                cat "\$f"
            else
                gzip -1 -c "\$f"
            fi
        done > ${sample}.merged.fastq.gz
        """
}
