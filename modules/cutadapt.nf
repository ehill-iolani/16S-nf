// Reverse complement of an IUPAC primer sequence (U read as T)
def reverseComplement(String seq) {
    def comp = [A: 'T', C: 'G', G: 'C', T: 'A', U: 'A', M: 'K', R: 'Y', W: 'W', S: 'S', Y: 'R',
                K: 'M', V: 'B', H: 'D', D: 'H', B: 'V', N: 'N']
    seq.toUpperCase().reverse().collect { comp[it] }.join()
}

process CUTADAPT {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/cutadapt:4.9--py310h1fe012e_3"
    publishDir(path: { "${params.outdir}/${sample}/02_trimmed" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.trimmed.fastq.gz"), emit: trimmed

    script:
    // One linked adapter, "forward primer ... reverse complement of the
    // reverse primer", both required: a read is kept only if it carries the
    // whole amplicon, and is trimmed to what lies between the primers.
    // --revcomp also tries each read's reverse complement and keeps the
    // orientation that matched, so reads sequenced from either strand are
    // trimmed alike (and come out in the forward orientation).
    //
    // --rev_primer is the primer as ordered (5'->3'), so it's reverse-
    // complemented here. Before v0.3.2 it was passed as-is as a 3' adapter
    // (-a), which never matched a forward read's 3' end, and on a reverse-
    // oriented read -- which *starts* with it -- cut away everything after it,
    // leaving a stub that went on to clustering.
    def linked = (params.fwd_primer && params.rev_primer)
        ? "${params.fwd_primer};required...${reverseComplement(params.rev_primer)};required"
        : ""
    """
    if [ -n "${linked}" ]; then
        cutadapt -j ${task.cpus} -g '${linked}' --revcomp --discard-untrimmed \\
            -o ${sample}.trimmed.fastq.gz ${fastq}
    else
        cp ${fastq} ${sample}.trimmed.fastq.gz
    fi
    """
}
