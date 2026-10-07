process VSEARCH_CLUSTER {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/vsearch:2.30.6--h0bb26bb_0"
    publishDir(path: { "${params.outdir}/${sample}/03_clusters" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("clusters/*.fastq"), emit: clusters
    tuple val(sample), path("clusters.uc"), emit: cluster_uc

    script:
    /*
     * NOTE: cluster_id is vsearch's --id (fraction identity to the cluster
     * centroid) and is the most sensitivity-critical parameter in the whole
     * pipeline -- too loose merges species, too strict splits one species into
     * many small clusters that then fall under --min_cluster and are dropped.
     * ONT reads carry real sequencing error, so it has to sit comfortably
     * below (1 - per-read error rate); re-tune per primer set / basecaller.
     *
     * Greedy identity clustering (like decona's CD-HIT step) rather than
     * isONclust: isONclust was built for transcriptome gene families and does
     * not scale to millions of short amplicon reads (a 2.6M-read barcode ran
     * past the 8h task limit); vsearch clusters the same in minutes.
     */
    """
    set -o pipefail

    # Strip everything but the bare read id from each header. Modern basecaller
    # (Dorado) headers carry tab-separated tags after the read id (qs:f:.. mx:i:..
    # etc), and vsearch only truncates labels at the first space, not a tab.
    # (input is always gzip -- CUTADAPT emits .fastq.gz -- and this image's
    # busybox zcat has no -f plain-text passthrough)
    zcat ${fastq} | awk 'NR % 4 == 1 { sub(/[ \\t].*/, "", \$0) } { print }' > input.fastq
    if [ ! -s input.fastq ]; then
        echo "VSEARCH_CLUSTER: no reads in ${fastq}" >&2
        exit 1
    fi

    # fastq -> fasta for clustering (vsearch's fastq parser would also reject
    # modern Q-scores above its default --fastq_qmax of 41, so avoid it)
    awk 'NR % 4 == 1 { print ">" substr(\$0, 2) } NR % 4 == 2 { print }' input.fastq > reads.fasta

    # --strand both: amplicon reads come off the pore in either orientation
    vsearch --cluster_fast reads.fasta \\
        --id ${params.cluster_id} \\
        --strand both \\
        --threads ${task.cpus} \\
        --uc clusters.uc

    # Split the original fastq (qualities intact, needed by racon/medaka) by
    # cluster. clusters.uc records: S = centroid (new cluster), H = member hit;
    # field 2 = cluster number, field 5 = strand vs. the centroid ('+'/'-'),
    # field 9 = read id. Clusters under --min_cluster reads are dropped here.
    mkdir -p clusters
    awk -F'\\t' -v min=${params.min_cluster} '
        NR == FNR {
            if (\$1 == "S" || \$1 == "H") { cl[\$9] = \$2; st[\$9] = \$5; n[\$2]++ }
            next
        }
        FNR % 4 == 1 { id = substr(\$0, 2) }
        FNR % 4 == 2 { seq = \$0 }
        FNR % 4 == 0 {
            if (id in cl) {
                c = cl[id]
                if (n[c] >= min) {
                    out = (st[id] == "-") ? "minus.tsv" : "plus.tsv"
                    print c "\\t" id "\\t" seq "\\t" \$0 > out
                }
            }
        }
    ' clusters.uc input.fastq
    touch plus.tsv minus.tsv

    # members that matched the centroid on the reverse strand get reverse-
    # complemented (sequence) / reversed (qualities) so every read in a cluster
    # shares the centroid's orientation -- spoa needs that to build one draft
    cut -f1,2 minus.tsv > minus.meta
    cut -f3 minus.tsv | rev | tr 'ACGTacgt' 'TGCAtgca' > minus.seq
    cut -f4 minus.tsv | rev > minus.qual
    paste minus.meta minus.seq minus.qual > minus.oriented.tsv

    # append + close per record so a run with thousands of clusters never
    # holds more than one output file open at a time
    cat plus.tsv minus.oriented.tsv | awk -F'\\t' '{
        out = "clusters/" \$1 ".fastq"
        printf "@%s\\n%s\\n+\\n%s\\n", \$2, \$3, \$4 >> out
        close(out)
    }'

    if ! ls clusters/*.fastq > /dev/null 2>&1; then
        echo "VSEARCH_CLUSTER: no cluster reached --min_cluster ${params.min_cluster} reads at --cluster_id ${params.cluster_id}; see clusters.uc (lower either one?)" >&2
        exit 1
    fi
    echo "kept \$(ls clusters/*.fastq | wc -l) clusters (>= ${params.min_cluster} reads), \$(( \$(cat clusters/*.fastq | wc -l) / 4 )) of \$(( \$(wc -l < input.fastq) / 4 )) reads" >&2
    """
}
