#!/usr/bin/env nextflow
/*
 * 16S-nf: ONT 16S rRNA amplicon classification for microbial community
 * profiling. Same clustering + consensus + BLAST-taxonomy architecture as
 * edna-ont-nf (isONclust + spoa/racon/medaka + BLAST), loosely following
 * wf-metagenomics' stages (ingest -> filter -> classify -> abundance +
 * diversity report) but classifying via alignment against a targeted
 * reference (wf-metagenomics' "minimap2 mode", here via BLAST) rather than
 * Kraken2 -- a better fit once noisy per-read ONT data has been
 * consensus-polished, and it keeps the reference a plain FASTA a
 * non-specialist can swap out, same as edna-ont-nf's --taxdb.
 *
 * Emits the exact same abundance_table.tsv schema as edna-ont-nf, plus an
 * optional pcoa_coordinates.tsv for beta-diversity -- both consumed as-is by
 * the iolani-bioinformatics frontend's existing chart components.
 */

nextflow.enable.dsl = 2

include { SIXTEEN_S   } from './workflows/sixteen_s.nf'
include { MERGE_FASTQ } from './modules/merge_fastq.nf'

// ---- top-level params (override via -params-file or --flag) ----
// filtering/clustering/taxonomy defaults live in nextflow.config, not here --
// anything read directly inside an included module or workflow script
// (rather than only inside this file's own `workflow` block) must be set
// there to be reliably visible by the time that module/workflow script runs
params.input        = null   // path to samplesheet.csv (sample,fastq); fastq globs resolve relative to the launch dir, not the CSV's location
params.outdir       = "results"
params.fwd_primer   = null
params.rev_primer   = null
params.taxdb        = null   // path to a 16S reference sequences FASTA (e.g. a SILVA/NCBI 16S export); a BLAST db is built from this at runtime
params.metadata     = null   // optional CSV: sample,<arbitrary metadata columns> -- joined onto the PCoA output for coloring
params.help         = false

def helpMessage() {
    log.info """
    16S-nf: ONT 16S rRNA amplicon microbial community classification
    ------------------------------------------------------------------
    Usage:
      nextflow run main.nf --input samplesheet.csv --taxdb silva_16s.fasta --outdir results

    Required:
      --input       CSV: sample,fastq_path
      --taxdb       16S reference sequences FASTA (BLAST db is built from this each run)

    Key optional:
      --fwd_primer / --rev_primer   primer sequences for cutadapt trimming (e.g. 27F/1492R for full-length 16S)
      --min_len / --max_len / --min_qual   chopper filtering thresholds (defaults tuned for full-length ~1.5kb 16S)
      --cluster_id   isONclust similarity threshold (default ${params.cluster_id})
      --min_cluster  minimum reads to polish a cluster (default ${params.min_cluster})
      --enable_medaka   use medaka-polished consensus instead of racon consensus (default ${params.enable_medaka})
      --min_pident   BLAST %identity below which a hit is low-confidence (default ${params.min_pident})
      --metadata     CSV: sample,<metadata columns> -- joined onto the PCoA coordinates for coloring by group/site/etc.
    """.stripIndent()
}

def samplesheetToChannel(path) {
    Channel
        .fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // checkIfExists doesn't catch a glob that matches zero files, so check explicitly --
            // otherwise this silently produces an empty file list downstream instead of failing
            def fq = file(row.fastq, checkIfExists: true)
            def files = fq instanceof List ? fq : [fq]
            if (files.isEmpty()) {
                error "Samplesheet row '${row.sample}': no files matched '${row.fastq}' (glob paths resolve relative to the launch directory, not the samplesheet's location)"
            }
            tuple(row.sample, fq)
        }
}

workflow {
    if (params.help || !params.input || !params.taxdb) {
        helpMessage()
        exit 0
    }

    taxdb_fasta_ch = Channel.fromPath(params.taxdb, checkIfExists: true)

    // PCOA's metadata input is optional -- when unset, pass a placeholder
    // file instead of an empty channel so the process (which has a fixed
    // two-input signature) always receives something to stage; the script
    // checks the staged filename to tell the two cases apart.
    metadata_ch = params.metadata
        ? Channel.fromPath(params.metadata, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_METADATA")

    SIXTEEN_S(samplesheetToChannel(params.input), taxdb_fasta_ch, metadata_ch)
}

// dev/debug entry point -- run a single step in isolation, e.g.:
//   nextflow run main.nf -entry MERGE_ONLY --input samplesheet.csv -profile docker
workflow MERGE_ONLY {
    if (!params.input) {
        log.error "MERGE_ONLY requires --input samplesheet.csv"
        exit 1
    }
    MERGE_FASTQ(samplesheetToChannel(params.input))
}
