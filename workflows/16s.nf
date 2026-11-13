include { MAKEBLASTDB      } from '../modules/makeblastdb.nf'
include { MERGE_FASTQ      } from '../modules/merge_fastq.nf'
include { CHOPPER          } from '../modules/chopper.nf'
include { READ_STATS; READ_STATS_REPORT } from '../modules/read_stats.nf'
include { CUTADAPT         } from '../modules/cutadapt.nf'
include { VSEARCH_CLUSTER  } from '../modules/vsearch_cluster.nf'
include { SPOA_CONSENSUS   } from '../modules/spoa_consensus.nf'
include { MINIMAP2_ALIGN   } from '../modules/minimap2_align.nf'
include { RACON            } from '../modules/racon.nf'
include { MEDAKA           } from '../modules/medaka.nf'
include { MERGE_CONSENSUS  } from '../modules/merge_consensus.nf'
include { BLAST_TAX        } from '../modules/blast_tax.nf'
include { BUILD_REPORT     } from '../modules/report.nf'
include { PCOA             } from '../modules/pcoa.nf'

workflow WF_16S {

    take:
    reads_ch    // tuple(sample, fastq)
    ref_fasta   // path to 16S reference sequences fasta for BLAST taxonomy
    taxmap      // SILVA's taxonomy rank map, or the NO_TAXMAP placeholder
    metadata    // path to a metadata CSV, or the NO_METADATA placeholder

    main:
    // 1. build a BLAST db from the reference fasta (once per run, independent
    //    of per-sample steps below)
    MAKEBLASTDB(ref_fasta)

    // 2. merge multi-part fastq(.gz) files per barcode into one file per sample
    MERGE_FASTQ(reads_ch)

    // 3. length/quality filter
    CHOPPER(MERGE_FASTQ.out.merged)

    // 3b. read length / Q-score summary, before vs. after filtering -- a
    // side branch, nothing downstream depends on it. Publishes
    // final_report/read_stats.tsv, read_length_qscore_hist.tsv and
    // read_qc_summary.html, which the platform's Read QC tab charts.
    if (params.enable_read_stats) {
        READ_STATS(MERGE_FASTQ.out.merged.join(CHOPPER.out.filtered))
        READ_STATS_REPORT(
            READ_STATS.out.stats.map { sample, stats, hist -> stats }.collect(),
            READ_STATS.out.stats.map { sample, stats, hist -> hist }.collect()
        )
    }

    // 4. primer trimming (skipped internally if no primers supplied)
    CUTADAPT(CHOPPER.out.filtered)

    // 5. de novo identity clustering (vsearch) -- groups reads from the same
    // organism together before consensus, which is what makes BLAST
    // classification reliable against noisy raw ONT reads
    VSEARCH_CLUSTER(CUTADAPT.out.trimmed)

    // A sample with nothing to cluster (an empty blank, a failed barcode) is
    // skipped there rather than failing the run: it emits no clusters, so it
    // drops out of every step below, and BUILD_REPORT lists it in
    // final_report/skipped_samples.tsv. Only a run where every sample came out
    // empty stops here, since there is nothing left to report on.
    VSEARCH_CLUSTER.out.skipped.subscribe { f ->
        def (name, reads, largest, reason) = f.text.trim().split('\t')
        log.warn(reason == 'no_reads'
            ? "Sample ${name}: no reads left after primer trimming; skipped"
            : "Sample ${name}: no cluster reached --min_cluster ${params.min_cluster} reads (${reads} reads, largest cluster ${largest}); skipped")
    }
    // (the isAborted check: a run that already failed elsewhere also ends
    // with this channel empty, and shouldn't blame the samples for it)
    clusters_ch = VSEARCH_CLUSTER.out.clusters.ifEmpty {
        if (!workflow.session.isAborted()) error "No sample produced a cluster of at least --min_cluster ${params.min_cluster} reads, so there is nothing to identify. Each sample's 01_filtered and 02_trimmed folders show how many reads were left: check the primers, --min_len / --max_len and --cluster_id"
    }

    // Everything from here to BLAST_TAX runs once per sample, not once per
    // cluster: VSEARCH_CLUSTER emits all of a sample's cluster fastqs as one
    // list, and each step below loops over them inside its task. A task per
    // cluster meant thousands of Batch VMs (start-up, image pull, GCS staging)
    // each doing seconds of work. Channels are tuple(sample, [files]) throughout.

    // 6. draft consensus per cluster
    SPOA_CONSENSUS(clusters_ch)

    // 7. alignment-based refinement
    MINIMAP2_ALIGN(SPOA_CONSENSUS.out.draft.join(clusters_ch))
    RACON(MINIMAP2_ALIGN.out.aligned)

    // A cluster none of whose reads aligned to its draft is skipped there
    // rather than failing the run (racon has nothing to polish with), and
    // BUILD_REPORT lists it in final_report/skipped_clusters.tsv. A sample
    // whose every cluster went that way drops out like an empty one; only a
    // run with no polished cluster at all stops here.
    RACON.out.skipped.subscribe { f ->
        f.readLines().findAll { it.trim() }.each { line ->
            def (name, cid, reads, draft_len, reason) = line.split('\t')
            log.warn("Sample ${name}: cluster ${cid} (${reads} reads, ${draft_len} bp draft) skipped, none of its reads aligned to its draft")
        }
    }
    polished_ch = RACON.out.polished.ifEmpty {
        if (!workflow.session.isAborted()) error "No cluster in any sample could be polished: none of their reads aligned to their own draft consensus. Each sample's 02_trimmed and 04_draft folders show the read and draft lengths: check the primers and --min_len / --max_len"
    }

    // 8. ONT-specific polish -- opt-in via --enable_medaka; off by default, in
    // which case the racon output above is used as the consensus directly
    if (params.enable_medaka) {
        MEDAKA(polished_ch)
        consensus_ch = MEDAKA.out.consensus
    } else {
        consensus_ch = polished_ch
            .map { sample, racon_fastas, cluster_fastqs -> tuple(sample, racon_fastas) }
    }

    // 8b. second vsearch pass: fold clusters whose consensus sequences are
    // near-identical back together (one organism split across several read
    // clusters by raw ONT error) -- skipped with --merge_id 0
    if (params.merge_id) {
        MERGE_CONSENSUS(consensus_ch)
        // one fasta per merged group, named for its representative cluster
        consensus_ch = MERGE_CONSENSUS.out.merged
    }

    // 9. taxonomy assignment against the freshly built BLAST db -- the
    // "alignment against a targeted reference" classification approach
    // (wf-metagenomics' minimap2 mode), run here via BLAST since the
    // consensus step above already collapsed each cluster to one polished
    // sequence per organism. Each sample's consensus sequences (one per fasta)
    // are split into chunks of --blast_chunk_size, one blastn task per chunk:
    // at several seconds per query against SILVA, a soil sample's ~3k
    // sequences took over 4 h as one task, and spread over VMs they run side
    // by side. The hits are the same either way. A sample with fewer
    // sequences than one chunk still gets a single task.
    // .first() turns the (single-emission) db channels into value channels so
    // they're reused for every chunk instead of being consumed after one
    def chunk_size = params.blast_chunk_size as int
    blast_in = consensus_ch.flatMap { sample, fastas ->
        def files = [fastas].flatten()
        def chunks = chunk_size > 0 ? files.collate(chunk_size) : [files]
        chunks.withIndex().collect { chunk, i -> tuple(sample, String.format('%03d', i + 1), chunk) }
    }
    BLAST_TAX(blast_in, MAKEBLASTDB.out.db_files.first(), MAKEBLASTDB.out.db_name.first())

    // 10. per-run abundance table + QC report -- identical schema to
    // edna-ont-nf's, so the platform's existing abundance/rarefaction charts
    // work against this pipeline's output with no frontend changes. Also
    // gathers every consensus fasta into confident / low_confidence / no_hit
    // dirs, using the same best-hit call as the abundance table.
    BUILD_REPORT(
        BLAST_TAX.out.hits.map { sample, hits -> hits }.collect(),
        consensus_ch
            .flatMap { sample, fastas -> [fastas].flatten() }
            .collect(),
        taxmap,
        // [] when nothing was skipped, so BUILD_REPORT still runs
        VSEARCH_CLUSTER.out.skipped.collect().ifEmpty([]),
        RACON.out.skipped.collect().ifEmpty([])
    )

    // 11. beta diversity -- Bray-Curtis + PCoA across samples from that same
    // abundance table, optionally joined with user-supplied metadata for
    // coloring by group/site/etc. in the frontend
    PCOA(BUILD_REPORT.out.report, metadata)

    emit:
    consensus = consensus_ch
    report    = BUILD_REPORT.out.report
    pcoa      = PCOA.out.coordinates
}
