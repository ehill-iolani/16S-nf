include { MAKEBLASTDB      } from '../modules/makeblastdb.nf'
include { MERGE_FASTQ      } from '../modules/merge_fastq.nf'
include { CHOPPER          } from '../modules/chopper.nf'
include { CUTADAPT         } from '../modules/cutadapt.nf'
include { ISONCLUST        } from '../modules/isonclust.nf'
include { SPOA_CONSENSUS   } from '../modules/spoa_consensus.nf'
include { MINIMAP2_ALIGN   } from '../modules/minimap2_align.nf'
include { RACON            } from '../modules/racon.nf'
include { MEDAKA           } from '../modules/medaka.nf'
include { BLAST_TAX        } from '../modules/blast_tax.nf'
include { SORT_CONSENSUS   } from '../modules/sort_consensus.nf'
include { BUILD_REPORT     } from '../modules/report.nf'
include { PCOA             } from '../modules/pcoa.nf'

workflow SIXTEEN_S {

    take:
    reads_ch    // tuple(sample, fastq)
    ref_fasta   // path to 16S reference sequences fasta for BLAST taxonomy
    metadata    // path to a metadata CSV, or the NO_METADATA placeholder

    main:
    // 1. build a BLAST db from the reference fasta (once per run, independent
    //    of per-sample steps below)
    MAKEBLASTDB(ref_fasta)

    // 2. merge multi-part fastq(.gz) files per barcode into one file per sample
    MERGE_FASTQ(reads_ch)

    // 3. length/quality filter
    CHOPPER(MERGE_FASTQ.out.merged)

    // 4. primer trimming (skipped internally if no primers supplied)
    CUTADAPT(CHOPPER.out.filtered)

    // 5. quality-aware de novo clustering -- groups reads from the same
    // organism together before consensus, which is what makes BLAST
    // classification reliable against noisy raw ONT reads
    ISONCLUST(CUTADAPT.out.trimmed)

    // ISONCLUST emits one fastq per cluster per sample; flatten and tag
    clusters_ch = ISONCLUST.out.clusters
        .flatMap { sample, cluster_fastqs ->
            def files = cluster_fastqs instanceof List ? cluster_fastqs : [cluster_fastqs]
            files.collect { fq -> tuple(sample, fq.baseName, fq) }
        }

    // 6. draft consensus per cluster
    SPOA_CONSENSUS(clusters_ch)

    // 7. alignment-based refinement
    MINIMAP2_ALIGN(SPOA_CONSENSUS.out.draft.join(clusters_ch, by: [0, 1]))
    RACON(MINIMAP2_ALIGN.out.aligned)

    // 8. ONT-specific polish -- opt-in via --enable_medaka; off by default, in
    // which case the racon output above is used as the consensus directly
    if (params.enable_medaka) {
        MEDAKA(RACON.out.polished)
        consensus_ch = MEDAKA.out.consensus
    } else {
        consensus_ch = RACON.out.polished
            .map { sample, cluster_id, racon_fasta, cluster_fastq -> tuple(sample, cluster_id, racon_fasta) }
    }

    // 9. taxonomy assignment against the freshly built BLAST db -- the
    // "alignment against a targeted reference" classification approach
    // (wf-metagenomics' minimap2 mode), run here via BLAST since the
    // consensus step above already collapsed each cluster to one polished
    // sequence per organism
    // .first() turns the (single-emission) db channels into value channels so
    // they're reused for every cluster instead of being consumed after one
    BLAST_TAX(consensus_ch, MAKEBLASTDB.out.db_files.first(), MAKEBLASTDB.out.db_name.first())

    // 10. per-run abundance table + QC report -- identical schema to
    // edna-ont-nf's, so the platform's existing abundance/rarefaction charts
    // work against this pipeline's output with no frontend changes
    BUILD_REPORT(
        BLAST_TAX.out.hits.map { sample, cluster_id, hits -> hits }.collect(),
        consensus_ch.map { sample, cluster_id, fasta -> fasta }.collect()
    )

    // 11. beta diversity -- Bray-Curtis + PCoA across samples from that same
    // abundance table, optionally joined with user-supplied metadata for
    // coloring by group/site/etc. in the frontend
    PCOA(BUILD_REPORT.out.report, metadata)

    // 12. gather every consensus fasta into confident / low_confidence / no_hit
    // dirs, using the same best-hit classification as the abundance table
    SORT_CONSENSUS(consensus_ch.join(BLAST_TAX.out.hits, by: [0, 1]))

    emit:
    consensus = consensus_ch
    report    = BUILD_REPORT.out.report
    pcoa      = PCOA.out.coordinates
}
