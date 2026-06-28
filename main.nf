#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * =============================================================================
 *  Metagenome Pipeline v4.0 - Advanced Illumina Shotgun Metagenomics
 * =============================================================================
 *  QC: FastQC | fastp
 *  Assembly: MEGAHIT | QUAST
 *  Binning: MetaBAT2 | CONCOCT | MaxBin2 | DAS Tool
 *  MAG QA: CheckM
 *  MAG Taxonomy: GTDB-Tk
 *  Read Taxonomy: Kraken2/Bracken
 *  Gene Prediction: Prodigal
 *  Functional Annotation: EggNOG-mapper
 *  Functional Profiling: HUMAnN
 *  Diversity Analysis: Alpha (Shannon, Simpson, Chao1, Observed)
 *                      Beta (Bray-Curtis, Jaccard, UniFrac)
 *                      Ordination (PCA, PCoA, NMDS, UMAP, t-SNE)
 *  Differential Abundance: DESeq2 | ANCOM-BC | LEfSe
 *  Visualization: Krona | Sankey | Boxplots | Volcano | Cladograms
 *  Reporting: MultiQC
 * =============================================================================
 */

// ── Input Data ───────────────────────────────────────────────────────────────
params.reads         = "$projectDir/data/samplesheet.csv"
params.metadata      = "$projectDir/data/metadata.csv"
params.outdir        = "./results"

// ── Assembly ─────────────────────────────────────────────────────────────────
params.assembly_mem    = "10"
params.min_contig_len  = 200
params.run_quast       = false

// ── Binning ──────────────────────────────────────────────────────────────────
params.metabat_min_contig_len = 1500
params.run_concoct            = false
params.run_maxbin2            = false
params.run_das_tool           = false

// ── MAG Quality ──────────────────────────────────────────────────────────────
params.checkm_db      = "/path/to/checkm-db"
params.run_checkm     = false

// ── MAG Taxonomy ─────────────────────────────────────────────────────────────
params.gtdbtk_db      = "/path/to/gtdbtk-db"
params.run_gtdbtk     = false

// ── Read Taxonomy ────────────────────────────────────────────────────────────
params.read_length   = 150
params.kraken_db     = "${projectDir}/../k2_database"

// ── Functional ───────────────────────────────────────────────────────────────
params.eggnog_db     = "/path/to/eggnog-db"
params.run_eggnog    = false
params.humann_db     = "/path/to/humann-db"
params.run_humann    = false

// ── Diversity Analysis ───────────────────────────────────────────────────────
params.run_diversity = true

// ── Differential Abundance ───────────────────────────────────────────────────
params.run_diffabund = false

include { FASTQC          } from './modules/fastqc'
include { FASTP           } from './modules/fastp'
include { MEGAHIT         } from './modules/megahit'
include { QUAST           } from './modules/quast'
include { METABAT2        } from './modules/metabat2'
include { CONCOCT         } from './modules/concoct'
include { MAXBIN2         } from './modules/maxbin2'
include { DASTOOL         } from './modules/dastool'
include { CHECKM          } from './modules/checkm'
include { GTDBTK          } from './modules/gtdbtk'
include { KRAKEN2_BRACKEN } from './modules/kraken2'
include { PRODIGAL        } from './modules/prodigal'
include { EGGNOG          } from './modules/eggnog'
include { HUMANN          } from './modules/humann'
include { KRONA           } from './modules/krona'
include { DIVERSITY       } from './modules/diversity'
include { DIFFABUND       } from './modules/diffabund'
include { MULTIQC         } from './modules/multiqc'

workflow {
    ch_samples = Channel.fromPath(params.reads, checkIfExists: true)
        | splitCsv(header: true, sep: ',')
        | map { row ->
            def meta = [id: row.sample_id, single_end: row.single_end.toBoolean()]
            if (meta.single_end) {
                [meta, file(row.read1)]
            } else {
                [meta, [file(row.read1), file(row.read2)]]
            }
        }

    // ── Stage 1: Quality Control ────────────────────────────────────────
    FASTQC(ch_samples)

    // ── Stage 2: Trimming ────────────────────────────────────────────────
    FASTP(ch_samples)
    ch_trimmed = FASTP.out.reads

    // ── Stage 3: Assembly + QC ───────────────────────────────────────────
    MEGAHIT(ch_trimmed)
    ch_assemblies = MEGAHIT.out.assembly

    if (params.run_quast) {
        QUAST(ch_assemblies)
    }

    // ── Stage 4: Binning ─────────────────────────────────────────────────
    METABAT2(ch_assemblies.join(ch_trimmed))
    ch_bins = METABAT2.out.bins

    ch_all_bins = ch_bins

    if (params.run_concoct) {
        CONCOCT(ch_assemblies.join(ch_trimmed))
        ch_all_bins = ch_all_bins.mix(CONCOCT.out.bins)
    }

    if (params.run_maxbin2) {
        MAXBIN2(ch_assemblies.join(ch_trimmed))
        ch_all_bins = ch_all_bins.mix(MAXBIN2.out.bins)
    }

    if (params.run_das_tool) {
        DASTOOL(ch_assemblies.join(ch_all_bins.groupTuple()))
        ch_bins = DASTOOL.out.bins
    }

    // ── Stage 5: MAG quality + taxonomy ──────────────────────────────────
    if (params.run_checkm) {
        CHECKM(ch_bins)
    }

    if (params.run_gtdbtk) {
        GTDBTK(ch_bins)
    }

    // ── Stage 6: Read-based taxonomy ──────────────────────────────────────
    KRAKEN2_BRACKEN(ch_trimmed)

    // ── Stage 7: Gene prediction ─────────────────────────────────────────
    PRODIGAL(ch_assemblies)
    ch_proteins = PRODIGAL.out.proteins

    // ── Stage 8: Functional annotation ────────────────────────────────────
    if (params.run_eggnog) {
        EGGNOG(ch_proteins)
    }

    // ── Stage 9: Functional profiling ────────────────────────────────────
    if (params.run_humann) {
        HUMANN(ch_trimmed)
    }

    // ── Stage 10: Visualization ──────────────────────────────────────────
    KRONA(KRAKEN2_BRACKEN.out.kraken_output)

    // ── Stage 11: Diversity Analysis ─────────────────────────────────────
    if (params.run_diversity) {
        DIVERSITY(
            KRAKEN2_BRACKEN.out.bracken_s.map { meta, report -> report }.collect()
        )
    }

    // ── Stage 12: Differential Abundance ─────────────────────────────────
    if (params.run_diffabund) {
        ch_abundance = KRAKEN2_BRACKEN.out.bracken_s.map { meta, report -> report }.collect()
        DIFFABUND(
            file(params.metadata),
            ch_abundance
        )
    }

    // ── Stage 13: MultiQC report ─────────────────────────────────────────
    MULTIQC(
        FASTQC.out.html.map { meta, html -> html }.collect(),
        FASTQC.out.zip.map { meta, zip -> zip }.collect(),
        FASTP.out.html.map { meta, html -> html }.collect(),
        ch_assemblies.map { meta, assembly -> assembly }.collect()
    )
}

workflow.onComplete {
    log.info """
    ===========================================
       Metagenome Pipeline Complete!
       Output: ${params.outdir}
       Duration: ${workflow.duration}
       Completed: ${workflow.success}
    ===========================================
    """
}
