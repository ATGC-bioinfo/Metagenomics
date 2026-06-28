# Metagenomics Pipeline v4.0

Comprehensive Illumina shotgun metagenomics analysis pipeline built with Nextflow.

## Pipeline Overview

```
QC         → FastQC + fastp (trimming)
Assembly   → MEGAHIT + QUAST (quality assessment)
Binning    → MetaBAT2 / CONCOCT / MaxBin2 / DAS Tool
MAG QA     → CheckM
MAG Tax    → GTDB-Tk
Read Tax   → Kraken2 + Bracken (species/genus/phylum)
Gene Pred  → Prodigal
Func Ann   → EggNOG-mapper
Func Prof  → HUMAnN
Diversity  → Alpha (Shannon, Simpson, Chao1, Observed)
             Beta (Bray-Curtis, Jaccard)
             Ordination (PCA, PCoA, NMDS, UMAP, t-SNE)
Diff Abund → DESeq2 / ANCOM-BC / LEfSe
             Volcano / MA / Cladogram / Heatmap
Reporting  → MultiQC
Visual     → Krona / Sankey / Boxplots / Violin / ...
```

## Stages

| Stage | Process | Description |
|-------|---------|-------------|
| 1 | `FASTQC` | Raw read quality reports |
| 2 | `FASTP` | Adapter trimming & quality filtering |
| 3 | `MEGAHIT` | Metagenomic assembly |
| 3b | `QUAST` | Assembly quality statistics (optional) |
| 4 | `METABAT2` | Metagenomic binning |
| 5 | `CHECKM` | MAG completeness/contamination (optional) |
| 6 | `GTDB-Tk` | MAG taxonomic classification (optional) |
| 7 | `Kraken2/Bracken` | Read-based taxonomic profiling |
| 8 | `Prodigal` | Gene prediction on assemblies |
| 9 | `EggNOG-mapper` | Functional annotation of proteins (optional) |
| 10 | `HUMAnN` | Functional profiling of reads (optional) |
| 11 | `Krona` | Interactive taxonomic Krona charts |
| 12 | `DIVERSITY` | Alpha & beta diversity + ordination plots |
| 13 | `DIFFABUND` | Differential abundance (DESeq2, ANCOM-BC, LEfSe) |
| 14 | `MultiQC` | Aggregated HTML report |

## Quick Start

```bash
# Install dependencies (conda)
conda env create -f environment.yml
conda activate metagenome

# Prepare samplesheet
cat data/samplesheet.csv
# sample_id,read1,read2,single_end
# SRR39192342,/path/to/SRR39192342_R1.fastq.gz,/path/to/SRR39192342_R2.fastq.gz,false

# Prepare metadata (for differential abundance)
cat data/metadata.csv
# sample_id,group
# SRR39192342,treatment
# SRR39192343,control

# Run pipeline
nextflow run main.nf

# Resume after interruption
nextflow run main.nf -resume

# Stub run (test pipeline structure)
nextflow run main.nf -stub-run
```

## Output Structure

```
results/
├── qc/fastqc/           # FastQC HTML reports
├── trimmed/             # Trimmed FASTQ files
├── assembly/{sample}/   # MEGAHIT contigs, QUAST results
├── bins/{sample}/       # MAG bins from MetaBAT2 et al.
├── taxonomy/{sample}/   # Kraken2 reports, Bracken profiles, Krona HTML
├── annotation/{sample}/ # Prodigal genes/proteins, EggNOG annotations
├── functional/{sample}/ # HUMAnN gene families & pathways
├── diversity/           # Alpha/beta diversity tables + HTML plots
├── diffabund/           # DESeq2/ANCOM-BC/LEfSe results + visualizations
├── community/           # Community comparison plots
├── plots/{sample}/      # Per-sample Sankey & summary plots
└── multiqc/             # Aggregated MultiQC report
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | `data/samplesheet.csv` | Samplesheet path |
| `--metadata` | `data/metadata.csv` | Sample metadata (for diff. abundance) |
| `--outdir` | `./results` | Output directory |
| `--kraken_db` | `../k2_database` | Kraken2/Bracken database path |
| `--run_quast` | `false` | Enable QUAST assembly QC |
| `--run_eggnog` | `false` | Enable EggNOG-mapper annotation |
| `--run_humann` | `false` | Enable HUMAnN functional profiling |
| `--run_diversity` | `true` | Enable diversity analysis |
| `--run_diffabund` | `false` | Enable differential abundance |
| `--run_checkm` | `false` | Enable CheckM MAG evaluation |
| `--run_gtdbtk` | `false` | Enable GTDB-Tk classification |

## Dependencies

Core tools installed in the `metagenome` conda environment:

- **QC**: fastp, FastQC, MultiQC
- **Assembly**: MEGAHIT, QUAST
- **Binning**: MetaBAT2, CONCOCT, MaxBin2, DAS Tool
- **MAG**: CheckM, GTDB-Tk
- **Taxonomy**: Kraken2, Bracken
- **Annotation**: Prodigal, EggNOG-mapper, HUMAnN
- **Diversity**: Python (numpy, pandas, scipy, scikit-learn, plotly, umap-learn)
- **Diff. Abundance**: Python (PyDESeq2, scipy, statsmodels, plotly)

## Configuration

Edit `nextflow.config` to set:

- **`process.conda`** — path to your conda environment
- **Profiles** — `standard` (local), `slurm` (SLURM cluster), or `cluster` (PBS Pro)
- **Resource labels** — CPU/memory/time per process type

## License

MIT
