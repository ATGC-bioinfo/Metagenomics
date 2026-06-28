# Metagenomics Pipeline v4.0

Comprehensive Illumina shotgun metagenomics analysis pipeline built with Nextflow.

## Pipeline Overview

```
QC           → FastQC + fastp (trimming)
Assembly     → MEGAHIT + QUAST (quality assessment)
Binning      → MetaBAT2 / CONCOCT / MaxBin2 / DAS Tool
MAG QA       → CheckM
MAG Tax      → GTDB-Tk
Read Tax     → Kraken2 + Bracken (species/genus/phylum)
Gene Pred    → Prodigal
Func Ann     → EggNOG-mapper (optional)
Func Prof    → HUMAnN (optional)
Taxonomy CSV → Report → CSV conversion
Per-sample   → Sankey diagrams + taxonomic summary plots
Krona        → Interactive taxonomic charts
Community    → Cross-sample alpha/beta diversity, PCA, heatmap
Diversity    → Alpha (Shannon, Simpson, Chao1, Observed)
               Beta (Bray-Curtis, Jaccard)
               Ordination (PCA, PCoA, NMDS, UMAP, t-SNE)
               Boxplots / Violin / Heatmap / Dendrogram
Diff Abund   → DESeq2 / ANCOM-BC / LEfSe (optional)
               Volcano / MA / Cladogram / Heatmap
Reporting    → MultiQC
```

## Stages

| Stage | Process | Description |
|-------|---------|-------------|
| 1 | `FASTQC` | Raw read quality reports |
| 2 | `FASTP` | Adapter trimming & quality filtering |
| 3 | `MEGAHIT` | Metagenomic assembly |
| 3b | `QUAST` | Assembly quality statistics (optional) |
| 4 | `METABAT2` | Metagenomic binning |
| 4b | `CONCOCT` / `MaxBin2` / `DAS Tool` | Additional binning + refinement (optional) |
| 5 | `CHECKM` | MAG completeness/contamination (optional) |
| 6 | `GTDB-Tk` | MAG taxonomic classification (optional) |
| 7 | `Kraken2/Bracken` | Read-based taxonomic profiling |
| 8 | `Prodigal` | Gene prediction on assemblies |
| 9 | `EggNOG-mapper` | Functional annotation of proteins (optional) |
| 10 | `HUMAnN` | Functional profiling of reads (optional) |
| 11 | `PLOTS` | Per-sample Sankey diagrams + taxonomic summary PNG |
| 12 | `KRONA` | Interactive taxonomic Krona charts |
| 13 | `COMMUNITY` | Cross-sample alpha/beta bar, PCA, heatmap |
| 14 | `DIVERSITY` | Alpha & beta diversity tables + ordination plots |
| 15 | `DIFFABUND` | Differential abundance (optional) |
| 16 | `MultiQC` | Aggregated HTML report |

## Quick Start

```bash
# 1. Clone and enter the pipeline directory
git clone <repo-url> Metagenomics
cd Metagenomics

# 2. Create conda environment (if not already available)
conda create -n metagenome -c bioconda -c conda-forge \
  fastqc fastp megahit metabat2 kraken2 bracken prodigal multiqc \
  python=3.10 numpy pandas scipy scikit-learn plotly umap-learn
conda activate metagenome
pip install pydeseq2 kaleido

# 3. Point to your Kraken2/Bracken database
#    Edit main.nf or pass via --kraken_db flag

# 4. Prepare samplesheet
cat data/samplesheet.csv
# sample_id,read1,read2,single_end
# SRR39192342,/path/to/SRR39192342_R1.fastq.gz,/path/to/SRR39192342_R2.fastq.gz,false

# 5. Prepare metadata (for differential abundance)
cat data/metadata.csv
# sample_id,group
# SRR39192342,treatment
# SRR39192343,control

# 6. Update conda path in nextflow.config if needed

# 7. Run pipeline
nextflow run main.nf

# Resume after interruption
nextflow run main.nf -resume

# Stub run (test pipeline structure)
nextflow run main.nf -stub-run
```

## Output Structure

```
results/
├── qc/fastqc/               # FastQC HTML reports + zip
├── trimmed/                 # Trimmed FASTQ + fastp reports
├── assembly/{sample}/       # MEGAHIT contigs, QUAST results (optional)
├── bins/{sample}/           # MAG bins, depth profiles
├── taxonomy/{sample}/       # Kraken2 reports, Bracken profiles, CSV copies, Krona HTML
├── annotation/{sample}/     # Prodigal genes/proteins, EggNOG annotations (optional)
├── functional/{sample}/     # HUMAnN gene families & pathways (optional)
├── plots/{sample}/          # Sankey diagrams + taxonomic summary PNG (per sample)
├── community/               # Alpha/beta bar, heatmap, PCA PNG (across samples)
├── diversity/               # Alpha/beta diversity tables + ordination PNG plots
├── diffabund/               # DESeq2/ANCOM-BC/LEfSe results + volcano/MA/cladogram/heatmap PNG
├── multiqc/                 # Aggregated MultiQC report
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | `data/samplesheet.csv` | Samplesheet path |
| `--metadata` | `data/metadata.csv` | Sample metadata (for diff. abundance) |
| `--outdir` | `./results` | Output directory |
| `--kraken_db` | *(see main.nf)* | Kraken2/Bracken database path |
| `--assembly_mem` | `0.9` | MEGAHIT memory (fraction of RAM, 0-1) |
| `--min_contig_len` | `200` | Minimum contig length for MEGAHIT |
| `--read_length` | `150` | Read length for Bracken estimation |
| `--run_quast` | `false` | Enable QUAST assembly QC |
| `--run_eggnog` | `false` | Enable EggNOG-mapper annotation |
| `--run_humann` | `false` | Enable HUMAnN functional profiling |
| `--run_diversity` | `true` | Enable alpha/beta diversity analysis |
| `--run_diffabund` | `true` | Enable differential abundance analysis |

### Optional module flags

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_concoct` | `false` | Enable CONCOCT binning |
| `--run_maxbin2` | `false` | Enable MaxBin2 binning |
| `--run_das_tool` | `false` | Enable DAS Tool bin refinement |
| `--run_checkm` | `false` | Enable CheckM MAG evaluation |
| `--run_gtdbtk` | `false` | Enable GTDB-Tk classification |
| `--checkm_db` | `/path/to/checkm-db` | CheckM database path |
| `--gtdbtk_db` | `/path/to/gtdbtk-db` | GTDB-Tk database path |
| `--eggnog_db` | `/path/to/eggnog-db` | EggNOG-mapper database path |
| `--humann_db` | `/path/to/humann-db` | HUMAnN database path |

## Diversity Outputs

When `--run_diversity` is enabled (default), the following are generated:

```
results/diversity/
├── alpha_diversity.tsv        # Shannon, Simpson, Chao1, Observed per sample
├── beta_diversity.tsv         # Bray-Curtis & Jaccard pairwise matrices
├── alpha_boxplots.png         # Boxplots of all alpha metrics
├── alpha_violin.png           # Violin plots of all alpha metrics
├── beta_bray_curtis.png       # Bray-Curtis dissimilarity heatmap
├── beta_jaccard.png           # Jaccard dissimilarity heatmap
├── beta_dendrogram.png        # Hierarchical clustering dendrogram
├── ordination_pca.png         # Principal Component Analysis
├── ordination_pcoa.png        # Principal Coordinate Analysis (MDS)
├── ordination_nmds.png        # Non-metric Multidimensional Scaling
├── ordination_umap.png        # UMAP embedding (requires umap-learn)
└── ordination_tsne.png        # t-SNE embedding
```

## Differential Abundance Outputs

When `--run_diffabund` is enabled with a valid `--metadata` file:

```
results/diffabund/
├── deseq2_results.tsv         # DESeq2 results table
├── ancombc_results.tsv        # ANCOM-BC results table
├── lefse_results.tsv          # LEfSe results table
├── deseq2_volcano.png         # Volcano plot (log2FC vs -log10 p-value)
├── deseq2_ma.png              # MA plot (mean abundance vs log2FC)
├── lefse_cladogram.png        # Biomarker cladogram
└── sig_heatmap.png            # Significant taxa heatmap
```

## Dependencies

Core tools available in the `metagenome` conda environment:

- **QC**: fastp, FastQC, MultiQC
- **Assembly**: MEGAHIT, QUAST
- **Binning**: MetaBAT2, CONCOCT, MaxBin2, DAS Tool
- **MAG**: CheckM, GTDB-Tk
- **Taxonomy**: Kraken2, Bracken
- **Annotation**: Prodigal, EggNOG-mapper, HUMAnN
- **Diversity**: Python (numpy, pandas, scipy, scikit-learn, plotly, umap-learn)
- **Diff. Abundance**: Python (PyDESeq2, scipy, statsmodels, plotly)

## Configuration

Edit `nextflow.config`:

- **`process.conda`** — path to your conda environment (update to match your system)
- **Profiles** — `standard` (local), `slurm` (SLURM), or `cluster` (PBS Pro)
- **Resource labels** — CPU/memory/time per process type

## License

MIT
