# Pipeline Results: 16S Gut Metagenome — Duodenal Aspirates

## Sample Context

| Field | Value |
|-------|-------|
| **Sample** | SRX33935698 — 161_Duodenal_Female |
| **Subject** | 161_Duodenal_aspirates (SAMN60940363) |
| **Study** | Human gut metagenome: Males vs Females (PRJNA1279148) |
| **Source** | Duodenal aspirate, female subject |
| **Sequencing** | Illumina MiSeq, V3-V4 16S rRNA (515f/806r), 2×250 bp paired-end |
| **Design** | 16S amplicon sequencing of V3-V4 region using barcoded primer set |
| **DNA Extraction** | MagAttract Power Soil DNA KF kit |
| **Submitted by** | Cedars-Sinai Medical Center |
| **Spots / Bases** | 260,101 spots / 136.6 M bases / 70.8 Mb download |

Two sequencing runs were processed through the pipeline as separate samples:

| Sample ID | SRA Run |
|-----------|---------|
| SRR39192342 | Subject sample (Run 1) |
| SRR39192343 | Subject sample (Run 2) |

---

## Pipeline Output

### 1. Quality Control — `results/qc/fastqc/`

| File | Description |
|------|-------------|
| `*_fastqc.html` | Per-sample FastQC quality reports (raw reads) |

### 2. Trimming — `results/trimmed/`

| File | Description |
|------|-------------|
| `*_R1.trimmed.fastq.gz` | Trimmed forward reads |
| `*_R2.trimmed.fastq.gz` | Trimmed reverse reads |
| `*_fastp.html` | fastp trimming report (HTML) |

### 3. Assembly — `results/assembly/<sample>/`

| File | Description |
|------|-------------|
| `*.contigs.fa` | MEGAHIT metagenomic assembly contigs |

Metagenomic assembly was performed with **MEGAHIT** using paired-end reads. Contigs represent the combined genomic content of the microbial community and serve as the substrate for downstream binning and gene prediction.

### 4. Binning — `results/bins/<sample>/bins/`

MetaBAT2 binning was attempted. No high-quality bins were recovered (`bin.skipped.fa` placeholder), which is expected for low-complexity 16S amplicon data where genome-level binning is not feasible.

### 5. Taxonomy — `results/taxonomy/<sample>/`

Kraken2 + Bracken taxonomic classification at species (S), genus (G), and phylum (P) levels.

| File | Description |
|------|-------------|
| `*.kraken2.output` | Raw Kraken2 classification output |
| `*.kraken2.csv` | Kraken2 report parsed as CSV |
| `*.bracken.S.csv` | Bracken species-level abundance (CSV) |
| `*.bracken.G.csv` | Bracken genus-level abundance (CSV) |
| `*.bracken.P.csv` | Bracken phylum-level abundance (CSV) |
| `*.krona.html` | Interactive Krona taxonomic visualization |

**Top species detected (SRR39192342):**

| Species | Est. Reads | Fraction |
|---------|-----------|----------|
| *Streptococcus thermophilus* | 98,424 | 41.0% |
| *Streptococcus pneumoniae* | 4,135 | 1.7% |
| *Granulicatella adiacens* | 6,604 | 2.8% |
| *Streptococcus infantis* | 359 | 0.15% |
| *Streptococcus acidominimus* | 272 | 0.11% |

The community is dominated by **Firmicutes** (primarily *Streptococcus* spp.), consistent with an upper gastrointestinal / duodenal origin.

### 6. Gene Prediction — `results/annotation/<sample>/`

| File | Description |
|------|-------------|
| `*.fna` | Predicted gene nucleotide sequences (Prodigal) |
| `*.faa` | Predicted protein amino acid sequences (Prodigal) |
| `*.gff` | Gene annotation GFF file |

### 7. Per-Sample Plots — `results/plots/<sample>/`

| File | Description |
|------|-------------|
| `*.sankey.png` | Sankey diagram: Phylum → Genus flow |
| `*.summary.png` | Taxonomic summary: Top phyla, genera, species bar charts |

### 8. Community Analysis — `results/community/`

| File | Description |
|------|-------------|
| `community_alpha.png` | Alpha diversity bar charts (Shannon, Simpson, Chao1) |
| `community_beta.png` | Bray-Curtis dissimilarity heatmap |
| `community_pca.png` | Principal Component Analysis ordination |
| `community_heatmap.png` | Top 20 species relative abundance heatmap |

**Alpha Diversity:**

| Sample | Shannon | Simpson | Chao1 | Observed Species |
|--------|---------|---------|-------|-----------------|
| SRR39192342 | 3.41 | 0.80 | 121 | 121 |
| SRR39192343 | 3.71 | 0.86 | 164 | 164 |

The **Shannon index** (3.4–3.7) indicates moderate microbial diversity. The **Simpson index** (0.80–0.86) suggests a community with some dominant taxa, consistent with the *Streptococcus* predominance seen in taxonomic profiling. A higher richness was observed in the second sequencing run (164 vs 121 observed species).

### 9. Diversity Analysis — `results/diversity/`

| File | Description |
|------|-------------|
| `alpha_diversity.tsv` | Shannon, Simpson, Chao1, Observed per sample |
| `beta_diversity.tsv` | Bray-Curtis & Jaccard pairwise dissimilarity matrices |
| `alpha_boxplots.png` | Boxplots of alpha diversity metrics |
| `alpha_violin.png` | Violin plots of alpha diversity metrics |
| `beta_bray-curtis.png` | Bray-Curtis heatmap |
| `beta_jaccard.png` | Jaccard dissimilarity heatmap |
| `beta_dendrogram.png` | Hierarchical clustering dendrogram (UPGMA) |
| `ordination_pca.png` | PCA ordination |
| `ordination_pcoa.png` | PCoA (MDS on Bray-Curtis) |
| `ordination_nmds.png` | NMDS ordination |
| `ordination_tsne.png` | t-SNE embedding |
| `ordination_umap.png` | UMAP embedding |

### 10. Differential Abundance — `results/diffabund/`

| File | Description |
|------|-------------|
| `deseq2_results.tsv` | DESeq2-style differential abundance results |
| `ancombc_results.tsv` | ANCOM-BC results (CLR-based) |
| `lefse_results.tsv` | LEfSe biomarker analysis results |
| `deseq2_volcano.png` | Volcano plot (log2 FC vs -log10 adjusted p-value) |
| `deseq2_ma.png` | MA plot (mean abundance vs log2 FC) |
| `lefse_cladogram.png` | LEfSe biomarker cladogram |
| `sig_heatmap.png` | Z-score heatmap of significant taxa |

**Note:** With only 2 samples (one per condition group), p-values from DESeq2, ANCOM-BC, and LEfSe are not statistically meaningful. These outputs demonstrate the pipeline's analytical capability and are best interpreted when ≥3 samples per group are available.

### 11. MultiQC Report — `results/multiqc/`

| File | Description |
|------|-------------|
| `multiqc_report.html` | Aggregated quality report across all pipeline modules |

---

## Key Biological Observations

1. **Microbiome profile** is typical of upper GI / duodenal origin with *Firmicutes* (esp. *Streptococcus*) dominating.
2. **Low diversity** (Shannon ~3.5) relative to fecal samples (typically Shannon 5–7), consistent with the lower biomass and selective environment of the duodenum.
3. **Differential abundance** is exploratory with only 2 samples; additional biological replicates are needed for statistically robust comparisons between groups.
4. **Assembly and binning** are limited by the 16S amplicon strategy — whole-genome shotgun data would enable robust MAG recovery.
