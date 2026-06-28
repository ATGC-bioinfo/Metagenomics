#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# Verify all required tools and databases are available
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

echo "=== Checking Tool Dependencies ==="

tools=(
    "fastqc"
    "fastp"
    "megahit"
    "bwa"
    "metabat2"
    "jgi_summarize_bam_contig_depths"
    "checkm2"
    "kraken2"
    "bracken"
    "prodigal"
    "multiqc"
)

missing=0
for tool in "${tools[@]}"; do
    if command -v "$tool" &>/dev/null; then
        echo "  [✓] $tool"
    else
        echo "  [✗] $tool -- NOT FOUND"
        missing=$((missing + 1))
    fi
done

echo ""
echo "=== Database Paths ==="
echo "  Kraken2 DB : ${params_kraken_db:-NOT SET}"
echo "  CheckM2 DB : ${params_checkm_db:-NOT SET}"

if [ $missing -gt 0 ]; then
    echo ""
    echo "WARNING: $missing tool(s) missing. Install via conda/mamba:"
    echo "  conda install -c bioconda fastqc fastp megahit bwa metabat2 checkm2 kraken2 bracken prodigal multiqc"
    exit 1
fi

echo ""
echo "All dependencies satisfied."
