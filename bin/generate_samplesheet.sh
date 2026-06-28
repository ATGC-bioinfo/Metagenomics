#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# Generate samplesheet.csv from a directory of FASTQ files
# Usage: generate_samplesheet.sh <fastq_dir> > samplesheet.csv
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <fastq_dir>" >&2
    exit 1
fi

fastq_dir="$1"
echo "sample_id,read1,read2,single_end"

# Collect all sample prefixes by stripping _R1/_R2 or .fastq suffixes
for f in "$fastq_dir"/*.fastq.gz; do
    base=$(basename "$f" .fastq.gz)

    # Strip _R1 / _R2 to get sample ID
    sample="${base%_R1}"
    sample="${sample%_R2}"

    # Skip if we already processed this sample
    if [[ "$base" == *_R1.fastq.gz ]]; then
        r1="$f"
        r2="${f/_R1/_R2}"
        if [ -f "$r2" ]; then
            echo "${sample},${r1},${r2},false"
        else
            echo "${sample},${f},,true"
        fi
    elif [[ "$base" != *_R2* ]]; then
        echo "${sample},${f},,true"
    fi
done | sort -u
