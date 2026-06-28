process METABAT2 {
    tag    "${meta.id}"
    label  'metabat2'
    publishDir "${params.outdir}/bins/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly), path(reads)

    output:
    tuple val(meta), path("bins/*.fa"), emit: bins
    path("depth.txt"),                  emit: depth

    script:
    def (r1, r2) = reads
    def min_contig = 1500  // MetaBAT2 hard minimum
    """
    # Count contigs long enough for MetaBAT2 (hard minimum 1500 bp)
    long_count=\$(awk '!/^>/ {len += length(\$0)} /^>/ {if (len >= ${min_contig}) count++; len=0} END {if (len >= ${min_contig}) count++; print count+0}' ${assembly})

    if [ "\$long_count" -eq 0 ]; then
        echo "No contigs >= ${min_contig} bp found. Skipping MetaBAT2 for ${meta.id}."
        mkdir -p bins
        echo "No bins produced (assembly too fragmented)" > bins/bin.skipped.fa
        touch depth.txt
    else
        # Map reads to assembly for abundance profiling
        bwa index ${assembly}
        bwa mem -t ${task.cpus} ${assembly} ${r1} ${r2} > aln.sam

        # Sort SAM into BAM (required by jgi_summarize_bam_contig_depths)
        samtools sort -@ ${task.cpus} aln.sam -o aln.sorted.bam

        # Generate depth file
        jgi_summarize_bam_contig_depths --outputDepth depth.txt aln.sorted.bam

        # Run MetaBAT2
        mkdir -p bins
        metabat2 -i ${assembly}                           \\
                 -a depth.txt                             \\
                 -o bins/bin                              \\
                 -t ${task.cpus}                          \\
                 -m ${params.metabat_min_contig_len}       \\
                 --unbinned

        # Clean up large intermediates
        rm -f aln.sam aln.sorted.bam ${assembly}.amb ${assembly}.ann ${assembly}.bwt ${assembly}.pac ${assembly}.sa
    fi
    """

    stub:
    """
    mkdir -p bins
    touch bins/bin.1.fa depth.txt
    """
}
