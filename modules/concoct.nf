process CONCOCT {
    tag    "${meta.id}"
    label  'concoct'
    publishDir "${params.outdir}/bins/${meta.id}", mode: 'copy', pattern: "concoct_bins/*.fa"

    input:
    tuple val(meta), path(assembly), path(reads)

    output:
    tuple val(meta), path("concoct_bins/*.fa"), emit: bins

    script:
    def (r1, r2) = reads
    """
    # Map reads to assembly
    bwa index ${assembly}
    bwa mem -t ${task.cpus} ${assembly} ${r1} ${r2} > aln.sam
    samtools sort -@ ${task.cpus} aln.sam -o aln.sorted.bam
    rm -f aln.sam

    # Calculate coverage
    mkdir -p concoct_input
    samtools index aln.sorted.bam
    concoct_coverage_table.py concoct_input/coverage.tsv aln.sorted.bam

    # Run CONCOCT
    mkdir -p concoct_bins
    concoct --composition_file ${assembly} \\
            --coverage_file concoct_input/coverage.tsv \\
            --threads ${task.cpus} \\
            -b concoct_output/ 2> concoct.log

    # Extract bins as FASTA
    mkdir -p concoct_bins
    cut_up_fasta.py ${assembly} \\
        -c 10000 -o 0 \\
        --merge_last -b contigs_10K.bed > contigs_10K.fa
    extract_clusters.py contigs_10K.fa concoct_output/clustering_gt.csv > /dev/null 2>&1
    mv concoct_output/clustering_gt.csv/*.fa concoct_bins/ 2>/dev/null || true

    # Cleanup
    rm -f aln.sorted.bam aln.sorted.bam.bai
    """

    stub:
    """
    mkdir -p concoct_bins
    touch concoct_bins/placeholder.fa
    """
}
