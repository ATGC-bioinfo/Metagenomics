process MAXBIN2 {
    tag    "${meta.id}"
    label  'maxbin2'
    publishDir "${params.outdir}/bins/${meta.id}", mode: 'copy', pattern: "maxbin2_bins/*.fasta"

    input:
    tuple val(meta), path(assembly), path(reads)

    output:
    tuple val(meta), path("maxbin2_bins/*.fasta"), emit: bins

    script:
    def (r1, r2) = reads
    """
    # Map reads
    bwa index ${assembly}
    bwa mem -t ${task.cpus} ${assembly} ${r1} ${r2} > aln.sam
    samtools sort -@ ${task.cpus} aln.sam -o aln.sorted.bam
    rm -f aln.sam

    # Calculate read depths
    mkdir -p maxbin2_input
    samtools depth aln.sorted.bam > maxbin2_input/depth.txt

    # Run MaxBin2
    mkdir -p maxbin2_bins
    run_MaxBin.pl -contig ${assembly} \\
                  -out maxbin2_bins/bin \\
                  -abund maxbin2_input/depth.txt \\
                  -thread ${task.cpus} 2> maxbin2.log

    rm -f aln.sorted.bam aln.sorted.bam.bai
    """

    stub:
    """
    mkdir -p maxbin2_bins
    touch maxbin2_bins/placeholder.fasta
    """
}
