process PRODIGAL {
    tag    "${meta.id}"
    label  'prodigal'
    publishDir "${params.outdir}/annotation/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.faa"),  emit: proteins
    tuple val(meta), path("*.fna"),  emit: nucleotide
    tuple val(meta), path("*.gff"),  emit: gff

    script:
    """
    prodigal -i ${assembly}                              \
             -a ${meta.id}.proteins.faa                  \
             -d ${meta.id}.genes.fna                     \
             -f gff -o ${meta.id}.annotations.gff        \
             -p meta

    # Simple stats
    echo "Total proteins: \$(grep -c '^>' ${meta.id}.proteins.faa)" > prodigal_stats.txt
    """

    stub:
    """
    touch ${meta.id}.proteins.faa ${meta.id}.genes.fna ${meta.id}.annotations.gff prodigal_stats.txt
    """
}
