process HUMANN {
    tag    "${meta.id}"
    label  'humann'
    publishDir "${params.outdir}/functional/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_*"), emit: results
    path("humann.log"),                    emit: log

    script:
    def (r1, r2) = reads
    """
    humann \\
        --input ${r1} \\
        --input ${r2} \\
        --output . \\
        --threads ${task.cpus} \\
        --metaphlan ${params.kraken_db}/metaphlan \\
        --bowtie-db ${params.humann_db}/chocophlan \\
        --protein-db ${params.humann_db}/uniref \\
        --diamond-db ${params.humann_db}/diamond \\
        2> humann.log
    """

    stub:
    """
    touch ${meta.id}_genefamilies.tsv ${meta.id}_pathabundance.tsv humann.log
    """
}
