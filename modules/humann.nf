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
    if (meta.single_end) {
        """
        humann \\
            --input ${reads} \\
            --output . \\
            --threads ${task.cpus} \\
            --metaphlan ${params.humann_db}/metaphlan \\
            --bowtie-db ${params.humann_db}/chocophlan \\
            --protein-db ${params.humann_db}/uniref \\
            2> humann.log
        """
    } else {
        def (r1, r2) = reads
        """
        cat ${r1} ${r2} > ${meta.id}.merged.fq
        humann \\
            --input ${meta.id}.merged.fq \\
            --output . \\
            --threads ${task.cpus} \\
            --metaphlan ${params.humann_db}/metaphlan \\
            --bowtie-db ${params.humann_db}/chocophlan \\
            --protein-db ${params.humann_db}/uniref \\
            2> humann.log
        """
    }

    stub:
    """
    touch ${meta.id}_genefamilies.tsv ${meta.id}_pathabundance.tsv humann.log
    """
}
