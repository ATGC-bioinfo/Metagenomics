process QUAST {
    tag    "${meta.id}"
    label  'quast'
    publishDir "${params.outdir}/assembly/${meta.id}", mode: 'copy', pattern: "quast*"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("quast_results"), emit: results
    path("quast.log"),                     emit: log

    script:
    """
    quast.py ${assembly} -o quast_results --threads ${task.cpus} --meta 2> quast.log
    """

    stub:
    """
    mkdir -p quast_results
    touch quast.log
    """
}
