process CHECKM {
    tag    "${meta.id}"
    label  'checkm'
    publishDir "${params.outdir}/checkm/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("checkm_results"), emit: results
    path("checkm.log"),                     emit: log

    script:
    """
    checkm lineage_wf -x fa \\
        --threads ${task.cpus} \\
        --pplacer_threads ${task.cpus} \\
        --tab_table \\
        -f checkm_results/quality_report.tsv \\
        ${bins} checkm_results 2> checkm.log

    checkm qa \\
        -o 2 \\
        --threads ${task.cpus} \\
        checkm_results/lineage.ms checkm_results \\
        --tab_table -f checkm_results/detailed.tsv 2>> checkm.log || true
    """

    stub:
    """
    mkdir -p checkm_results
    touch checkm_results/quality_report.tsv checkm.log
    """
}
