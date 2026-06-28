process GTDBTK {
    tag    "${meta.id}"
    label  'gtdbtk'
    publishDir "${params.outdir}/gtdbtk/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("gtdbtk_results"), emit: results
    path("gtdbtk.log"),                     emit: log

    script:
    """
    GTDB-Tk classify_wf \\
        --genome_dir ${bins} \\
        --out_dir gtdbtk_results \\
        --cpus ${task.cpus} \\
        --pplacer_cpus ${task.cpus} \\
        -x fa 2>&1 | tee gtdbtk.log
    """

    stub:
    """
    mkdir -p gtdbtk_results
    touch gtdbtk_results/gtdbtk.summary.tsv gtdbtk.log
    """
}
