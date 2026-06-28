process PLOTS {
    tag    "${meta.id}"
    label  'plots'
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy', pattern: "*.png"

    input:
    tuple val(meta), path(bracken_p), path(bracken_g), path(bracken_s)

    output:
    tuple val(meta), path("${meta.id}.sankey.png"),  emit: sankey
    tuple val(meta), path("${meta.id}.summary.png"),  emit: summary

    script:
    """
    python3 ${projectDir}/bin/per_sample_plots.py "${meta.id}" "${bracken_p}" "${bracken_g}" "${bracken_s}"
    """

    stub:
    """
    touch ${meta.id}.sankey.png ${meta.id}.summary.png
    """
}
