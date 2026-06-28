process DIFFABUND {
    tag    "DifferentialAbundance"
    label  'diffabund'
    publishDir "${params.outdir}/diffabund", mode: 'copy', pattern: "*.{tsv,png}"

    input:
    path(metadata_file)
    path(bracken_reports)

    output:
    path("deseq2_results.tsv"),       emit: deseq2_table
    path("ancombc_results.tsv"),      emit: ancombc_table
    path("lefse_results.tsv"),        emit: lefse_table
    path("deseq2_volcano.png"),      emit: volcano_plot
    path("deseq2_ma.png"),           emit: ma_plot
    path("lefse_cladogram.png"),     emit: cladogram
    path("sig_heatmap.png"),         emit: heatmap

    script:
    """
    python3 ${projectDir}/bin/diff_abundance.py "${metadata_file}"
    """

    stub:
    """
    touch deseq2_results.tsv ancombc_results.tsv lefse_results.tsv \
          deseq2_volcano.png deseq2_ma.png \
          lefse_cladogram.png sig_heatmap.png
    """
}
