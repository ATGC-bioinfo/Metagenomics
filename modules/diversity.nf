process DIVERSITY {
    tag    "DiversityAnalysis"
    label  'diversity'
    publishDir "${params.outdir}/diversity", mode: 'copy', pattern: "*.{png,tsv}"

    input:
    path(bracken_reports)

    output:
    path("alpha_diversity.tsv"),     emit: alpha_table
    path("beta_diversity.tsv"),      emit: beta_table
    path("alpha_*.png"),             emit: alpha_plots
    path("beta_*.png"),              emit: beta_plots
    path("ordination_*.png"),        emit: ordination_plots

    script:
    """
    python3 ${projectDir}/bin/diversity_analysis.py
    """

    stub:
    """
    touch alpha_diversity.tsv beta_diversity.tsv \
          alpha_boxplots.png alpha_violin.png \
          beta_bray_curtis.png beta_jaccard.png beta_dendrogram.png \
          ordination_pca.png ordination_pcoa.png ordination_nmds.png \
          ordination_umap.png ordination_tsne.png
    """
}
