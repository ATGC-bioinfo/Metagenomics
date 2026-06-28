process COMMUNITY {
    tag    "CommunityAnalysis"
    label  'community'
    publishDir "${params.outdir}/community", mode: 'copy', pattern: "*.png"

    input:
    path(bracken_reports)

    output:
    path("community_alpha.png"),    emit: alpha
    path("community_beta.png"),     emit: beta
    path("community_heatmap.png"),  emit: heatmap
    path("community_pca.png"),      emit: pca

    script:
    """
    python3 ${projectDir}/bin/community_analysis.py
    """

    stub:
    """
    touch community_alpha.png community_beta.png community_heatmap.png community_pca.png
    """
}
