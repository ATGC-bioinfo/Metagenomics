process MULTIQC {
    tag    "MultiQC"
    label  'multiqc'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(fastqc_html)
    path(fastqc_zips)
    path(fastp_html)
    path(assemblies)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"),        emit: data

    script:
    """
    multiqc .
    """

    stub:
    """
    echo '<html>MultiQC Stub</html>' > multiqc_report.html
    mkdir -p multiqc_data
    """
}
