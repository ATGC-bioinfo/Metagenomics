process KRONA {
    tag    "${meta.id}"
    label  'krona'
    publishDir "${params.outdir}/taxonomy/${meta.id}", mode: 'copy', pattern: "*.html"

    input:
    tuple val(meta), path(kraken_output)

    output:
    tuple val(meta), path("${meta.id}.krona.html"), emit: html

    script:
    """
    # Extract taxonomy IDs of classified reads for Krona chart
    awk '\$3 > 0 {print \$3}' ${kraken_output} > taxa.txt
    ktImportTaxonomy -t 1 \\
        -o ${meta.id}.krona.html \\
        -n "${meta.id}" taxa.txt 2>&1
    rm -f taxa.txt
    """

    stub:
    """
    touch ${meta.id}.krona.html
    """
}
