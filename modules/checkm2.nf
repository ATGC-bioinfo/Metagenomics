process CHECKM2 {
    tag    "${meta.id}"
    label  'checkm2'
    publishDir "${params.outdir}/checkm2/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("quality_report.tsv"), emit: quality_report
    path("checkm2.log"),                        emit: log

    script:
    """
    checkm2 predict --threads ${task.cpus}        \
                    --input bins                   \
                    --output-directory ./          \
                    --db-path ${params.checkm_db}  \
                    --extension .fa                \
                    2> checkm2.log
    """

    stub:
    """
    touch quality_report.tsv checkm2.log
    """
}
