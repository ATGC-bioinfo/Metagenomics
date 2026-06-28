process EGGNOG {
    tag    "${meta.id}"
    label  'eggnog'
    publishDir "${params.outdir}/annotation/${meta.id}", mode: 'copy', pattern: "*.emapper.*"

    input:
    tuple val(meta), path(proteins)

    output:
    tuple val(meta), path("${meta.id}.emapper.annotations"), emit: annotations
    tuple val(meta), path("${meta.id}.emapper.hits"),       emit: hits
    tuple val(meta), path("${meta.id}.emapper.seed_orthologs"), emit: seed

    script:
    """
    emapper.py -i ${proteins}                               \\
               -o ${meta.id}                                \\
               --output_dir .                               \\
               --data_dir ${params.eggnog_db}               \\
               --cpu ${task.cpus}                           \\
               --tax_scope auto                             \\
               --target_orthologs all                       \\
               --go_evidence non-electronic                 \\
               --pfam_realign realign                       \\
               --report_orthologs                           \\
                --decorate_gff ${meta.id}.decorated.gff      \\
                2> ${meta.id}.eggnog.log

    # Summary stats
    grep -c '^#' ${meta.id}.emapper.annotations > ${meta.id}.eggnog_stats.txt 2>/dev/null || echo 0 > ${meta.id}.eggnog_stats.txt
    """

    stub:
    """
    touch ${meta.id}.emapper.annotations ${meta.id}.emapper.hits ${meta.id}.emapper.seed_orthologs ${meta.id}.eggnog.log
    """
}
