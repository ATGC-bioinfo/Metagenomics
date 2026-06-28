process FASTQC {
    tag    "${meta.id}"
    label  'fastqc'
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip"),  emit: zip

    script:
    if (meta.single_end) {
        """
        fastqc -t ${task.cpus} --noextract ${reads}
        """
    } else {
        def (r1, r2) = reads
        """
        fastqc -t ${task.cpus} --noextract ${r1} ${r2}
        """
    }

    stub:
    """
    touch ${meta.id}_fastqc.html ${meta.id}_fastqc.zip
    """
}
