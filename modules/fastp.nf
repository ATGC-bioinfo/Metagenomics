process FASTP {
    tag    "${meta.id}"
    label  'fastp'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"),             emit: json
    tuple val(meta), path("*.html"),             emit: html

    script:
    if (meta.single_end) {
        """
        fastp -i ${reads} \
              -o ${meta.id}.trimmed.fastq.gz \
              --json ${meta.id}_fastp.json \
              --html ${meta.id}_fastp.html \
              --qualified_quality_phred 20 \
              --length_required 50 \
              --cut_front --cut_tail \
              --cut_mean_quality 20 \
              --thread ${task.cpus}
        """
    } else {
        def (r1, r2) = reads
        """
        fastp -i ${r1} -I ${r2} \
              -o ${meta.id}_R1.trimmed.fastq.gz \
              -O ${meta.id}_R2.trimmed.fastq.gz \
              --json ${meta.id}_fastp.json \
              --html ${meta.id}_fastp.html \
              --qualified_quality_phred 20 \
              --length_required 50 \
              --cut_front --cut_tail \
              --cut_mean_quality 20 \
              --detect_adapter_for_pe \
              --thread ${task.cpus}
        """
    }

    stub:
    """
    touch ${meta.id}_R1.trimmed.fastq.gz ${meta.id}_fastp.json ${meta.id}_fastp.html
    """
}
