process MEGAHIT {
    tag    "${meta.id}"
    label  'megahit'
    publishDir "${params.outdir}/assembly/${meta.id}", mode: 'copy', pattern: "${meta.id}.contigs.fa"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa"), emit: assembly
    path("megahit.log"),                           emit: log

    script:
    def min_len = params.min_contig_len ? "--min-contig-len ${params.min_contig_len}" : ''
    if (meta.single_end) {
        """
        rm -rf megahit_out
        megahit -r ${reads}                               \
                -o megahit_out                             \
                --num-cpu-threads ${task.cpus}             \
                --memory ${params.assembly_mem}            \
                ${min_len}                                 \
                2> megahit.log
        mv megahit_out/final.contigs.fa ${meta.id}.contigs.fa
        """
    } else {
        def (r1, r2) = reads
        """
        rm -rf megahit_out
        megahit -1 ${r1} -2 ${r2}                         \
                -o megahit_out                             \
                --num-cpu-threads ${task.cpus}             \
                --memory ${params.assembly_mem}            \
                ${min_len}                                 \
                2> megahit.log
        mv megahit_out/final.contigs.fa ${meta.id}.contigs.fa
        """
    }

    stub:
    """
    touch ${meta.id}.contigs.fa megahit.log
    """
}
