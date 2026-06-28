process KRAKEN2_BRACKEN {
    tag    "${meta.id}"
    label  'kraken2'
    publishDir "${params.outdir}/taxonomy/${meta.id}", mode: 'copy', pattern: "*.{report,output,classified*,unclassified*,csv}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.kraken2.report"),        emit: report
    tuple val(meta), path("*.kraken2.output"),        emit: kraken_output
    tuple val(meta), path("*.bracken.S.report"),      emit: bracken_s
    tuple val(meta), path("*.bracken.G.report"),      emit: bracken_g
    tuple val(meta), path("*.bracken.P.report"),      emit: bracken_p
    tuple val(meta), path("*.kraken2.csv"),           emit: kraken_csv
    tuple val(meta), path("*.bracken.S.csv"),         emit: bracken_s_csv
    tuple val(meta), path("*.bracken.G.csv"),         emit: bracken_g_csv
    tuple val(meta), path("*.bracken.P.csv"),         emit: bracken_p_csv
    tuple val(meta), path("*.classified*"),           emit: classified
    tuple val(meta), path("*.unclassified*"),         emit: unclassified

    script:
    def read_len = params.read_length ?: 150
    if (meta.single_end) {
        """
        kraken2 --db ${params.kraken_db}                 \
                --threads ${task.cpus}                    \
                --report ${meta.id}.kraken2.report        \
                --output ${meta.id}.kraken2.output        \
                --classified-out ${meta.id}.classified.fq \
                --unclassified-out ${meta.id}.unclassified.fq \
                ${reads}

        # Bracken: abundance estimation at multiple ranks
        bracken -d ${params.kraken_db}                    \
                -i ${meta.id}.kraken2.report              \
                -o ${meta.id}.bracken.S.report            \
                -r ${read_len} -l S

        bracken -d ${params.kraken_db}                    \
                -i ${meta.id}.kraken2.report              \
                -o ${meta.id}.bracken.G.report            \
                -r ${read_len} -l G

        bracken -d ${params.kraken_db}                    \
                -i ${meta.id}.kraken2.report              \
                -o ${meta.id}.bracken.P.report            \
                -r ${read_len} -l P

        python3 ${projectDir}/bin/report_to_csv.py "${meta.id}"
        """
    } else {
        def (r1, r2) = reads
        """
        kraken2 --db ${params.kraken_db}                  \
                --threads ${task.cpus}                     \
                --report ${meta.id}.kraken2.report         \
                --output ${meta.id}.kraken2.output         \
                --classified-out ${meta.id}.classified#.fq \
                --unclassified-out ${meta.id}.unclassified#.fq \
                --paired ${r1} ${r2}

        bracken -d ${params.kraken_db}                    \
                -i ${meta.id}.kraken2.report              \
                -o ${meta.id}.bracken.S.report            \
                -r ${read_len} -l S

        bracken -d ${params.kraken_db}                    \
                -i ${meta.id}.kraken2.report              \
                -o ${meta.id}.bracken.G.report            \
                -r ${read_len} -l G

        bracken -d ${params.kraken_db}                    \
                -i ${meta.id}.kraken2.report              \
                -o ${meta.id}.bracken.P.report            \
                -r ${read_len} -l P

        python3 ${projectDir}/bin/report_to_csv.py "${meta.id}"
        """
    }

    stub:
    """
    touch ${meta.id}.kraken2.report ${meta.id}.kraken2.output \
          ${meta.id}.bracken.S.report ${meta.id}.bracken.G.report ${meta.id}.bracken.P.report \
          ${meta.id}.kraken2.csv \
          ${meta.id}.bracken.S.csv ${meta.id}.bracken.G.csv ${meta.id}.bracken.P.csv \
          ${meta.id}.classified.fq ${meta.id}.unclassified.fq
    """
}
