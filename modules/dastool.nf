process DASTOOL {
    tag    "${meta.id}"
    label  'dastool'
    publishDir "${params.outdir}/bins/${meta.id}", mode: 'copy', pattern: "das_tool_bins/*.fa"

    input:
    tuple val(meta), path(assembly), path(binners)

    output:
    tuple val(meta), path("das_tool_bins/*.fa"), emit: bins

    script:
    """
    # Prepare bin sets for DAS Tool
    mkdir -p das_tool_input das_tool_bins

    # Symlink all bins from different binners
    for f in ${binners}; do
        ln -sf "../\$f" das_tool_input/
    done

    # Generate scaffold2bin tables
    Fasta_to_Scaffolds2Bin.sh -e fa \
        -i das_tool_input \
        -o das_tool_input/scaffolds2bin.tsv 2>/dev/null || true

    # Run DAS Tool
    DAS_Tool -i das_tool_input/scaffolds2bin.tsv \
             -l metabat2,concoct,maxbin2 \
             -c ${assembly} \
             -o das_tool_bins/DASTool \
             --search_engine diamond \
             --threads ${task.cpus} \
             --write_bins \
             --score_threshold 0.5 2> dastool.log

    mv das_tool_bins/DASTool_DASTool_bins/*.fa das_tool_bins/ 2>/dev/null || true
    """

    stub:
    """
    mkdir -p das_tool_bins
    touch das_tool_bins/placeholder.fa
    """
}
