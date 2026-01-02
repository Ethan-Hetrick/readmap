process CUSTOMSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker://staphb/pandas:2.3.3'

    input:

    tuple val(meta), path(stats)
    tuple val(meta2), path(coverage)

    output:
    tuple val(meta), path("*_customstats.txt"), emit: customstats
    tuple val("${task.process}"), val('customstats'), eval("python3 --version"), topic: versions, emit: versions_customstats

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sam_stats.py $stats $coverage > ${prefix}_customstats.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    
    touch ${prefix}.txt
    """
}
