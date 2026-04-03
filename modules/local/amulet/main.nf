process AMULET {
    tag "$meta.id"
    container 'quay.io/cellgeni/amulet:1.1'

    input:
    tuple val(meta), path(fragments), path(metrics)
    path autosomes
    path blacklist


    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.csv"), emit: csv
    tuple val("${task.process}"), val('amulet'), eval("echo v1.1"), topic: versions, emit: versions_amulet

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--qthresh 0.01'
    """
    amulet.sh \
        ${fragments} \
        ${metrics} \
        ${autosomes} \
        ${blacklist} \
        . \
        ${args}
    cat MultipletProbabilities.txt | tr '\t' ',' > MultipletProbabilities.csv
    """
}
