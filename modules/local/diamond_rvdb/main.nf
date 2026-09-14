process DIAMOND_RVDB {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(contigs)
    path db

    output:
    tuple val(meta), path("${meta.id}.rvdb.tsv"), emit: rvdb

    script:
    """
    diamond blastx -d ${db} -q ${contigs} -o ${meta.id}.rvdb.tsv --more-sensitive --evalue 1e-5 --top 5 --query-cover 70 --threads ${task.cpus}
    """
}
