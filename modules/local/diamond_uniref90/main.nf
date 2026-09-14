process DIAMOND_UNIREF90 {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(contigs)
    path db

    output:
    tuple val(meta), path("${meta.id}.uniref90.tsv"), emit: uniref90

    script:
    """
    diamond blastx -d ${db} -q ${contigs} -o ${meta.id}.uniref90.tsv --more-sensitive --evalue 1e-5 --max-target-seqs 5 --query-cover 70 --threads ${task.cpus} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp
    """
}
