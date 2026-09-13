// Screen filtered contigs against the viral nucleotide database with blastn.
// The database lives at <db>/blastn/<name>.* (built with makeblastdb
// -dbtype nucl -parse_seqids, e.g. viral_genomic.*); the whole directory is
// staged so the prefix files arrive intact.
process BLASTN {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(contigs)
    path blastdb
    val blastn_db_name

    output:
    tuple val(meta), path("${meta.id}.blastn.tsv"), emit: hits

    script:
    """
    blastn \\
        -db ${blastdb}/${blastn_db_name} \\
        -query ${contigs} \\
        -out ${meta.id}.blastn.tsv \\
        -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs slen sstart send staxids sscinames" \\
        -evalue ${params.blastn_evalue} \\
        -max_target_seqs ${params.blastn_max_target_seqs} \\
        -num_threads ${task.cpus}
    """
}