// Filter SPAdes contigs by length and mean coverage, and prefix every header
// with the sample id. One pass: seqkit fx2tab then awk parses the
// 'NODE_1_length_20085_cov_12194.803368' header and re-emits FASTA.
process FILTER {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.contigs.fasta"), emit: contigs

    script:
    """
    seqkit fx2tab ${contigs} \\
      | awk -F'\\t' -v minlen=${params.contig_min_length} -v mincov=${params.contig_min_cov} -v id='${meta.id}' '
          { split(\$1, a, "_"); if (a[4] >= minlen && a[6] >= mincov) print ">" id "_" substr(\$1, 6) "\\n" \$2 }' \\
      > ${meta.id}.contigs.fasta
    """
}