// Map rRNA-depleted reads back to filtered contigs for per-contig support.
process BOWTIE2 {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(r1), path(r2), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam
    path "${meta.id}.bowtie2.log", emit: log

    script:
    """
    bowtie2-build ${contigs} ${meta.id}_idx
    # Unaligned reads would dominate the SAM (most reads are host) and are
    # useless to the downstream per-contig coverage step; stderr still reports
    # the total alignment rate.
    bowtie2 -x ${meta.id}_idx -1 ${r1} -2 ${r2} --threads ${task.cpus} --no-unal -S ${meta.id}.sam 2> ${meta.id}.bowtie2.log
    """
}
