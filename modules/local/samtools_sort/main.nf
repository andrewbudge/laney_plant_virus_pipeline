// Sort and index bowtie2 alignments for downstream per-contig coverage.
process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam
    tuple val(meta), path("${meta.id}.sorted.bam.bai"), emit: bai
    // Per-contig collapse of the alignment (numreads/covbases/coverage/meandepth) —
    // the joinable unit for the evidence table; raw per-base depth stays recomputable from the BAM.
    tuple val(meta), path("${meta.id}.coverage.tsv"), emit: coverage

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam ${sam}
    samtools index ${meta.id}.sorted.bam
    samtools coverage ${meta.id}.sorted.bam > ${meta.id}.coverage.tsv
    """
}
