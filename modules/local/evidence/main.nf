// Merge geNomad and bowtie2 map-back results into one per-contig evidence table.
process EVIDENCE {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(contigs), path(virus_summary), path(coverage)

    output:
    tuple val(meta), path("${meta.id}.evidence.tsv"), emit: evidence

    script:
    """
    python3 ${projectDir}/bin/aggregate_evidence.py ${contigs} ${virus_summary} ${coverage} > ${meta.id}.evidence.tsv
    """
}
