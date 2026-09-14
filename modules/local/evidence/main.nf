// Merge geNomad and bowtie2 map-back results into one per-contig evidence table with R.
process EVIDENCE {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(contigs), path(virus_summary), path(coverage)

    output:
    tuple val(meta), path("${meta.id}.evidence.tsv"), emit: evidence

    script:
    """
    Rscript ${projectDir}/bin/aggregate_evidence.R ${contigs} ${virus_summary} ${coverage} > ${meta.id}.evidence.tsv
    """
}
