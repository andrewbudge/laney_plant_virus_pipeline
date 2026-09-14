// Classify filtered contigs with geNomad. Plasmid_* outputs also exist but are
// not declared: V1 is virus-focused only.
process GENOMAD {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(contigs)
    path db

    output:
    tuple val(meta), path("${meta.id}_genomad/${meta.id}_summary/${meta.id}_virus_summary.tsv"), emit: virus_summary
    tuple val(meta), path("${meta.id}_genomad/${meta.id}_summary/${meta.id}_virus.fna"), emit: virus_contigs
    tuple val(meta), path("${meta.id}_genomad/${meta.id}_summary/${meta.id}_virus_proteins.faa"), emit: virus_proteins
    path "${meta.id}_genomad/${meta.id}_summary/${meta.id}_summary.json", emit: summary_json
    path "${meta.id}.genomad.log", emit: log

    script:
    """
    cp ${contigs} ${meta.id}.fna
    genomad end-to-end --cleanup --threads ${task.cpus} ${meta.id}.fna ${meta.id}_genomad ${db} 2> ${meta.id}.genomad.log
    """
}
