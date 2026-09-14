// Merge all evidence legs into one descriptive per-sample TSV.
process EVIDENCE {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(rvdb), path(uniref90), path(virus_summary), path(coverage)
    path rvdb_taxmap
    path uniref90_taxmap
    path aggregate_script

    output:
    tuple val(meta), path("${meta.id}.evidence.tsv"), emit: evidence

    script:
    """
    mkdir -p inputs/05_genomad/${meta.id} inputs/06_bowtie2/${meta.id} inputs/07_diamond/${meta.id} db/blastx
    cp ${virus_summary} inputs/05_genomad/${meta.id}/${meta.id}_virus_summary.tsv
    cp ${coverage} inputs/06_bowtie2/${meta.id}/${meta.id}.coverage.tsv
    cp ${rvdb} inputs/07_diamond/${meta.id}/${meta.id}.rvdb.tsv
    cp ${uniref90} inputs/07_diamond/${meta.id}/${meta.id}.uniref90.tsv
    ln -s ../../${rvdb_taxmap} db/blastx/U-RVDBv32.0-prot.taxmap.tsv
    ln -s ../../${uniref90_taxmap} db/blastx/uniref90.taxmap.tsv
    cp ${aggregate_script} aggregate_evidence.sh
    chmod +x aggregate_evidence.sh
    ./aggregate_evidence.sh -s ${meta.id} -o inputs -d db
    mv inputs/evidence/${meta.id}.evidence.tsv ${meta.id}.evidence.tsv
    """
}
