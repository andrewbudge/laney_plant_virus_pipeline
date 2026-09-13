// De novo RNA viral assembly from rRNA-depleted reads (rnaviralSPAdes).
// Output lands in <outdir>/spades/<sample>/ per sample; 'spades_out' inside the
// work dir keeps the fixed SPAdes file names apart from the inputs.
process SPADES {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    // rnaviralSPAdes writes genome-style outputs (contigs.fasta +
    // scaffolds.fasta), verified empirically -- the rnaSPAdes-family
    // 'transcripts.fasta' name does NOT appear here.
    tuple val(meta), path("spades_out/contigs.fasta"),   emit: contigs
    tuple val(meta), path("spades_out/scaffolds.fasta"), emit: scaffolds
    tuple val(meta), path("spades_out/spades.log"),      emit: log

    script:
    """
    spades.py \\
        --rnaviral \\
        -1 ${r1} \\
        -2 ${r2} \\
        -o spades_out \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()}
    """
}
