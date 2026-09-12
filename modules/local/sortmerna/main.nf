// rRNA depletion: strips rRNA-aligned reads from trimmed pairs.
process SORTMERNA {
    tag "${meta.id}"
    label 'qc'

    input:
    tuple val(meta), path(r1), path(r2)
    path ref
    path idx

    output:
    // '--other X.clean' is a *prefix*: SortMeRNA appends the extension and
    // writes X.clean_fwd.fq.gz / X.clean_rev.fq.gz. Renaming the output:
    // declarations produces a missing-output error, not a rename.
    tuple val(meta), path("${meta.id}.clean_fwd.fq.gz"),
                     path("${meta.id}.clean_rev.fq.gz"), emit: clean
    tuple val(meta), path("${meta.id}.rrna_fwd.fq.gz"),
                     path("${meta.id}.rrna_rev.fq.gz"),  emit: rrna
    path "${meta.id}.rrna.log",                          emit: log

    script:
    """
    # '--index 0' reuses the shared index and errors out if it is missing,
    # rather than silently rebuilding it once per sample.
    # '--paired_out' keeps a pair together in the non-rRNA output when only one
    # mate aligns, so R1/R2 stay in lockstep for the assembler.
    sortmerna \\
        --ref ${ref} \\
        --reads ${r1} \\
        --reads ${r2} \\
        --workdir \$PWD/smr_work \\
        --idx-dir ${idx} \\
        --index 0 \\
        --aligned \$PWD/${meta.id}.rrna \\
        --other \$PWD/${meta.id}.clean \\
        --fastx \\
        --paired_out \\
        --out2 \\
        --threads ${task.cpus}
    """
}