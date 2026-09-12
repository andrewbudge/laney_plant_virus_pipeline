#!/usr/bin/env nextflow
//
// plant-virus QC — adapter/quality trim, rRNA depletion, aggregate report.
//
//   fastp ──► sortmerna ──► (clean paired reads, for downstream assembly)
//     └──────────┴────────► multiqc
//
//   nextflow run main.nf -profile slurm --input samplesheet.csv --outdir <dir>
//
nextflow.enable.dsl = 2

// SortMeRNA is pinned to 4.3.4 and must not be bumped: 4.3.6, 4.3.7 and the
// 7.0.0 binary are built for x86-64-v3 (AVX2+BMI2+FMA) and die with SIGILL on
// lonepeak's Sandy Bridge nodes. 4.3.4 is the last baseline-SSE build. Bumping
// it means moving SORTMERNA to notchpeak or granite — see conf/base.config.
def FASTP_IMAGE     = 'docker://quay.io/biocontainers/fastp:0.23.4--h125f33a_5'
def SORTMERNA_IMAGE = 'docker://quay.io/biocontainers/sortmerna:4.3.4--h9ee0642_0'
def MULTIQC_IMAGE   = 'docker://quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0'
def SPADES_IMAGE    = 'docker://quay.io/biocontainers/spades:4.3.0--hde4eca7_1'

process FASTP {
    tag "$sample"
    container FASTP_IMAGE
    publishDir "${params.outdir}/fastp/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_1.fastp.fastq.gz"),
                       path("${sample}_2.fastp.fastq.gz"), emit: reads
    path "${sample}.fastp.json",                           emit: json
    path "${sample}.fastp.html"

    script:
    """
    fastp \\
        -i ${r1} -I ${r2} \\
        -o ${sample}_1.fastp.fastq.gz -O ${sample}_2.fastp.fastq.gz \\
        --detect_adapter_for_pe \\
        -q ${params.fastp_quality} \\
        -h ${sample}.fastp.html \\
        -j ${sample}.fastp.json \\
        -R "${sample} fastp report" \\
        -w ${task.cpus}
    """
}


// Built once and shared by every sample; skipped when --sortmerna_idx is set.
process SORTMERNA_INDEX {
    tag "${ref.name}"
    container SORTMERNA_IMAGE
    publishDir "${params.outdir}/sortmerna_index", mode: 'copy'

    input:
    path ref

    output:
    path "idx", emit: idx

    script:
    """
    mkdir -p idx
    # '--index 1' means index and exit. It still prints "Missing required flag:
    # reads" before doing so — that message is not an error here.
    sortmerna \\
        --ref ${ref} \\
        --idx-dir \$PWD/idx \\
        --workdir \$PWD/index_work \\
        --index 1 \\
        -m ${params.sortmerna_index_mem}
    """
}


process SORTMERNA {
    tag "$sample"
    container SORTMERNA_IMAGE
    publishDir "${params.outdir}/sortmerna/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)
    path ref
    path idx

    output:
    // '--other X.clean' is a *prefix*: SortMeRNA appends the extension and
    // writes X.clean_fwd.fq.gz / X.clean_rev.fq.gz. Renaming these gives a
    // missing-output error, not a renamed file.
    tuple val(sample), path("${sample}.clean_fwd.fq.gz"),
                       path("${sample}.clean_rev.fq.gz"), emit: clean
    tuple val(sample), path("${sample}.rrna_fwd.fq.gz"),
                       path("${sample}.rrna_rev.fq.gz"),  emit: rrna
    path "${sample}.rrna.log",                            emit: log

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
        --aligned \$PWD/${sample}.rrna \\
        --other \$PWD/${sample}.clean \\
        --fastx \\
        --paired_out \\
        --out2 \\
        --threads ${task.cpus}
    """
}



process MULTIQC {
    container MULTIQC_IMAGE
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path report_files

    output:
    path "multiqc_report.html"
    path "multiqc_report_data"

    script:
    """
    multiqc --force -n multiqc_report.html .
    """
}

// TODO: assembly. Commented out rather than deleted -- a process with no
// script block fails compilation, so the stub blocked every run while it sat
// here uncommented. Fill in input/output/script and uncomment.
//
// process SPADESASSEMBLY {
//     container SPADES_IMAGE
// }

// Absolute paths are taken as given; relative ones resolve against a root.
// Relative-to-launch-dir is deliberately not supported: a samplesheet whose
// meaning depends on the caller's cwd is not reproducible, and it is the one
// failure mode that looks like it works until someone cds somewhere else.
def resolveIn(root, rootFlag, path, what) {
    if (!path) error "${what} is not set"
    def p = path.toString().trim()
    if (p.startsWith('/')) return file(p, checkIfExists: true)
    if (!root)
        error "${what} is relative ('${p}') but ${rootFlag} is not set — pass ${rootFlag}, or use an absolute path"
    return file("${root.toString().replaceAll('/$', '')}/${p}", checkIfExists: true)
}

// Reads live with the project. References are shared across projects and are
// far too large to copy per project, so they get their own root -- which falls
// back to --data_dir for the common case where both sit under one directory.
def resolveData(path, what) { resolveIn(params.data_dir, '--data_dir', path, what) }
def resolveDb(path, what) {
    resolveIn(params.db_dir ?: params.data_dir, '--db_dir (or --data_dir)', path, what)
}


workflow {

    if (params.data_dir && !file(params.data_dir).exists())
        error "--data_dir does not exist: ${params.data_dir}"
    if (params.db_dir && !file(params.db_dir).exists())
        error "--db_dir does not exist: ${params.db_dir}"

    // Validation sits in the data path, not on a side branch: an operator chain
    // with no consumer never runs, so a dangling validation branch would pass
    // silently. toList/flatMap forces it before FASTP sees a row.
    ch_samples = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2) {
                error "samplesheet ${params.input}: every row needs sample, fastq_1, fastq_2 — got: ${row}"
            }
            tuple(row.sample,
                  resolveData(row.fastq_1, "${row.sample}.fastq_1"),
                  resolveData(row.fastq_2, "${row.sample}.fastq_2"))
        }
        .toList()
        .map { rows ->
            if (!rows) error "samplesheet ${params.input} has no data rows"
            // Duplicates would collide in publishDir and overwrite each other.
            def dupes = rows.collect { it[0] }.countBy { it }.findAll { it.value > 1 }.keySet()
            if (dupes) error "duplicate sample names in ${params.input}: ${dupes}"
            rows
        }
        .flatMap()

    ch_ref = Channel.value(resolveDb(params.sortmerna_ref, '--sortmerna_ref'))

    // Value channel so the one index broadcasts to every sample instead of
    // being consumed by the first one.
    ch_idx = params.sortmerna_idx
        ? Channel.value(resolveDb(params.sortmerna_idx, '--sortmerna_idx'))
        : SORTMERNA_INDEX(ch_ref).idx.first()

    FASTP(ch_samples)
    SORTMERNA(FASTP.out.reads, ch_ref, ch_idx)

    MULTIQC(
        FASTP.out.json.mix(SORTMERNA.out.log).collect()
    )
}
