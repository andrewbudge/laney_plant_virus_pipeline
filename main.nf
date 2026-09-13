#!/usr/bin/env nextflow
//
// plant-virus QC + assembly + viral screening: adapter/quality trim (fastp),
// rRNA depletion (SortMeRNA), de novo RNA viral assembly (SPAdes --rnaviral),
// contig filtering, blastn against the viral nucleotide database, plus a
// MultiQC report per run.
//
// Processes live in modules/local/; their containers, publishDir and tool
// arguments are configured in conf/modules.config. Resources/retry live in
// conf/base.config, the HPC executor in conf/slurm.config.
//
//   nextflow run main.nf -profile local --input samplesheet.tsv --data_dir <dir> --db <dir>
//   nextflow run main.nf -profile slurm --input samplesheet.tsv --data_dir <dir> --db <dir>
//
include { FASTP          } from './modules/local/fastp/main.nf'
include { SORTMERNA      } from './modules/local/sortmerna/main.nf'
include { SPADES         } from './modules/local/spades/main.nf'
include { FILTER         } from './modules/local/filter/main.nf'
include { BLASTN         } from './modules/local/blastn/main.nf'
include { MULTIQC        } from './modules/local/multiqc/main.nf'

// Relative paths resolve against a caller-supplied root; absolute paths pass
// through untouched.
def resolve(path, root, what) {
    if (!path) error "${what} is not set"
    def p = path.toString().trim()
    if (p.startsWith('/')) return file(p, checkIfExists: true)
    if (!root) error "${what} is relative ('${p}') but its root is not set"
    return file("${root.toString().replaceAll('/$', '')}/${p}", checkIfExists: true)
}

workflow {
    if (params.data_dir && !file(params.data_dir).exists())
        error "--data_dir does not exist: ${params.data_dir}"

    db = resolve(params.db, params.data_dir, '--db')
    if (!db.isDirectory()) error "--db is not a directory: ${db}"

    sortmerna_ref = file("${db}/sortmerna/smr_v4.3_default_db.fasta", checkIfExists: true)
    sortmerna_idx = file("${db}/sortmerna/index", checkIfExists: true)
    if (!sortmerna_idx.isDirectory())
        error "SortMeRNA index is not a directory: ${sortmerna_idx}"

    blastn_dir = file("${db}/blastn")
    if (!blastn_dir.isDirectory())
        error "Blast database directory is not a directory: ${blastn_dir}"

    ch_samples = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2) {
                error "samplesheet ${params.input}: every row needs sample, fastq_1, fastq_2 — got: ${row}"
            }
            tuple([id: row.sample],
                  [resolve(row.fastq_1, params.data_dir, "${row.sample}.fastq_1"),
                   resolve(row.fastq_2, params.data_dir, "${row.sample}.fastq_2")])
        }
        .toList()
        .map { rows ->
            if (!rows) error "samplesheet ${params.input} has no data rows"
            // Duplicates would collide in publishDir and overwrite each other.
            def dupes = rows.collect { it[0].id }.countBy { it }.findAll { it.value > 1 }.keySet()
            if (dupes) error "duplicate sample names in ${params.input}: ${dupes}"
            rows
        }
        .flatMap()

    ch_ref = Channel.value(sortmerna_ref)
    ch_idx = Channel.value(sortmerna_idx)

    FASTP(ch_samples)
    SORTMERNA(FASTP.out.reads, ch_ref, ch_idx)
    SPADES(SORTMERNA.out.clean)
    FILTER(SPADES.out.contigs)
    BLASTN(FILTER.out.contigs, Channel.value(blastn_dir), Channel.value(params.blastn_db))

    MULTIQC(
        FASTP.out.json.mix(SORTMERNA.out.log).collect()
    )
}
