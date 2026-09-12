# plant-virus QC + assembly (Nextflow)

Adapter/quality trimming and rRNA depletion for plant virus RNA-seq, then de
novo RNA viral assembly (SPAdes `--rnaviral`), plus one MultiQC report per run.

```
fastp ──► sortmerna ──► spades --rnaviral ──► scaffolds.fasta
  └──────────┴────────────► multiqc_report.html
```

This is the QC + assembly stage of a viral-detection pipeline; viral screening
comes next.

## Layout

```
main.nf                workflow wiring + samplesheet/path validation
conf/base.config       per-process resources + retry policy
conf/modules.config    per-process containers, publishDir, tool args
conf/slurm.config      CHPC executor (profile slurm)
modules/local/         one process per file: fastp, sortmerna, spades, multiqc
bin/tagfind            biocontainer tag lookup helper
run.sh / run_local.sh  HPC and local launchers
```

## Run it

FASTQ paths in the TSV resolve against `--data_dir`. Every reference and
formatted search database lives under one `--db` root. `--db` may be absolute
or relative to `--data_dir`.

```
<data_dir>/
  raw/                            input reads
  databases/                      pass as --db databases
    sortmerna/
      reference.fasta
      index/                      prebuilt SortMeRNA index
    blastn/                       future viral-screening databases
```

The TSV is tab-separated and has exactly these columns:

```text
sample	fastq_1	fastq_2
Laney18	raw/Laney18_1.fastq.gz	raw/Laney18_2.fastq.gz
```

### On CHPC SLURM (production)

```bash
module load nextflow/25.04 apptainer/1.4.0
PVP=/scratch/general/vast/abudge/plant_virus_pipeline_storage
export NXF_APPTAINER_CACHEDIR="$PVP/apptainer_nf"

nextflow run andrewbudge/laney_plant_virus_pipeline \
  -r main \
  -profile slurm \
  -work-dir "$PVP/nf_work" \
  --data_dir "$PVP" \
  --input "$HOME/samplesheet.tsv" \
  --outdir "$PVP/qc" \
  --db databases
```

This is a fresh run; it deliberately omits `-resume`. Pin `-r` to a release tag
instead of `main` when reproducibility matters. Public GitHub repositories need
no credentials. This repository is currently private, so export a read-only
GitHub PAT as `GITHUB_TOKEN` on the login node (or configure
`$HOME/.nextflow/scm`); never put the token in this repository or a run script.

### Locally (development)

The local box has apptainer, same as CHPC, so behavior matches. Copy a small
slice of reads plus the `databases/` dir into a local instance dir, then:

```bash
PVP_STORAGE=~/plant_virus_data ./run_local.sh
```

To skip image pulls, copy the HPC `apptainer_nf/` cache over to the same dir.
For speed during iteration, subset the reads first (e.g. `zcat
raw/Laney18_1.fastq.gz | head -n 800000 | gzip > raw/sub_1.fastq.gz`) and put
the subset paths in `samplesheet.tsv`.

## Output

```
<outdir>/
  fastp/<sample>/       trimmed reads, .html, .json
  sortmerna/<sample>/   .clean_fwd/.clean_rev (assembly input), .rrna_*, .rrna.log
  spades/<sample>/      contigs.fasta, scaffolds.fasta, spades.log
  multiqc/              multiqc_report.html
  pipeline_info/        timeline, report, trace, dag
```

`spades/<sample>/scaffolds.fasta` is the assembled viral contigs (contigs.fasta
is the unfiltered set).

`pipeline_info/report.html` is what to send along when a run misbehaves — it
carries per-task exit codes, peak memory and runtime.

## Things that will bite you

**SortMeRNA is pinned to 4.3.4. Do not bump it.** 4.3.6, 4.3.7 and the 7.0.0
binary are built for x86-64-v3 (AVX2+BMI2+FMA) and die with `Illegal
instruction (core dumped)` on lonepeak's Sandy Bridge nodes. 4.3.4 is the last
baseline-SSE build.

**SortMeRNA names its own outputs.** `--other X.clean` produces
`X.clean_fwd.fq.gz` and `X.clean_rev.fq.gz` — the argument is a *prefix*, and
the extension is appended by the tool. Renaming the `output:` declarations
produces a missing-output error, not a renamed file.

**The SortMeRNA index must match its reference.** Keep the prebuilt index at
`<db>/sortmerna/index/` beside `<db>/sortmerna/reference.fasta`; the pipeline
fails before submitting work if either is absent.

## Requirements

- CHPC: `nextflow/25.04` and `apptainer/1.4.0` (loaded by `run.sh`), plus a
  CHPC allocation (`conf/slurm.config` carries account/partition/qos)
- Local: nextflow 23.10+, apptainer; the storage dir layout above
- Everything else ships in biocontainers, pulled automatically on first run
