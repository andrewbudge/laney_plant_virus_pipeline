# plant-virus QC (Nextflow)

Adapter/quality trimming and rRNA depletion for plant virus RNA-seq, producing
paired clean reads ready for assembly plus one MultiQC report per run.

```
fastp ──► sortmerna ──► {sample}.clean_fwd.fq.gz + {sample}.clean_rev.fq.gz
  └──────────┴────────► multiqc_report.html
```

## Run it

```bash
./run_nf.sh
```

That is the whole thing. It loads the modules, points at scratch, and resumes
an interrupted run automatically.

To process different samples, edit `samplesheet.csv` — one row per sample,
nothing else to change:

```csv
sample,fastq_1,fastq_2
Laney18,raw/Laney18_1.fastq.gz,raw/Laney18_2.fastq.gz
```

## Where the data lives

Read paths in the samplesheet, and `--sortmerna_ref` / `--sortmerna_idx`, are
relative to `--data_dir`. That one flag is the only thing that changes when the
storage moves or when you run this somewhere else:

```bash
# CHPC
--data_dir /scratch/general/vast/abudge/plant_virus_pipeline_storage
# laptop
--data_dir ~/plant_virus_data
```

The samplesheet itself does not change between the two. `--input` is the
exception and stays a real path, because the sheet lives with the code rather
than with the data.

Reference databases get a second root, `--db_dir`, which defaults to
`--data_dir` when unset. They are split because they have a different lifecycle:
reads belong to one project and get archived with it, while references are
shared across projects and people, and are far too large to copy per project —
the sortmerna index alone is 1.8 GB, and a Kraken2 or DIAMOND database is tens
of GB. Once the references live somewhere shared:

```bash
--data_dir /scratch/general/vast/abudge/plant_virus_pipeline_storage \
--db_dir   /path/to/shared/databases \
--sortmerna_ref smr_v4.3_default_db.fasta \
--sortmerna_idx idx
```

Absolute paths in the samplesheet still work and ignore `--data_dir`, so an
existing sheet keeps running. Paths relative to your *current directory* are
deliberately not supported — a sheet that means something different depending
on where you `cd` is the one broken case that looks fine until it isn't.

Useful flags (pass straight through `run_nf.sh`):

| Flag | Effect |
|---|---|
| `-preview` | resolve the workflow, run nothing |
| `-resume` | reuse completed tasks (already on by default) |
| `--input other.csv` | different samplesheet |
| `--outdir DIR` | different output location |
| `--data_dir DIR` | root for relative read paths in the samplesheet |
| `--db_dir DIR` | root for relative reference/index paths (defaults to `--data_dir`) |

## Output

```
<outdir>/
  fastp/<sample>/       trimmed reads, .html, .json
  sortmerna/<sample>/   .clean_fwd/.clean_rev (use these), .rrna_*, .rrna.log
  multiqc/              multiqc_report.html
  pipeline_info/        timeline, report, trace, dag
```

`.clean_fwd.fq.gz` and `.clean_rev.fq.gz` are the reads to assemble, as an
ordinary paired-end library:

```bash
spades.py -1 <s>.clean_fwd.fq.gz -2 <s>.clean_rev.fq.gz ...
```

`pipeline_info/report.html` is what to send along when a run misbehaves — it
carries per-task exit codes, peak memory and runtime.

## Things that will bite you

**SortMeRNA is pinned to 4.3.4. Do not bump it.** Versions 4.3.6, 4.3.7 and the
7.0.0 release binary are compiled for x86-64-v3 (AVX2 + BMI2 + FMA). Lonepeak
nodes are Sandy Bridge E5-2680 — AVX only — so anything newer dies with
`Illegal instruction (core dumped)` after printing its version banner. 4.3.4 is
the last baseline-SSE build. Using a newer one means moving the SORTMERNA
process to notchpeak or granite; `conf/base.config` marks where that override
goes.

**SortMeRNA names its own outputs.** `--other X.clean` produces
`X.clean_fwd.fq.gz` and `X.clean_rev.fq.gz` — the argument is a *prefix*, and
the extension is appended by the tool. Renaming the `output:` declarations
produces a missing-output error, not a renamed file.

**The rRNA index is built once and reused.** `--sortmerna_idx` points at a
prebuilt index directory and skips the ~14 min / 4 GB `SORTMERNA_INDEX` step.
The index depends only on the reference FASTA, so it is reusable indefinitely;
drop the flag only if the reference changes.

**Editing a comment inside a `"""` script block costs a re-run.** Nextflow
hashes the whole script body, `#` comments included, so changing one
invalidates `-resume` for that process. Comments above `process`, or in the
`script:` section before the opening `"""`, are free. Worth knowing before
tidying a comment on a step that takes an hour and a half.

**Read merging was tried and parked.** Median insert here is 142 bp against
2x151 bp reads, so ~59% of pairs overlap and merging them is tempting. BBMerge
(`bbmap:39.81`) died instantly with `java.lang.AssertionError` at
`stream.WriterFactory.makeWriter` on this argument set:

```
bbmerge.sh -Xmx6g in1=R1 in2=R2 out=m.fq.gz outu1=u1.fq.gz outu2=u2.fq.gz ihist=h.txt
```

Worse than the crash: the main thread threw but BBTools' non-daemon worker
threads kept the JVM alive, so the job never exited and burned 80 minutes of
wall clock looking healthy to SLURM. The cause was most likely the argument
combination (`outu1=`/`outu2=` vs `outu=`/`outu2=`) rather than the tool, but
that was never isolated. If revisiting: `fastp --merge`, `vsearch
--fastq_mergepairs` and NGmerge are the alternatives; seqkit does not do
overlap merging at all. Merging is an optimization, not a repair -- fastp
already strips the adapter -- so assembling merged vs. unmerged and comparing
contig N50 is the way to decide whether it is worth any of this.

**Any process that runs a JVM needs a timeout mindset.** A Java tool whose main
thread dies can leave non-daemon threads holding the process open indefinitely.
SLURM sees RUNNING, Nextflow sees nothing, and the wall-clock limit is the only
thing that ends it -- at which point exit 143 lands in the retry list and it
happens twice more. Check `.command.log` and the task's own log rather than
trusting job state.

**Work directory grows.** Nextflow keeps every intermediate under
`$SCRATCH/nf_work` until told otherwise. After a run is confirmed good:

```bash
nextflow clean -f -before <run-name>    # run names: nextflow log
```

## Requirements

- `nextflow/25.04` and `apptainer/1.4.0` (both CHPC modules, loaded by `run_nf.sh`)
- A CHPC allocation; `conf/slurm.config` carries account/partition/qos
- Everything else ships in biocontainers, pulled automatically on first run
