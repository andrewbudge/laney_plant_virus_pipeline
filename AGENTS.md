# AGENTS.md — laney_plant_virus_pipeline (project notes for the coding agent)

Plant virus DETECTION pipeline for RNA-seq, Nextflow (nf-core-lite layout).
Owner: Andrew Budge. README.md is .gitignored until release; treat this file
as the live design record.

## Pipeline (current)

```
fastp ─► sortmerna ─► spades --rnaviral ─► filter (len/cov)
  └──────────┴──────────► multiqc_report.html
         filter/contigs ─► diamond RVDB-prot · diamond UniRef90 · geNomad
                            plus bowtie2 read mapping
```

- `modules/local/` — one process per file: fastp, sortmerna, sortmerna_index
  (removed; index now prebuilt), spades, filter, multiqc. blastn is deferred.
- Conf layout: `conf/base.config` (resources/retry), `conf/modules.config`
  (container/publishDir/ext.args), `conf/slurm.config` (CHPC executor).
- `main.nf` validates `--db` layout up front: sortmerna ref + index +
  genomad/genomad_db must exist before any process is submitted.

## Interface

```
nextflow run andrewbudge/laney_plant_virus_pipeline -profile slurm|local \
  --data_dir <root> --input samplesheet.tsv --outdir <dir> --db <dir>
```

- TSV columns: `sample  fastq_1  fastq_2`. Relative FASTQ paths resolve against
  `--data_dir`; absolute paths pass through. Relative `--db` also resolves
  against `--data_dir`.
- No `-resume` in the launchers (deliberate, documented). Pin `-r` to a tag for
  reproducibility.
- Repo is PUBLIC on GitHub; `nextflow run owner/repo` needs no token. `-latest`
  fetches latest main; run with `-r main` normally.

## DB layout (contract)

```
<db>/
  sortmerna/smr_v4.3_default_db.fasta   # MUST keep this exact basename
  sortmerna/index/                      # prebuilt index; name-encodes the ref basename
  blastn/U-RVDBv32.0.*                  # deferred/optional; makeblastdb -dbtype nucl -parse_seqids
  blastx/                               # RVDB-prot.fasta, uniref90.fasta.gz (raw)
  blastx/uniref90.dmnd                  # Diamond UniRef90 database
  genomad/genomad_db/                   # genomad download-database output
```

Index must match reference basename (Nextflow stages symlinked refs under the
link name → SortMeRNA index mismatch). Do NOT symlink the reference.

## Pinned images (AVX constraints — lonepeak is Sandy Bridge, AVX only)

| tool | image | why |
|---|---|---|
| fastp | fastp:0.23.4--h125f33a_5 | |
| sortmerna | sortmerna:4.3.4--h9ee0642_0 | 4.3.6+/7.0.0 need x86-64-v3 (AVX2+BMI2+FMA) → SIGILL |
| spades | spades:4.0.0--haf24da9_4 | 4.2+ needs AVX2 → SIGILL; 4.0.0 supports --rnaviral |
| seqkit | seqkit:2.13.0--he881be0_0 | Go binary, AVX-safe |
| blast | blast:2.17.0--hb02a186_1 | blastn deferred; 2.17.0 verified working on lonepeak |
| bowtie2 | bowtie2:2.5.5--ha27dd3b_0 | SSE-era code, AVX-safe |
| genomad | genomad:1.12.0--pyhdfd78af_0 | SIGILL risk: bundled mmseqs2/TF need modern CPUs; verify before production, may need non-lonepeak partition |
| multiqc | multiqc:1.25.2--pyhdfd78af_0 | |

Test any new container with `apptainer exec <img> <tool> --version` before
trusting it on lonepeak. tagfind (`./bin/tagfind <tool> [prefix]`) lists
biocontainer tags. diamond:2.2.6 flagged AVX2-risk — verify before use.

## Screening strategy (V1 LOCKED with user)

Goal: **curate evidence for experts to sift** — no prefiltration, accrue info.
Host-agnostic; can't guarantee completeness. Three active parallel legs on
FILTER's contigs, each feeds one raw per-leg TSV; blastn is deferred.

| leg | tool → DB | what it answers |
|---|---|---|
| 1 | blastn → U-RVDB (nucleotide) | (deferred) known-virus relative? (homology) |
| 2 | diamond blastx → RVDB-prot | divergent-virus homology (protein) (implemented) |
| 3 | diamond blastx → UniRef90 | unbiased: is the best hit viral? (confirmation, implemented) |
| 4 | geNomad | viral by sequence/marker signal, no homology needed (implemented, confirmed) |
| support | bowtie2 map-back reads → contigs | per-contig read support and confidence for calls |

Decisions made:
- **UniRef90 only — no nr.** UniRef90 ⊂ nr for any known virus; nr is ~500M seqs
  (hours to format, huge disk) and adds nothing except for viruses absent from
  UniRef90. Confirmation leg = contigs themselves vs UniRef90, NOT the viral
  hit sequences.
- **geNomad only for V1 — no VirSorter2.** geNomad was retrained for RNA +
  giant viruses (99.9% of RdRP contigs flagged vs ~44% for others) and beats
  all tools in its benchmark (MCC 95.3%); handles <3kb contigs; ~26x faster
  than VirSorter2. Skip k-mer CNN tools (PPR-Meta/DeepVirFinder/Seeker) — they
  false-positive on eukaryotic seqs at 2-3x rate. VirSorter2 (`--include-groups
   RNA`) is the agreed backup if geNomad underperforms.
- **blastn is deferred while geNomad + bowtie2 map-back are the focus.**
- Evidence aggregation deferred until diamond + geNomad outputs compared; schema
  to be developed from real data.
- **Map-back support is a confidence leg.** rRNA-depleted reads mapped back to
  filtered contigs support calls such as the RNA mitovirus; BAM enables
  `samtools depth`/`idxstats` for the deferred evidence table.

DBs needed (sources already staged in `<db>/blastx/` + user downloading):
- U-RVDB v32.0 (rvdb.dbi.udel.edu), deferred → `makeblastdb -dbtype nucl -parse_seqids`
- RVDB-prot (Institut Pasteur, rvdb-prot.pasteur.fr) → `diamond makedb`
- UniRef90 (uniref90.fasta.gz on HPC) → `diamond makedb`
- geNomad DB → `genomad download-database` on HPC

Note: RVDB includes endogenous retrovirus / LTR-retrotransposon sequences —
plant hosts will hit them. That's expected; UniRef90 best-hit sorts it.
RVDB "Putative Non Viral Annotation" file exists if we want to tag hits later.

Notes: blastx >> blastn for sensitivity on divergent viruses (protein diverges
~10x slower). "Both legs viral" = confidence flag, NOT a hard filter. Diamond =
100-1000x faster than blastx. Diamond 2.2.6 AVX2-risk on lonepeak — verify
`diamond --version` before use.

## Output

```
<outdir>/
  fastp/<sample>/  sortmerna/<sample>/  spades/<sample>/
  filter/<sample>/<sample>.contigs.fasta   # len>=contig_min_length(1000), cov>=contig_min_cov(10)
  diamond/<sample>/{<sample>.rvdb.tsv,<sample>.uniref90.tsv}
  genomad/<sample>/{virus_summary.tsv,virus.fna,virus_proteins.faa,summary.json,*.genomad.log}
  bowtie2/<sample>/{*.sorted.bam,*.sorted.bam.bai,*.coverage.tsv,*.bowtie2.log}
  multiqc/  pipeline_info/{timeline,report,trace,dag}
```

blastn remains deferred; its planned outfmt is deliberately 12 columns.
qcovs = fraction of contig matched, `length` vs `slen` = does the contig span
the whole genome. staxids/sscinames stay blank without `makeblastdb -taxid_map`.

## Verification workflow

Full local runtime tests with synthetic data under /tmp (proven pattern):
random 3kb genome → sortmerna ref + index; second random 3kb genome → 500 FR
151bp reads (p in [0,L-302], r1=B[p:p+151], r2=rc(B[p+151:p+302])) + assembly
target. Validate publishDir closures + output declarations only at runtime
(preview doesn't catch them). Always: `nextflow run -preview`, `git diff --check`.

## HPC

- lonepeak SLURM, account ogden. Images cached in $PVP/apptainer_nf via
  NXF_APPTAINER_CACHEDIR. PVP=/scratch/general/vast/abudge/plant_virus_pipeline_storage.
- Work dir $PVP/nf_work, outdir $PVP/qc_assembly. Production run:
  `nextflow run andrewbudge/laney_plant_virus_pipeline -latest -r main -resume
  -profile slurm -work-dir "$PVP/nf_work" --data_dir "$PVP"
  --input "$HOME/samplesheet.tsv" --outdir "$PVP/qc_assembly" --db databases`
- Resources measured on Laney18 (18.2M pairs): FASTP 2m/3GB, SORTMERNA 1h24m/4.5GB.
  SPADES seeds: 8 cpus/32GB/12h. Log real runtimes in base.config comments.
