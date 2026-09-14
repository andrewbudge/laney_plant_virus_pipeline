#!/usr/bin/env python3
import csv
import re
import sys


HEADER = [
    "seq_name", "length", "assembly_cov", "genomad_virus_score",
    "genomad_fdr", "genomad_n_hallmarks", "genomad_marker_enrichment",
    "genomad_taxonomy", "map_numreads", "map_covbases", "map_coverage_pct",
    "map_meandepth", "viral_by_genomad",
]
VIRUS_COLUMNS = ["virus_score", "fdr", "n_hallmarks", "marker_enrichment", "taxonomy"]
COVERAGE_COLUMNS = ["numreads", "covbases", "coverage", "meandepth"]


def rows(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def contigs(path):
    result = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                seq_name = line[1:].split(None, 1)[0]
                length = re.search(r"length_(\d+)", seq_name)
                assembly_cov = re.search(r"cov_([0-9.]+)", seq_name)
                result.append((seq_name, length.group(1) if length else "NA",
                               assembly_cov.group(1) if assembly_cov else "NA"))
    return result


def load(path, key, columns):
    result = {}
    for row in rows(path):
        result[row[key]] = {column: row.get(column, "") for column in columns}
    return result


def main():
    fasta, virus_path, coverage_path = sys.argv[1:]
    virus = load(virus_path, "seq_name", VIRUS_COLUMNS)
    coverage = load(coverage_path, "#rname", COVERAGE_COLUMNS)
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(HEADER)
    for seq_name, length, assembly_cov in contigs(fasta):
        virus_row = virus.get(seq_name, {})
        coverage_row = coverage.get(seq_name, {})
        writer.writerow([
            seq_name, length, assembly_cov,
            *(virus_row.get(column, "NA") or "NA" for column in VIRUS_COLUMNS),
            *(coverage_row.get(column, "NA") or "NA" for column in COVERAGE_COLUMNS),
            1 if seq_name in virus else 0,
        ])


if __name__ == "__main__":
    main()
