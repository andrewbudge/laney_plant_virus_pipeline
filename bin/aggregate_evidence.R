#!/usr/bin/env Rscript

header <- c(
    "seq_name", "length", "assembly_cov", "genomad_virus_score",
    "genomad_fdr", "genomad_n_hallmarks", "genomad_marker_enrichment",
    "genomad_taxonomy", "map_numreads", "map_covbases", "map_coverage_pct",
    "map_meandepth", "viral_by_genomad"
)

args <- commandArgs(trailingOnly = TRUE)
fasta <- args[[1]]
virus_path <- args[[2]]
coverage_path <- args[[3]]

fasta_lines <- readLines(fasta)
fasta_headers <- fasta_lines[startsWith(fasta_lines, ">")]
seq_name <- sub("^>(\\S+).*", "\\1", fasta_headers)

length_match <- regexpr("length_([0-9]+)", seq_name, perl = TRUE)
length_values <- rep(NA_character_, length(seq_name))
length_values[length_match > 0] <- regmatches(seq_name, length_match)[length_match > 0]
length_values[length_match > 0] <- sub("^length_", "", length_values[length_match > 0])

cov_match <- regexpr("cov_([0-9.]+)", seq_name, perl = TRUE)
cov_values <- rep(NA_character_, length(seq_name))
cov_values[cov_match > 0] <- regmatches(seq_name, cov_match)[cov_match > 0]
cov_values[cov_match > 0] <- sub("^cov_", "", cov_values[cov_match > 0])

master <- data.frame(
    seq_name = seq_name,
    length = length_values,
    assembly_cov = cov_values,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

virus <- read.delim(
    virus_path, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE
)
virus <- virus[, c("seq_name", "virus_score", "fdr", "n_hallmarks", "marker_enrichment", "taxonomy"), drop = FALSE]
names(virus)[-1] <- paste0("genomad_", names(virus)[-1])
virus$viral_by_genomad <- 1

coverage <- read.delim(
    coverage_path, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, comment.char = ""
)
coverage <- coverage[, c("#rname", "numreads", "covbases", "coverage", "meandepth"), drop = FALSE]
names(coverage) <- c("seq_name", "map_numreads", "map_covbases", "map_coverage_pct", "map_meandepth")

merged <- merge(master, virus, by = "seq_name", all.x = TRUE, sort = FALSE)
merged <- merge(merged, coverage, by = "seq_name", all.x = TRUE, sort = FALSE)
merged <- merged[match(seq_name, merged$seq_name), , drop = FALSE]
merged$viral_by_genomad[is.na(merged$viral_by_genomad)] <- 0
merged <- merged[, header, drop = FALSE]

write.table(
    merged, file = "", sep = "\t", row.names = FALSE, col.names = TRUE,
    quote = FALSE, na = "NA"
)
