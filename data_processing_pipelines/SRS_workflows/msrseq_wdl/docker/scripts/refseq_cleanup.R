#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: refseq_cleanup.R <in_fa> <out_fa>\n")
  quit(status = 1)
}

in_fa  <- args[1]
out_fa <- args[2]

lines <- readLines(in_fa, warn = FALSE)
is_header <- startsWith(lines, ">")
lines[is_header] <- paste0(">", sub("=.*$", "", substring(lines[is_header], 2)))
writeLines(lines, out_fa)
message("Wrote cleaned FASTA to: ", normalizePath(out_fa, winslash = "/", mustWork = FALSE))
