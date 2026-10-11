# Added 2026-10-11 for the code release: collects the per-locus technical files written by aim1_coloc_full.R
# into the signal-mapping table and technical log used downstream (reconstructs the step originally run interactively).
# Usage: Rscript scripts/collect_aim1_tech.R out/aim1
suppressPackageStartupMessages(library(data.table))
OUT <- commandArgs(TRUE)[1]
fs <- list.files(file.path(OUT, "tech"), pattern = "\\.rds$", full.names = TRUE)
x <- lapply(fs, readRDS)
fwrite(rbindlist(lapply(x, `[[`, "signals"), fill = TRUE), file.path(OUT, "aim1_signal_mapping.csv"))
fwrite(rbindlist(lapply(x, `[[`, "tech"), fill = TRUE), file.path(OUT, "aim1_technical_log.csv"))
