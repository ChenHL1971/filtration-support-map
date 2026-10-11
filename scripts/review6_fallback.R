# Post hoc (review round 6): fallback coloc.abf calibration checks.
# (1) class distribution by method and credible-set size; (2) p12 varied for fallback signals only, direct signals fixed.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
F <- fread("out/aim1/f_s_signals.csv"); S <- fread("out/aim1/aim1_signal_mapping.csv")
# cs_size already in f_s_signals.csv
F <- F[!is.na(f_s)]
F[, grp := ifelse(method == "fallback_abf", "fallback", "direct")]
F[, cls := fifelse(f_s >= 0.8, "consistent", fifelse(f_s >= 0.2, "intermediate", "low"))]
F[, size_bin := cut(cs_size, c(0, 1, 5, 20, Inf), labels = c("1", "2-5", "6-20", ">20"))]
by_method <- dcast(F[, .N, by = .(grp, cls)], grp ~ cls, value.var = "N", fill = 0)
by_size <- dcast(F[, .N, by = .(grp, size_bin, cls)], grp + size_bin ~ cls, value.var = "N", fill = 0)
size_q <- F[, .(median = as.numeric(median(cs_size)), q1 = as.numeric(quantile(cs_size, .25)), q3 = as.numeric(quantile(cs_size, .75))), by = grp]
resc <- function(pp, k) { pp <- as.numeric(pp); k * pp[5] / (sum(pp[1:4]) + k * pp[5]) }
grid <- c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
fb <- F[grp == "fallback"]
vals <- t(sapply(seq_len(nrow(fb)), function(i) {
  r <- fb[i]; z <- readRDS(file.path("out/aim1/sealed", paste0(r$locus, ".rds")))
  pp <- z$fallback[[r$signal_id]]$summary[c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")]
  v <- sapply(grid, function(g) resc(pp, g / 1e-5)); if (isFALSE(r$same_dir)) v[] <- 0; v }))
ndir <- sum(F[grp == "direct"]$f_s >= 0.8)
fb_only <- data.table(p12 = grid, fallback_consistent = colSums(vals >= 0.8), direct_consistent = ndir,
                      total_consistent = colSums(vals >= 0.8) + ndir)
res <- list(n = F[, .N, by = grp], cs_size = size_q, by_method = by_method, by_size = by_size, fallback_only_p12 = fb_only)
write_json(res, "out/review6_fallback.json", auto_unbox = TRUE, pretty = TRUE, digits = NA); print(res)
