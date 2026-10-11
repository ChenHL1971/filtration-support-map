# Unseal (approved 2026-10-10) and compute f_s per S8 signal per locked rules
suppressPackageStartupMessages(library(data.table))
S <- fread("out/aim1/aim1_signal_mapping.csv")
res <- list()
for (loc in unique(S$locus)) {
  z <- readRDS(file.path("out/aim1/sealed", paste0(loc, ".rds")))
  for (i in which(S$locus == loc)) {
    r <- S[i]; f <- NA_real_; pp4 <- NA_real_; dir <- NA; note <- ""
    if (r$method %in% c("overlap", "LD")) {
      lead <- z$lead_z$lead[z$lead_z$cs == r$susie_cs]
      sm <- z$susie[hit1 == lead]
      if (length(lead) && nrow(sm)) {
        j <- which.max(sm$PP.H4.abf); pp4 <- sm$PP.H4.abf[j]
        lz <- z$lead_z[cs == r$susie_cs]; dir <- sign(lz$z_crea) == sign(lz$z_cys)
      } else note <- "susie pair trimmed"
    } else if (r$method == "fallback_abf") {
      fb <- z$fallback[[r$signal_id]]; pp4 <- fb$summary[["PP.H4.abf"]]; dir <- sign(fb$z_crea) == sign(fb$z_cys)
    }
    if (!is.na(pp4)) f <- if (isTRUE(dir)) pp4 else if (isFALSE(dir)) 0 else NA_real_
    if (!is.na(f) && r$cs_coverage < 0.8) { f <- NA_real_; note <- "coverage < 0.8" }
    res[[length(res) + 1]] <- data.table(signal_id = r$signal_id, locus = loc, method = r$method, cs_size = r$cs_size,
                                         cs_coverage = r$cs_coverage, pp4 = pp4, same_dir = dir, f_s = f, note = note)
  }
}
F <- rbindlist(res); fwrite(F, "out/aim1/f_s_signals.csv")
cat("signals", nrow(F), " classifiable", sum(!is.na(F$f_s)), "\n"); print(F[, .(n = .N, classifiable = sum(!is.na(f_s)), median_fs = round(median(f_s, na.rm = TRUE), 3),
  fs_ge_0.8 = sum(f_s >= 0.8, na.rm = TRUE), opposite_dir = sum(same_dir == FALSE, na.rm = TRUE)), by = method])
print(table(F$note)); cat("overall f_s quantiles:", round(quantile(F$f_s, c(.1,.25,.5,.75,.9), na.rm=TRUE),3), "\n")
cat("f_s>=0.8:", sum(F$f_s >= 0.8, na.rm=TRUE), " f_s<0.2:", sum(F$f_s < 0.2, na.rm=TRUE), "\n")
