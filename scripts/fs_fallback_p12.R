# Post hoc (2026-10-10, after independent audit): fallback coloc.abf used default p12 = 1e-5;
# recompute fallback PP.H4 under p12 = 5e-6 (same as coloc-SuSiE) by rescaling H4 prior odds.
suppressPackageStartupMessages(library(data.table))
F <- fread("out/aim1/f_s_signals.csv")
F[, f_s_p12 := f_s]
for (i in which(F$method == "fallback_abf" & !is.na(F$f_s))) {
  z <- readRDS(file.path("out/aim1/sealed", paste0(F$locus[i], ".rds")))
  s <- z$fallback[[F$signal_id[i]]]$summary
  pp <- as.numeric(s[c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")])
  new4 <- 0.5 * pp[5] / (sum(pp[1:4]) + 0.5 * pp[5])
  if (isTRUE(F$f_s[i] > 0)) F$f_s_p12[i] <- new4
}
fwrite(F, "out/aim1/f_s_signals_p12.csv")
c <- F[!is.na(f_s)]
cat("fallback ge0.8 old/new:", sum(c$method=="fallback_abf" & c$f_s>=0.8), sum(c$method=="fallback_abf" & c$f_s_p12>=0.8), "\n")
cat("all ge0.8 / lt0.2 new:", sum(c$f_s_p12>=0.8), sum(c$f_s_p12<0.2), " median", median(c$f_s_p12), "\n")
