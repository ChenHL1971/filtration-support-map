# Post hoc (review round 5): PP.H4 of direction-rule signals; class counts under alternative p12 (1e-6 to 1e-4).
# Rescales H4 prior odds by k = p12_new / p12_used (coloc-SuSiE used 5e-6; fallback coloc.abf used 1e-5).
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
F <- fread("out/aim1/f_s_signals.csv"); S <- fread("out/aim1/aim1_signal_mapping.csv")
F <- merge(F, S[, .(signal_id, susie_cs)], by = "signal_id", all.x = TRUE)
resc <- function(pp, k) { pp <- as.numeric(pp); k * pp[5] / (sum(pp[1:4]) + k * pp[5]) }
grid <- c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
rows <- list()
for (i in seq_len(nrow(F))) {
  r <- F[i]; if (is.na(r$f_s)) next
  z <- readRDS(file.path("out/aim1/sealed", paste0(r$locus, ".rds")))
  if (r$method == "fallback_abf") {
    pp <- z$fallback[[r$signal_id]]$summary[c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")]; used <- 1e-5
    vals <- sapply(grid, function(g) resc(pp, g / used))
  } else {
    lead <- z$lead_z$lead[z$lead_z$cs == r$susie_cs]; sm <- as.data.frame(z$susie)
    sm <- sm[sm$hit1 == lead, , drop = FALSE]; used <- 5e-6
    vals <- sapply(grid, function(g) max(apply(sm[, c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")], 1, resc, k = g / used)))
  }
  if (isFALSE(r$same_dir)) vals[] <- 0
  rows[[length(rows) + 1]] <- data.table(signal_id = r$signal_id, method = r$method, t(setNames(vals, paste0("p12_", grid))))
}
R <- rbindlist(rows)
cnt <- R[, lapply(.SD, function(v) c(consistent = sum(v >= 0.8), intermediate = sum(v >= 0.2 & v < 0.8), low = sum(v < 0.2))), .SDcols = patterns("p12_")]
cnt[, class := c("consistent", "intermediate", "low")]
opp <- F[!is.na(f_s) & same_dir == FALSE]
q <- quantile(opp$pp4, c(0, .25, .5, .75, 1), na.rm = TRUE)
res <- list(direction_rule = list(n = nrow(opp), pp4_quantiles = round(q, 4), n_pp4_ge_0.8 = sum(opp$pp4 >= 0.8, na.rm = TRUE),
                                  n_pp4_ge_0.5 = sum(opp$pp4 >= 0.5, na.rm = TRUE), n_pp4_lt_0.2 = sum(opp$pp4 < 0.2, na.rm = TRUE)),
            p12_class_counts = cnt, check_original_counts = c(consistent = sum(F$f_s >= 0.8, na.rm = TRUE), low = sum(F$f_s < 0.2, na.rm = TRUE)))
write_json(res, "out/review5_coloc_checks.json", auto_unbox = TRUE, pretty = TRUE, digits = NA); print(toJSON(res, auto_unbox = TRUE, pretty = TRUE))
