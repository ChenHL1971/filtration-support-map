# Map eGFRcrea SuSiE credible sets to Liu S8 signals (rule drafted 2026-10-10; no PP.H4 used)
suppressPackageStartupMessages({library(data.table); library(readxl)})
args <- commandArgs(TRUE); SUPP <- args[1]; SUS <- args[2]; REF <- args[3]; OUT <- args[4]
src <- readLines("scripts/coloc_pilot.R"); a <- grep("^read_bed <- function", src); b <- grep("read_bed\\(REF\\)", src)
eval(parse(text = src[a:(b - 1)])); ref <- read_bed(REF)
s8 <- as.data.table(read_excel(SUPP, sheet = "S8 Finemap.Multi", skip = 1))
s8 <- s8[!is.na(RSID) & grepl("^chr", `Merged GWAS loci`)]
s8[, signal_id := paste0(`Merged GWAS loci`, "|", `CARMA Credible Set`)]
s8[, pip := as.numeric(`CARMA posterior inclusion probability (PIP)`)]
S <- readRDS(SUS); rows <- list()
r2 <- function(v1, v2) { i <- match(c(v1, v2), ref$bim$rsid); if (anyNA(i)) return(NA_real_)
  x <- ref$G[, i]; x <- x[complete.cases(x), , drop = FALSE]; if (nrow(x) < 50 || any(apply(x, 2, sd) == 0)) return(NA_real_); cor(x)[1, 2]^2 }
for (loc in names(S)) {
  X <- S[[loc]]; sig <- s8[`Merged GWAS loci` == loc]
  for (sid in unique(sig$signal_id)) {
    v <- sig[signal_id == sid]; top <- v$RSID[which.max(v$pip)]
    if (is.null(X$alpha) || nrow(X$alpha) == 0) { rows[[length(rows) + 1]] <- data.table(locus = loc, signal_id = sid, match = "no eGFRcrea credible set"); next }
    ov <- apply(X$alpha, 1, function(al) sum(al[X$snps %in% v$RSID]))
    k <- which.max(ov)
    if (ov[k] >= 0.5) { rows[[length(rows) + 1]] <- data.table(locus = loc, signal_id = sid, susie_cs = k, overlap = round(ov[k], 3), match = "PIP overlap >= 0.5"); next }
    leads <- X$snps[apply(X$alpha, 1, which.max)]
    rr <- sapply(leads, function(l) r2(top, l)); kk <- which.max(rr)
    if (length(kk) && !is.na(rr[kk]) && rr[kk] >= 0.8) rows[[length(rows) + 1]] <- data.table(locus = loc, signal_id = sid, susie_cs = kk, overlap = round(ov[kk], 3), r2_lead = round(rr[kk], 3), match = "LD r2 >= 0.8")
    else rows[[length(rows) + 1]] <- data.table(locus = loc, signal_id = sid, overlap = round(ov[k], 3), r2_lead = if (length(kk)) round(rr[kk], 3) else NA, match = "unmatched")
  }
}
M <- rbindlist(rows, fill = TRUE)
M[match != "unmatched" & !is.na(susie_cs), shared := .N > 1, by = .(locus, susie_cs)]
fwrite(M, OUT); print(M[, .N, by = match]); print(M[shared == TRUE, .N]); print(M)
