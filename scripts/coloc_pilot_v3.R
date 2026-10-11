# Aim 1 coloc-SuSiE pilot (technical feasibility only; Plan-5 locked).
# Checks input validity, LD/allele-order compatibility, numerical stability, runtime, memory.
# PP.H4 values are written to a sealed file and NOT printed or summarized.
suppressPackageStartupMessages({library(data.table); library(susieR); library(coloc)})
set.seed(20261008)
args <- commandArgs(TRUE); GDIR <- args[1]; REF <- args[2]; LOCI <- args[3]; OUT <- args[4]

## ---- PLINK .bed reader (SNP-major); dosage = count of bim allele 2 ----
read_bed <- function(prefix) {
  bim <- fread(paste0(prefix, ".bim"), header = FALSE, col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
  fam <- fread(paste0(prefix, ".fam"), header = FALSE)
  n <- nrow(fam); nb <- ceiling(n / 4)
  con <- file(paste0(prefix, ".bed"), "rb"); raw <- readBin(con, "raw", 3 + nb * nrow(bim)); close(con)
  stopifnot(identical(raw[1:3], as.raw(c(0x6c, 0x1b, 0x01))))
  lut <- c(0, NA, 1, 2)                       # 00 hom a1, 01 missing, 10 het, 11 hom a2
  G <- matrix(NA_real_, n, nrow(bim))
  for (j in seq_len(nrow(bim))) {
    b <- as.integer(raw[3 + (j - 1) * nb + seq_len(nb)])
    codes <- as.vector(rbind(b %% 4, (b %/% 4) %% 4, (b %/% 16) %% 4, (b %/% 64) %% 4))[seq_len(n)]
    G[, j] <- lut[codes + 1]
  }
  list(bim = bim, G = G)
}
message(Sys.time(), " reading bed"); ref <- read_bed(REF); message(Sys.time(), " bed done")
comp <- c(A = "T", T = "A", C = "G", G = "C")

## ---- GWAS inputs ----
liu <- fread(paste0(GDIR, "/liu_crea.gz"), select = c("MarkerName", "CHR", "POS", "REF", "ALT", "MAF", "Zscore", "N"))
setnames(liu, c("rsid", "chr", "pos", "nea", "ea", "maf", "z", "n"))
cys <- fread(paste0(GDIR, "/cys.gz"), select = c("RSID", "Allele1", "Allele2", "n", "Freq1", "Effect", "P.value", "Direction"))
cys[, sgn := fifelse(Effect > 0, 1, fifelse(Effect < 0, -1,
             fifelse(grepl("^\\++$", Direction), 1, fifelse(grepl("^-+$", Direction), -1, NA_real_))))]
cys[, z := sgn * qnorm(P.value / 2, lower.tail = FALSE)]
cys <- cys[is.finite(z), .(rsid = RSID, ea = toupper(Allele1), nea = toupper(Allele2), maf = pmin(Freq1, 1 - Freq1), z, n)]
liu[, `:=`(ea = toupper(ea), nea = toupper(nea))]

align <- function(d, bim) {        # returns z oriented to bim a2 (dosage allele); drops mismatches and ambiguous palindromes
  m <- merge(d, bim[, .(rsid, a1, a2)], by = "rsid")
  pal <- (paste0(m$a1, m$a2) %in% c("AT", "TA", "CG", "GC")) & m$maf >= 0.42
  same <- m$ea == m$a2 & m$nea == m$a1; swap <- m$ea == m$a1 & m$nea == m$a2
  ispal <- paste0(m$a1, m$a2) %in% c("AT", "TA", "CG", "GC")
  fl <- !ispal & comp[m$ea] == m$a2 & comp[m$nea] == m$a1; fls <- !ispal & comp[m$ea] == m$a1 & comp[m$nea] == m$a2
  fl[is.na(fl)] <- FALSE; fls[is.na(fls)] <- FALSE
  m[, z := fifelse(same, z, fifelse(swap, -z, fifelse(fl, z, fifelse(fls, -z, NA_real_))))]
  m[!pal & is.finite(z)]
}
std_beta <- function(z, maf, n) { v <- 1 / (2 * maf * (1 - maf) * (n + z^2)); list(beta = z * sqrt(v), varbeta = v) }

loci <- fread(LOCI, colClasses = list(character = "chr"))
log <- list(); sealed <- list(); susie_keep <- list()
for (i in seq_len(nrow(loci))) {
  message(Sys.time(), " start ", loci$locus[i]); L <- loci[i]; t0 <- proc.time()[["elapsed"]]; gc(reset = TRUE)
  rec <- list(locus = L$locus, n_signals_S8 = L$n_signals, status = "ok", note = "")
  res <- tryCatch({
    bim <- ref$bim[chr == L$chr & pos >= L$start & pos <= L$end]
    a <- align(liu[chr == L$chr & pos >= L$start & pos <= L$end], bim)
    b <- align(cys[rsid %in% bim$rsid], bim)
    snps <- Reduce(intersect, list(a$rsid, b$rsid, bim$rsid))
    idx <- match(snps, ref$bim$rsid); X <- ref$G[, idx, drop = FALSE]
    keep <- apply(X, 2, function(x) sd(x, na.rm = TRUE) > 0 && mean(is.na(x)) < 0.05)
    snps <- snps[keep]; X <- X[, keep, drop = FALSE]
    X <- apply(X, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
    R <- cor(X); dimnames(R) <- list(snps, snps)
    a <- a[match(snps, rsid)]; b <- b[match(snps, rsid)]
    rec$variants_liu_in_region <- nrow(liu[chr == L$chr & pos >= L$start & pos <= L$end])
    rec$variants_used <- length(snps)
    rec$ld_symmetric <- isSymmetric(R)
    ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    rec$ld_min_eigen <- signif(min(ev), 3)
    rec$ld_mismatch_s_liu <- signif(estimate_s_rss(a$z, R, n = median(a$n)), 3)
    rec$ld_mismatch_s_cys <- signif(estimate_s_rss(b$z, R, n = median(b$n)), 3)
    ba <- std_beta(a$z, a$maf, a$n); bb <- std_beta(b$z, b$maf, b$n)
    D1 <- list(beta = setNames(ba$beta, snps), varbeta = setNames(ba$varbeta, snps), N = median(a$n), sdY = 1,
               type = "quant", snp = snps, LD = R)
    D2 <- list(beta = setNames(bb$beta, snps), varbeta = setNames(bb$varbeta, snps), N = median(b$n), sdY = 1,
               type = "quant", snp = snps, LD = R)
    check_dataset(D1, req = "LD"); check_dataset(D2, req = "LD")
    w <- character(0)
    S1 <- withCallingHandlers(runsusie(D1, repeat_until_convergence = FALSE, maxit = 1000), warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
    S2 <- withCallingHandlers(runsusie(D2, repeat_until_convergence = FALSE, maxit = 1000), warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
    rec$susie_converged_liu <- isTRUE(S1$converged); rec$susie_converged_cys <- isTRUE(S2$converged)
    rec$credible_sets_liu <- length(S1$sets$cs); rec$credible_sets_cys <- length(S2$sets$cs)
    rec$warnings <- paste(unique(substr(w, 1, 120)), collapse = " | ")
    susie_keep[[L$locus]] <- list(alpha = S1$alpha[S1$sets$cs_index, , drop = FALSE], cs = S1$sets$cs, snps = snps, z = setNames(a$z, snps))
    if (rec$credible_sets_liu > 0 && rec$credible_sets_cys > 0) {
      C <- coloc.susie(S1, S2); rec$coloc_method <- "coloc.susie"
      rec$coloc_pairs_tested <- nrow(C$summary); sealed[[L$locus]] <- C$summary
    } else if (rec$credible_sets_liu > 0) {
      # rule (2026-10-10): eGFRcys has no credible set -> each eGFRcrea SuSiE signal vs eGFRcys single-signal ABF
      bf1 <- S1$lbf_variable[S1$sets$cs_index, , drop = FALSE]; colnames(bf1) <- snps
      W2 <- 0.15^2; r <- W2 / (W2 + D2$varbeta); bf2 <- setNames(0.5 * (log(1 - r) + r * (D2$beta^2 / D2$varbeta)), snps)
      C <- coloc.bf_bf(bf1, bf2); rec$coloc_method <- "coloc.bf_bf (eGFRcrea SuSiE x eGFRcys ABF)"
      rec$coloc_pairs_tested <- nrow(C$summary); sealed[[L$locus]] <- C$summary
    } else {
      rec$coloc_method <- "none"; rec$coloc_pairs_tested <- 0; rec$note <- "no eGFRcrea credible set"
    }
    TRUE
  }, error = function(e) { rec$status <<- "error"; rec$note <<- conditionMessage(e); FALSE })
  rec$runtime_sec <- round(proc.time()[["elapsed"]] - t0, 1)
  rec$max_mem_mb <- round(sum(gc()[, 6]), 0)
  log[[i]] <- as.data.table(rec)
}
L <- rbindlist(log, fill = TRUE)
fwrite(L, file.path(OUT, "coloc_pilot_v3_technical_log.csv"))
saveRDS(susie_keep, file.path(OUT, "coloc_pilot_v3_susie_crea.rds")); saveRDS(sealed, file.path(OUT, "coloc_pilot_v3_SEALED_results.rds"))   # not opened during the pilot
si <- sessionInfo()
writeLines(c(R.version.string, sprintf("susieR %s", packageVersion("susieR")), sprintf("coloc %s", packageVersion("coloc")),
             sprintf("data.table %s", packageVersion("data.table")), "seed 20261008"), file.path(OUT, "coloc_pilot_v3_versions.txt"))
print(L[, !c("warnings"), with = FALSE])
cat("\nWarnings:\n"); print(L[, .(locus, warnings)])
