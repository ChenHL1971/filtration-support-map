# Aim 1 coloc-SuSiE pilot v2: UK Biobank LD (Weissbrod), rules fixed 2026-10-10 before data arrived. PP.H4 sealed.
suppressPackageStartupMessages({library(data.table); library(susieR); library(coloc)})
set.seed(20261008)
args <- commandArgs(TRUE); GDIR <- args[1]; UKB <- args[2]; REF <- args[3]; LOCI <- args[4]; OUT <- args[5]
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
ref <- read_bed(REF)
comp <- c(A = "T", T = "A", C = "G", G = "C")
liu <- fread(paste0(GDIR, "/liu_crea.gz"), select = c("MarkerName", "REF", "ALT", "MAF", "Zscore", "N"))
setnames(liu, c("rsid", "nea", "ea", "maf", "z", "n")); liu[, `:=`(ea = toupper(ea), nea = toupper(nea))]
cys <- fread(paste0(GDIR, "/cys.gz"), select = c("RSID", "Allele1", "Allele2", "n", "Freq1", "Effect", "P.value", "Direction"))
cys[, sgn := fifelse(Effect > 0, 1, fifelse(Effect < 0, -1,
             fifelse(grepl("^\\++$", Direction), 1, fifelse(grepl("^-+$", Direction), -1, NA_real_))))]
cys[, z := sgn * qnorm(P.value / 2, lower.tail = FALSE)]
cys <- cys[is.finite(z), .(rsid = RSID, ea = toupper(Allele1), nea = toupper(Allele2), maf = pmin(Freq1, 1 - Freq1), z, n)]

orient <- function(ea, nea, A1, A2) {     # +1 if (ea,nea) == (A2,A1), -1 if swapped; exact allele match is tested first,
  same <- ea == A2 & nea == A1; swap <- ea == A1 & nea == A2        # strand flip only for non-palindromic SNPs
  pal <- paste0(A1, A2) %in% c("AT", "TA", "CG", "GC")
  fs <- !pal & comp[ea] == A2 & comp[nea] == A1; fw <- !pal & comp[ea] == A1 & comp[nea] == A2
  same[is.na(same)] <- FALSE; swap[is.na(swap)] <- FALSE; fs[is.na(fs)] <- FALSE; fw[is.na(fw)] <- FALSE
  fifelse(same | fs, 1, fifelse(swap | fw, -1, NA_real_))
}
align <- function(d, snps) {              # z oriented to UKB allele2
  m <- merge(d, snps[, .(rsid, A1 = allele1, A2 = allele2)], by = "rsid")
  pal <- paste0(m$A1, m$A2) %in% c("AT", "TA", "CG", "GC") & m$maf >= 0.42
  f <- orient(m$ea, m$nea, m$A1, m$A2)
  m[, z := z * f]; m[!pal & is.finite(z)]
}
std_beta <- function(z, maf, n) { v <- 1 / (2 * maf * (1 - maf) * (n + z^2)); list(beta = z * sqrt(v), varbeta = v) }
mk <- function(d, R, snps) { b <- std_beta(d$z, d$maf, d$n)
  list(beta = setNames(b$beta, snps), varbeta = setNames(b$varbeta, snps), N = median(d$n), sdY = 1,
       type = "quant", snp = snps, LD = R, MAF = setNames(d$maf, snps)) }

loci <- fread(LOCI, colClasses = list(character = "chr"))
log <- list(); sealed <- list()
for (i in seq_len(nrow(loci))) {
  L <- loci[i]; t0 <- proc.time()[["elapsed"]]; message(Sys.time(), " start ", L$locus)
  rec <- list(locus = L$locus, n_signals_S8 = L$n_signals, status = "ok", note = "")
  tryCatch({
    U <- readRDS(file.path(UKB, paste0(L$locus, ".rds"))); R0 <- U$R; us <- as.data.table(U$snps)
    dup <- duplicated(us$rsid); R0 <- R0[!dup, !dup]; us <- us[!dup]
    a <- align(liu, us); b <- align(cys, us)
    snps <- intersect(intersect(us$rsid, a$rsid), b$rsid)
    ix <- match(snps, us$rsid); R <- R0[ix, ix]; dimnames(R) <- list(snps, snps)
    a <- a[match(snps, rsid)]; b <- b[match(snps, rsid)]
    rec$ukb_snps_in_locus <- nrow(us); rec$variants_used <- length(snps)
    ## rule 1: orientation check against 1000 Genomes (pairs with |r| >= 0.3 in both)
    kb <- match(snps, ref$bim$rsid); okb <- !is.na(kb)
    fk <- orient(ref$bim$a2[kb[okb]], ref$bim$a1[kb[okb]], us$allele1[ix][okb], us$allele2[ix][okb])
    X <- ref$G[, kb[okb], drop = FALSE]; X <- apply(X, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
    keepk <- apply(X, 2, sd) > 0 & !is.na(fk)
    if (sum(keepk) > 50) {
      Rk <- cor(X[, keepk]); Rk <- Rk * outer(fk[keepk], fk[keepk])
      Ru <- R[okb, okb][keepk, keepk]; ut <- upper.tri(Rk) & abs(Rk) >= 0.3 & abs(Ru) >= 0.3
      rec$orientation_pairs <- sum(ut); rec$orientation_sign_agree_pct <- round(100 * mean(sign(Rk[ut]) == sign(Ru[ut])), 2)
    }
    if (!is.null(rec$orientation_sign_agree_pct) && rec$orientation_sign_agree_pct < 95) stop("orientation check failed")
    ## LD mismatch before QC
    na <- median(a$n); nb <- median(b$n)
    rec$s_liu_before <- signif(estimate_s_rss(a$z, R, n = na), 3); rec$s_cys_before <- signif(estimate_s_rss(b$z, R, n = nb), 3)
    ## rule 2: kriging_rss outliers removed once
    ka <- kriging_rss(a$z, R, n = na)$conditional_dist; kb2 <- kriging_rss(b$z, R, n = nb)$conditional_dist
    out <- (ka$logLR > 2 & abs(a$z) > 2) | (kb2$logLR > 2 & abs(b$z) > 2)
    rec$qc_removed <- sum(out)
    snps <- snps[!out]; R <- R[!out, !out]; a <- a[!out]; b <- b[!out]
    rec$s_liu_after <- signif(estimate_s_rss(a$z, R, n = na), 3); rec$s_cys_after <- signif(estimate_s_rss(b$z, R, n = nb), 3)
    ## rule 3: SuSiE
    D1 <- mk(a, R, snps); D2 <- mk(b, R, snps); check_dataset(D1, req = "LD"); check_dataset(D2, req = "LD")
    w <- character(0); hw <- function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") }
    run <- function(D) tryCatch(withCallingHandlers(runsusie(D, repeat_until_convergence = FALSE, maxit = 1000, L = 10,
                                 estimate_residual_variance = FALSE), warning = hw), error = function(e) e)
    S1 <- run(D1); S2 <- run(D2)
    rec$converged_liu <- !inherits(S1, "error") && isTRUE(S1$converged)
    rec$converged_cys <- !inherits(S2, "error") && isTRUE(S2$converged)
    rec$cs_liu <- if (inherits(S1, "error")) NA else length(S1$sets$cs)
    rec$cs_cys <- if (inherits(S2, "error")) NA else length(S2$sets$cs)
    rec$warnings <- paste(unique(substr(w, 1, 100)), collapse = " | ")
    if (rec$converged_liu && rec$converged_cys && rec$cs_liu > 0 && rec$cs_cys > 0) {
      C <- coloc.susie(S1, S2); rec$coloc_method <- "coloc.susie"; rec$coloc_pairs <- nrow(C$summary); sealed[[L$locus]] <- C$summary
    } else if (L$n_signals == 1 && !(rec$converged_liu && rec$converged_cys)) {
      C <- coloc.abf(D1, D2); rec$coloc_method <- "coloc.abf (fallback)"; rec$coloc_pairs <- 1; sealed[[L$locus]] <- C$summary
    } else {
      rec$coloc_method <- "none"; rec$coloc_pairs <- 0
      rec$note <- if (!(rec$converged_liu && rec$converged_cys)) "non-convergence at multi-signal locus (recorded as failed)" else "no credible set in at least one trait"
    }
    ## rule 5: technical pass
    rec$technical_pass <- rec$s_liu_after < 0.3 && rec$s_cys_after < 0.3 && rec$converged_liu && rec$converged_cys
  }, error = function(e) { rec$status <<- "error"; rec$note <<- conditionMessage(e) })
  rec$runtime_sec <- round(proc.time()[["elapsed"]] - t0, 1)
  log[[i]] <- as.data.table(rec)
}
L <- rbindlist(log, fill = TRUE)
fwrite(L, file.path(OUT, "coloc_pilot_v2_technical_log.csv"))
saveRDS(sealed, file.path(OUT, "coloc_pilot_v2_SEALED_results.rds"))
message("done")
