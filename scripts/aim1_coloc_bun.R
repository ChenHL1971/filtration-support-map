# Aim 1: signal-level colocalization of eGFRcrea vs eGFRcys at all S8 merged loci (Plan-5; rules logged 2026-10-10)
# LD: 1000 Genomes EUR (MAGMA g1000_eur). PP.H4 and directions are written to per-locus SEALED files, never printed.
suppressPackageStartupMessages({library(data.table); library(susieR); library(coloc); library(readxl); library(parallel)})
set.seed(20261008)
args <- commandArgs(TRUE); GDIR <- args[1]; REF <- args[2]; SUPP <- args[3]; OUT <- args[4]; CORES <- as.integer(args[5])
dir.create(file.path(OUT, "sealed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "tech"), recursive = TRUE, showWarnings = FALSE)
MAXV <- 4000; FLANK <- 250000

## ---- reference genotypes: whole .bed held as raw bytes, columns decoded per locus ----
bim <- fread(paste0(REF, ".bim"), header = FALSE, col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
nind <- nrow(fread(paste0(REF, ".fam"), header = FALSE)); nb <- ceiling(nind / 4)
con <- file(paste0(REF, ".bed"), "rb"); BED <- readBin(con, "raw", 3 + nb * nrow(bim)); close(con)
stopifnot(identical(BED[1:3], as.raw(c(0x6c, 0x1b, 0x01))))
lut <- c(0, NA, 1, 2)
geno <- function(idx) vapply(idx, function(j) {
  b <- as.integer(BED[3 + (j - 1) * nb + seq_len(nb)])
  lut[as.vector(rbind(b %% 4, (b %/% 4) %% 4, (b %/% 16) %% 4, (b %/% 64) %% 4))[seq_len(nind)] + 1]
}, numeric(nind))

## ---- GWAS ----
comp <- c(A = "T", T = "A", C = "G", G = "C")
liu <- fread(paste0(GDIR, "/liu_crea.gz"), select = c("MarkerName", "CHR", "POS", "REF", "ALT", "MAF", "Zscore", "N"))
setnames(liu, c("rsid", "chr", "pos", "nea", "ea", "maf", "z", "n")); liu[, `:=`(chr = as.character(chr), ea = toupper(ea), nea = toupper(nea))]
cys <- fread(paste0(GDIR, "/bun.gz"), select = c("RSID", "Allele1", "Allele2", "n", "Freq1", "Effect", "P.value", "Direction"))
cys[, sgn := fifelse(Effect > 0, 1, fifelse(Effect < 0, -1,
             fifelse(grepl("^\\++$", Direction), 1, fifelse(grepl("^-+$", Direction), -1, NA_real_))))]
cys[, z := sgn * qnorm(P.value / 2, lower.tail = FALSE)]
cys <- cys[is.finite(z), .(rsid = RSID, ea = toupper(Allele1), nea = toupper(Allele2), maf = pmin(Freq1, 1 - Freq1), z, n)]
setkey(liu, rsid); setkey(cys, rsid)

orient <- function(ea, nea, A1, A2) {   # +1 if (ea,nea)==(A2,A1); exact match/swap first; strand flip only for non-palindromic
  same <- ea == A2 & nea == A1; swap <- ea == A1 & nea == A2
  pal <- paste0(A1, A2) %in% c("AT", "TA", "CG", "GC")
  fs <- !pal & comp[ea] == A2 & comp[nea] == A1; fw <- !pal & comp[ea] == A1 & comp[nea] == A2
  same[is.na(same)] <- FALSE; swap[is.na(swap)] <- FALSE; fs[is.na(fs)] <- FALSE; fw[is.na(fw)] <- FALSE
  fifelse(same | fs, 1, fifelse(swap | fw, -1, NA_real_))
}
align <- function(d, b) {
  m <- merge(d, b[, .(rsid, A1 = a1, A2 = a2, bidx)], by = "rsid")
  pal <- paste0(m$A1, m$A2) %in% c("AT", "TA", "CG", "GC") & m$maf >= 0.42
  m[, z := z * orient(ea, nea, A1, A2)]; m[!pal & is.finite(z) & maf > 0]
}
std <- function(d, snps, R) { v <- 1 / (2 * d$maf * (1 - d$maf) * (d$n + d$z^2))
  list(beta = setNames(d$z * sqrt(v), snps), varbeta = setNames(v, snps), N = median(d$n), sdY = 1,
       type = "quant", snp = snps, LD = R, MAF = setNames(d$maf, snps)) }
abf_lbf <- function(D) { r <- 0.15^2 / (0.15^2 + D$varbeta); 0.5 * (log(1 - r) + r * D$beta^2 / D$varbeta) }

## ---- S8 signals ----
s8 <- as.data.table(read_excel(SUPP, sheet = "S8 Finemap.Multi", skip = 1))
s8 <- s8[!is.na(RSID) & grepl("^chr", `Merged GWAS loci`)]
s8[, `:=`(locus = `Merged GWAS loci`, signal_id = paste0(`Merged GWAS loci`, "|", `CARMA Credible Set`),
          pip = as.numeric(`CARMA posterior inclusion probability (PIP)`), pos = as.numeric(Position), chr = as.character(Chromosome))]
loci <- unique(s8$locus); if (length(args) >= 6) loci <- intersect(loci, strsplit(args[6], ",")[[1]])

run_locus <- function(loc) {
  tf <- file.path(OUT, "tech", paste0(loc, ".rds")); if (file.exists(tf)) return(invisible(NULL))
  t0 <- proc.time()[["elapsed"]]; tech <- list(locus = loc, status = "ok", note = ""); sig_rows <- list(); sealed <- list()
  res <- tryCatch({
    sg <- s8[locus == loc]; ch <- sg$chr[1]
    rng <- as.numeric(strsplit(loc, "_")[[1]][2:3])
    lo <- max(rng[1], min(sg$pos) - FLANK); hi <- min(rng[2], max(sg$pos) + FLANK)
    b <- bim[chr == ch & pos >= lo & pos <= hi]; b[, bidx := match(rsid, bim$rsid)]
    a1 <- align(liu[chr == ch & pos >= lo & pos <= hi], b); c1 <- align(cys[rsid %in% b$rsid], b)
    snps <- intersect(a1$rsid, c1$rsid)
    a1 <- a1[match(snps, rsid)]; c1 <- c1[match(snps, rsid)]
    G <- geno(a1$bidx); okv <- apply(G, 2, function(x) mean(is.na(x)) < 0.05 && sd(x, na.rm = TRUE) > 0)
    if (sum(okv) > MAXV) { mf <- colMeans(G, na.rm = TRUE) / 2; okv <- okv & pmin(mf, 1 - mf) >= 0.01; tech$maf_filter_applied <- TRUE }
    if (sum(okv) > MAXV) { keepi <- which(okv); okv[] <- FALSE; okv[keepi[order(-abs(a1$z[keepi]))[seq_len(MAXV)]]] <- TRUE; tech$top_z_cap_applied <- TRUE }
    snps <- snps[okv]; a1 <- a1[okv]; c1 <- c1[okv]; G <- G[, okv, drop = FALSE]
    G <- apply(G, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
    R <- cor(G); dimnames(R) <- list(snps, snps); rm(G)
    tech$region <- paste0(ch, ":", lo, "-", hi); tech$variants_used <- length(snps)
    tech$s_liu <- signif(estimate_s_rss(a1$z, R, n = median(a1$n)), 3); tech$s_cys <- signif(estimate_s_rss(c1$z, R, n = median(c1$n)), 3)
    D1 <- std(a1, snps, R); D2 <- std(c1, snps, R)
    run <- function(D) tryCatch(suppressWarnings(suppressMessages(runsusie(D, repeat_until_convergence = FALSE, maxit = 1000, L = 10))), error = function(e) e)
    S1 <- run(D1); S2 <- run(D2)
    tech$conv_liu <- !inherits(S1, "error") && isTRUE(S1$converged); tech$conv_cys <- !inherits(S2, "error") && isTRUE(S2$converged)
    k1 <- if (tech$conv_liu) length(S1$sets$cs) else 0; k2 <- if (tech$conv_cys) length(S2$sets$cs) else 0
    tech$cs_liu <- k1; tech$cs_cys <- k2
    ## SuSiE-level colocalization (sealed)
    if (k1 > 0 && k2 > 0) { invisible(capture.output(C <- coloc.susie(S1, S2))); sealed$susie <- C$summary; sealed$method <- "coloc.susie" }
    else if (k1 > 0) { bf1 <- S1$lbf_variable[S1$sets$cs_index, , drop = FALSE]; colnames(bf1) <- snps
      invisible(capture.output(C <- coloc.bf_bf(bf1, setNames(abf_lbf(D2), snps)))); sealed$susie <- C$summary; sealed$method <- "coloc.bf_bf" }
    alpha <- if (k1 > 0) S1$alpha[S1$sets$cs_index, , drop = FALSE] else NULL
    leads <- if (k1 > 0) snps[apply(alpha, 1, which.max)] else character(0)
    sealed$lead_z <- data.table(cs = seq_along(leads), lead = leads, z_crea = a1$z[match(leads, snps)], z_cys = c1$z[match(leads, snps)])
    ## map each S8 signal; fallback = coloc.abf within its credible set
    for (sid in unique(sg$signal_id)) {
      v <- sg[signal_id == sid]; top <- v$RSID[which.max(v$pip)]
      cov <- mean(v$RSID %in% snps); row <- list(locus = loc, signal_id = sid, cs_size = nrow(v), cs_coverage = round(cov, 3))
      m <- "unmatched"; kk <- NA_integer_
      if (k1 > 0) {
        ov <- apply(alpha, 1, function(al) sum(al[snps %in% v$RSID])); k <- which.max(ov); row$overlap <- round(ov[k], 3)
        if (ov[k] >= 0.5) { m <- "overlap"; kk <- k } else if (top %in% snps) {
          rr <- R[top, leads]^2; j <- which.max(rr); row$r2_lead <- round(rr[j], 3)
          if (rr[j] >= 0.8) { m <- "LD"; kk <- j } }
      }
      if (m == "unmatched") {
        inc <- snps %in% v$RSID
        if (cov >= 0.8 && sum(inc) >= 1) {
          sub <- function(D) list(beta = D$beta[inc], varbeta = D$varbeta[inc], N = D$N, sdY = 1, type = "quant", snp = snps[inc], MAF = D$MAF[inc])
          invisible(capture.output(A <- coloc.abf(sub(D1), sub(D2)))); m <- "fallback_abf"
          zt <- if (top %in% snps) c(a1$z[snps == top], c1$z[snps == top]) else c(NA, NA)
          sealed$fallback[[sid]] <- list(summary = A$summary, top = top, z_crea = zt[1], z_cys = zt[2])
        } else m <- "missing_coverage"
      }
      row$method <- m; row$susie_cs <- kk; row$flag_small_cs <- nrow(v) < 5
      sig_rows[[length(sig_rows) + 1]] <- as.data.table(row)
    }
    TRUE
  }, error = function(e) { tech$status <<- "error"; tech$note <<- substr(conditionMessage(e), 1, 200); FALSE })
  tech$runtime_sec <- round(proc.time()[["elapsed"]] - t0, 1)
  saveRDS(sealed, file.path(OUT, "sealed", paste0(loc, ".rds")))
  saveRDS(list(tech = as.data.table(tech), signals = rbindlist(sig_rows, fill = TRUE)), tf)
  invisible(NULL)
}

message(Sys.time(), " start: ", length(loci), " loci, ", CORES, " cores")
invisible(mclapply(loci, function(l) { run_locus(l); gc(); NULL }, mc.cores = CORES, mc.preschedule = FALSE))
message(Sys.time(), " done")
