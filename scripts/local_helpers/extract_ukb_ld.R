# Extract UK Biobank LD (Weissbrod et al., PolyFun; GRCh37, British ancestry, CC BY 4.0) for the 10 coloc-pilot loci
# For each locus: download one 3-Mb region (.gz SNP list + ~1.2 GB .npz sparse matrix), check it,
# keep only the pilot SNPs inside the locus, save a small file, then delete the big download.
# Needs in ~/Desktop/GWAS: coloc_pilot_rsids.txt, pilot_ukb_windows.csv. Uses only base R + data.table.
# Safe to re-run: finished loci are skipped. Temporary disk space needed: about 6 GB.

library(data.table)
options(timeout = 7200)
URL <- "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/"
win <- fread("pilot_ukb_windows.csv", colClasses = list(character = "chr"))
ids <- readLines("coloc_pilot_rsids.txt")
dir.create("ukb_ld_pilot", showWarnings = FALSE)

## ---- minimal .npy reader (seek-based, so the big data array is never fully loaded) ----
npy_open <- function(path) {
  con <- file(path, "rb")
  magic <- readBin(con, "raw", 6); stopifnot(identical(magic[2:6], charToRaw("NUMPY")))
  ver <- as.integer(readBin(con, "raw", 2))
  hlen <- if (ver[1] == 1) readBin(con, "integer", 1, size = 2, signed = FALSE, endian = "little")
          else readBin(con, "integer", 1, size = 4, endian = "little")
  hdr <- rawToChar(readBin(con, "raw", hlen))
  descr <- sub(".*'descr': *'([^']+)'.*", "\\1", hdr)
  shp <- sub(".*'shape': *\\(([^)]*)\\).*", "\\1", hdr)
  n <- prod(as.numeric(strsplit(gsub(" ", "", shp), ",")[[1]]))
  list(con = con, descr = descr, n = ifelse(is.na(n), 1, n), offset = seek(con))
}
npy_read <- function(h, start = 0, count = h$n) {           # start is 0-based element index
  sz <- as.integer(sub("^[<|=]?[a-zA-Z]", "", h$descr)); ty <- sub("^[<|=]?([a-zA-Z]).*", "\\1", h$descr)
  seek(h$con, h$offset + start * sz)
  if (ty == "f") return(readBin(h$con, "numeric", count, size = sz, endian = "little"))
  if (ty %in% c("i", "u") && sz <= 4) return(readBin(h$con, "integer", count, size = sz, endian = "little"))
  if (ty %in% c("i", "u") && sz == 8) {
    v <- readBin(h$con, "integer", 2 * count, size = 4, endian = "little")
    lo <- v[c(TRUE, FALSE)]; hi <- v[c(FALSE, TRUE)]
    return(ifelse(lo < 0, lo + 2^32, lo) + hi * 2^32)
  }
  if (ty == "S") return(rawToChar(readBin(h$con, "raw", sz)))
  stop("Unsupported dtype ", h$descr)
}

for (k in seq_len(nrow(win))) {
  W <- win[k]; out <- file.path("ukb_ld_pilot", paste0(W$locus, ".rds"))
  if (file.exists(out)) { message("Already done: ", W$locus); next }
  rng <- as.numeric(strsplit(W$locus, "_")[[1]][2:3])
  ok <- FALSE
  for (pref in strsplit(W$ukb_windows, ";")[[1]]) {
    # reuse a complete earlier download of this same window (e.g. after an interrupted run)
    if (file.exists("tmp_ukb.gz") && file.exists("tmp_ukb.npz")) {
      chk <- tryCatch(fread("tmp_ukb.gz"), error = function(e) NULL)
      wp <- as.numeric(strsplit(pref, "_")[[1]][2:3]); wc <- sub("^chr", "", strsplit(pref, "_")[[1]][1])
      zipok <- tryCatch(nrow(unzip("tmp_ukb.npz", list = TRUE)) >= 4, error = function(e) FALSE)
      if (!is.null(chk) && zipok && all(as.character(chk$chromosome) == wc) &&
          min(chk$position) >= wp[1] && max(chk$position) <= wp[2]) {
        message(Sys.time(), "  ", W$locus, ": reusing downloaded ", pref); ok <- TRUE; break
      }
    }
    message(Sys.time(), "  ", W$locus, ": trying ", pref)
    g <- tryCatch({ download.file(paste0(URL, pref, ".gz"), "tmp_ukb.gz", mode = "wb", quiet = TRUE); TRUE }, error = function(e) FALSE)
    if (!g) next
    n <- tryCatch({ download.file(paste0(URL, pref, ".npz"), "tmp_ukb.npz", mode = "wb"); TRUE }, error = function(e) FALSE)
    if (n) { ok <- TRUE; break }
  }
  if (!ok) { message("  Could not download any window for ", W$locus, " - skipped"); next }

  snps <- fread("tmp_ukb.gz")
  sel <- which(snps$rsid %in% ids & snps$position >= rng[1] & snps$position <= rng[2])
  message(sprintf("  %d SNPs in window; %d pilot SNPs inside the locus", nrow(snps), length(sel)))
  md5 <- unname(tools::md5sum("tmp_ukb.npz"))
  unzip("tmp_ukb.npz", exdir = "tmp_npz")
  fmt <- { h <- npy_open("tmp_npz/format.npy"); x <- npy_read(h); close(h$con); x }
  shp <- { h <- npy_open("tmp_npz/shape.npy"); x <- npy_read(h); close(h$con); x }
  if (!fmt %in% c("csr", "csc", "coo")) stop("Unexpected sparse format: ", fmt)
  if (shp[1] != nrow(snps)) stop("Matrix size does not match SNP list")
  S <- matrix(0, length(sel), length(sel)); pos <- match(seq_len(nrow(snps)), sel)
  if (fmt %in% c("csr", "csc")) {
    hp <- npy_open("tmp_npz/indptr.npy"); ptr <- npy_read(hp); close(hp$con)
    hi <- npy_open("tmp_npz/indices.npy"); hd <- npy_open("tmp_npz/data.npy")
    for (a in seq_along(sel)) {                                 # one stored row (csr) or column (csc) at a time
      r <- sel[a]; st <- ptr[r]; cnt <- ptr[r + 1] - ptr[r]
      if (cnt == 0) next
      j <- npy_read(hi, st, cnt) + 1; v <- npy_read(hd, st, cnt)
      b <- pos[j]; keep <- !is.na(b)
      S[a, b[keep]] <- S[a, b[keep]] + v[keep]
    }
    close(hi$con); close(hd$con)
  } else {                                                      # coo: scan row/col/data in chunks of 10 million
    hr <- npy_open("tmp_npz/row.npy"); hc <- npy_open("tmp_npz/col.npy"); hd <- npy_open("tmp_npz/data.npy")
    nnz <- hr$n; step <- 1e7
    for (st in seq(0, nnz - 1, by = step)) {
      cnt <- min(step, nnz - st)
      ri <- pos[npy_read(hr, st, cnt) + 1]; ci <- pos[npy_read(hc, st, cnt) + 1]
      keep <- which(!is.na(ri) & !is.na(ci))
      if (length(keep)) {
        v <- npy_read(hd, st, cnt)[keep]
        idx <- cbind(ri[keep], ci[keep])
        S[idx] <- S[idx] + v
      }
      message(sprintf("    scanned %.0f%%", 100 * (st + cnt) / nnz))
    }
    close(hr$con); close(hc$con); close(hd$con)
  }
  R <- S + t(S)                                                 # stored as a triangle with diagonal 0.5
  info <- list(locus = W$locus, window = pref, npz_md5 = md5, format = fmt,
               diag_range = range(diag(R)), symmetric = isSymmetric(R), any_na = anyNA(R))
  saveRDS(list(R = R, snps = snps[sel], info = info), out)
  message(sprintf("  saved %s (diag %.3f-%.3f, NA: %s)", out, info$diag_range[1], info$diag_range[2], info$any_na))
  unlink(c("tmp_ukb.gz", "tmp_ukb.npz", "tmp_npz"), recursive = TRUE)
}
message("Done. Please upload the 10 .rds files in the ukb_ld_pilot folder.")
