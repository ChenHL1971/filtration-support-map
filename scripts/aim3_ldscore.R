suppressPackageStartupMessages({library(data.table); library(readxl)})
s8 <- as.data.table(read_excel("raw/science.adp4753_tables_s1_to_s31.xlsx", sheet = "S8 Finemap.Multi", skip = 1))
s8 <- s8[!is.na(RSID) & grepl("^chr", `Merged GWAS loci`)]
s8[, `:=`(signal_id = paste0(`Merged GWAS loci`, "|", `CARMA Credible Set`), pip = as.numeric(`CARMA posterior inclusion probability (PIP)`),
          pos = as.numeric(Position), chr = as.character(Chromosome))]
top <- s8[order(-pip)][!duplicated(signal_id), .(signal_id, chr, pos, rsid = RSID, cs_size = 0L)]
top[, cs_size := s8[, .N, by = signal_id][match(top$signal_id, signal_id), N]]
REF <- "raw/ldfull/full_ld_ref"
bim <- fread(paste0(REF, ".bim"), header = FALSE, col.names = c("chr","rsid","cm","pos","a1","a2")); bim[, chr := as.character(chr)]
nind <- 503; nb <- ceiling(nind/4); con <- file(paste0(REF, ".bed"), "rb"); BED <- readBin(con, "raw", 3 + nb*nrow(bim)); close(con)
lut <- c(0, NA, 1, 2)
geno <- function(idx) vapply(idx, function(j) { b <- as.integer(BED[3 + (j-1)*nb + seq_len(nb)]); lut[as.vector(rbind(b%%4,(b%/%4)%%4,(b%/%16)%%4,(b%/%64)%%4))[seq_len(nind)] + 1] }, numeric(nind))
top[, ldscore := NA_real_]
for (i in seq_len(nrow(top))) {
  j <- bim[chr == top$chr[i] & pos == top$pos[i], which = TRUE]; if (!length(j)) next
  w <- bim[chr == top$chr[i] & abs(pos - top$pos[i]) <= 5e5, which = TRUE]
  G <- geno(c(j[1], w)); x <- G[, 1]; G <- G[, -1, drop = FALSE]
  ok <- !is.na(x); r <- suppressWarnings(cor(x[ok], G[ok, , drop = FALSE], use = "pairwise.complete.obs"))
  top$ldscore[i] <- sum(r^2, na.rm = TRUE)
}
fwrite(top, "out/aim3_signal_leads.csv"); cat("leads", nrow(top), " with ldscore", sum(!is.na(top$ldscore)), "\n")
