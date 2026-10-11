# Map Titan et al. 2026 ST1 lead variants (hg19) to Liu S8 signals. Rule fixed 2026-10-10 before unsealing.
suppressPackageStartupMessages({library(data.table); library(readxl)})
t <- fread("out/titan_ST1.csv"); setnames(t, c("marker","rsid","chr","pos","class","cluster","cstart","cend"))
t[, chr := as.character(chr)]
s8 <- as.data.table(read_excel("raw/science.adp4753_tables_s1_to_s31.xlsx", sheet = "S8 Finemap.Multi", skip = 1))
s8 <- s8[!is.na(RSID) & grepl("^chr", `Merged GWAS loci`)]
s8[, `:=`(signal_id = paste0(`Merged GWAS loci`, "|", `CARMA Credible Set`), pip = as.numeric(`CARMA posterior inclusion probability (PIP)`),
          pos = as.numeric(Position), chr = as.character(Chromosome))]
top <- s8[order(-pip)][!duplicated(signal_id), .(signal_id, chr, top_pos = pos, top_rsid = RSID)]
REF <- "raw/ldfull/full_ld_ref"
bim <- fread(paste0(REF, ".bim"), header = FALSE, col.names = c("chr","rsid","cm","pos","a1","a2")); bim[, chr := as.character(chr)]
nind <- 503; nb <- ceiling(nind/4); con <- file(paste0(REF, ".bed"), "rb"); BED <- readBin(con, "raw", 3 + nb*nrow(bim)); close(con)
lut <- c(0, NA, 1, 2)
g1 <- function(j) { b <- as.integer(BED[3 + (j-1)*nb + seq_len(nb)]); lut[as.vector(rbind(b%%4,(b%/%4)%%4,(b%/%16)%%4,(b%/%64)%%4))[seq_len(nind)] + 1] }
r2 <- function(i, j) { x <- g1(i); y <- g1(j); ok <- !is.na(x) & !is.na(y); if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA); cor(x[ok], y[ok])^2 }
out <- list()
for (k in seq_len(nrow(t))) {
  T <- t[k]; hit <- unique(s8[chr == T$chr & pos == T$pos, signal_id])
  if (length(hit)) { out[[k]] <- data.table(marker = T$marker, class = T$class, signal_id = hit, how = "in credible set"); next }
  cand <- top[chr == T$chr & abs(top_pos - T$pos) <= 5e5]
  bi <- bim[chr == T$chr & pos == T$pos, which = TRUE]
  if (!nrow(cand) || !length(bi)) { out[[k]] <- data.table(marker = T$marker, class = T$class, signal_id = NA, how = if (!nrow(cand)) "no S8 signal within 500 kb" else "lead not in reference"); next }
  rr <- sapply(cand$top_pos, function(p) { bj <- bim[chr == T$chr & pos == p, which = TRUE]; if (!length(bj)) NA else r2(bi[1], bj[1]) })
  if (any(rr >= 0.8, na.rm = TRUE)) out[[k]] <- data.table(marker = T$marker, class = T$class, signal_id = cand$signal_id[which(rr >= 0.8)], how = "LD r2 >= 0.8")
  else out[[k]] <- data.table(marker = T$marker, class = T$class, signal_id = NA, how = "no S8 signal in LD")
}
M <- rbindlist(out); fwrite(M, "out/titan_to_s8_map.csv")
print(M[, .(loci = uniqueN(marker)), by = .(class, mapped = !is.na(signal_id))][order(class)])
print(M[!is.na(signal_id), .(signals = uniqueN(signal_id)), by = class])
cat("signals mapped to >1 class:", M[!is.na(signal_id), uniqueN(class), by = signal_id][V1 > 1, .N], "\n")
