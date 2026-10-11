# Extract the coloc-pilot variants from the 1000 Genomes EUR reference (MAGMA g1000_eur, GRCh37)
# Needs in ~/Desktop/GWAS: coloc_pilot_rsids.txt, and the unzipped g1000_eur.bed/.bim/.fam
# Output: pilot_ld_ref.bed/.bim/.fam (small; upload these three files)

if (!requireNamespace("bigsnpr", quietly = TRUE)) install.packages("bigsnpr")
library(bigsnpr)

prefix <- list.files(pattern = "\\.bed$", recursive = TRUE)
prefix <- prefix[grepl("g1000_eur", prefix)][1]
if (is.na(prefix)) stop("Cannot find g1000_eur.bed - please unzip g1000_eur.zip inside ~/Desktop/GWAS")
prefix <- sub("\\.bed$", "", prefix)

ids <- readLines("coloc_pilot_rsids.txt")
bim <- data.table::fread(paste0(prefix, ".bim"), header = FALSE)
keep <- which(bim$V2 %in% ids)
message(sprintf("Found %d of %d pilot variants in the reference", length(keep), length(ids)))

for (ext in c(".bk", ".rds")) if (file.exists(paste0("pilot_tmp", ext))) file.remove(paste0("pilot_tmp", ext))
rds <- snp_readBed2(paste0(prefix, ".bed"), backingfile = "pilot_tmp", ind.col = keep)
obj <- snp_attach(rds)
snp_writeBed(obj, "pilot_ld_ref.bed")
file.remove(c("pilot_tmp.bk", "pilot_tmp.rds"))
message("Done. Please upload pilot_ld_ref.bed, pilot_ld_ref.bim, pilot_ld_ref.fam")
