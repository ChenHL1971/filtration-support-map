# Download MGI_GenePheno.rpt and keep only rows with a renal-function phenotype (MP:0002136 and its 81 descendants)
setwd("~/Desktop/GWAS")
options(timeout = 3600)
url <- "https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt"
if (!file.exists("MGI_GenePheno.rpt")) download.file(url, "MGI_GenePheno.rpt", mode = "wb")
terms <- readLines("mgi_renal_function_terms.txt")
x <- readLines("MGI_GenePheno.rpt")
message("MGI_GenePheno.rpt rows: ", length(x))
keep <- x[grepl(paste(terms, collapse = "|"), x)]
writeLines(keep, "mgi_genepheno_renal_function.txt")
message("Kept ", length(keep), " rows -> mgi_genepheno_renal_function.txt (please upload this file)")
