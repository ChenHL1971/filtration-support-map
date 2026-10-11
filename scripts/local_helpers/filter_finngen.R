# Filter FinnGen R13 summary statistics to the kidney-locus variants (feasibility audit, Plan-5 locked)
# FinnGen R13 is GRCh38, while our windows are hg19, so the filter matches on rsID, not position.
# Needs in the working folder: kidney_loci_rsids.txt.gz and the 5 downloaded finngen_R13_*.gz files.
# The script only keeps rows and counts them; it prints no association results.

library(data.table)

keep_ids <- fread("kidney_loci_rsids.txt.gz")$rsid
endpoints <- c("N14_DIALYSIS", "Z21_KIDNEY_TRANSPLANT_STATUS", "N14_CHRONKIDNEYDIS",
               "N14_ACUTERENFAIL", "H7_VITROPACIT")

for (ep in endpoints) {
  infile  <- sprintf("finngen_R13_%s.gz", ep)
  outfile <- sprintf("fg_%s_kidneyloci.tsv.gz", ep)
  if (!file.exists(infile)) { message("Missing: ", infile); next }
  if (file.exists(outfile)) { message("Already done: ", outfile); next }
  dt <- fread(infile)
  if (!"rsids" %in% names(dt)) stop("No 'rsids' column in ", infile, ". Columns: ", paste(names(dt), collapse = ", "))
  first_id <- sub(",.*", "", dt$rsids)          # some rows list several rsIDs separated by commas
  hit <- first_id %in% keep_ids
  multi <- which(!hit & grepl(",", dt$rsids))
  if (length(multi)) hit[multi] <- vapply(strsplit(dt$rsids[multi], ","), function(x) any(x %in% keep_ids), logical(1))
  keep <- dt[hit]
  fwrite(keep, outfile, sep = "\t", compress = "gzip")
  message(sprintf("%s: kept %d / %d variants -> %s (%.0f MB)", ep, nrow(keep), nrow(dt), outfile, file.size(outfile) / 1e6))
  rm(dt, keep); gc()
}
message("Done. Please upload the 5 fg_*_kidneyloci.tsv.gz files.")
