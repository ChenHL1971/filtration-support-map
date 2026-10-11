# LD crosswalk for scorecard lead variants (feasibility audit, Plan-5 locked)
# For each of the 429 lead variants not found directly in an S8 credible set,
# compute r^2 (1000 Genomes EUR, GRCh37) with every credible-set variant within 500 kb.
# Needs: ld_query_input.csv in the working folder, and my_token set in the console first.
# Safe to re-run: finished lead variants are skipped.

if (!requireNamespace("LDlinkR", quietly = TRUE)) install.packages("LDlinkR")
library(LDlinkR)
library(data.table)

if (!exists("my_token") || nchar(my_token) < 10) stop("Please set my_token <- \"your token\" in the console first")

inp      <- fread("ld_query_input.csv")
out_file <- "ld_crosswalk_results.csv"
done     <- if (file.exists(out_file)) unique(fread(out_file)$lead_variant) else character(0)
leads    <- setdiff(unique(inp$lead_variant), done)
message(length(leads), " lead variants to query; ", length(done), " already done")

for (i in seq_along(leads)) {
  ld   <- leads[i]
  sub  <- inp[lead_variant == ld]
  snps <- unique(c(ld, sub$cs_rsid))
  res  <- tryCatch(
    suppressWarnings(LDmatrix(snps = snps, pop = "EUR", r2d = "r2",
                              token = my_token, genome_build = "grch37")),
    error = function(e) e)

  rec <- sub[, .(lead_variant, signal_id, cs_rsid, r2 = NA_real_, note = "")]
  if (inherits(res, "error")) {
    rec[, note := paste("error:", gsub("[\r\n,]", " ", conditionMessage(res)))]
  } else if (!is.data.frame(res) || !(ld %in% res$RS_number)) {
    rec[, note := "lead_not_in_1000G"]
  } else {
    row <- res[res$RS_number == ld, , drop = FALSE]
    for (k in seq_len(nrow(rec))) {
      v <- rec$cs_rsid[k]
      if (v == ld) {
        rec[k, `:=`(r2 = 1, note = "same_rsid")]
      } else if (v %in% names(row)) {
        rec[k, r2 := suppressWarnings(as.numeric(row[[v]][1]))]
      } else {
        rec[k, note := "cs_variant_not_in_1000G"]
      }
    }
  }
  fwrite(rec, out_file, append = file.exists(out_file))
  if (i %% 10 == 0) message(sprintf("%d / %d done", i, length(leads)))
  Sys.sleep(1)   # LDlink API: one request at a time
}
message("Finished. Please upload ld_crosswalk_results.csv")
