# One-off cleanup of the raw RData files handed off from the Shiny app.
# Keeps only the objects actually used by the default glmnet pipeline
# (see the plan / README for the audit). Original files are preserved
# alongside the cleaned ones with an "_original" suffix.

clean_one <- function(path, keep, backup_suffix = "_original") {
  ext <- tools::file_ext(path)
  backup_path <- sub(paste0("\\.", ext, "$"), paste0(backup_suffix, ".", ext), path)

  if (!file.exists(backup_path)) {
    file.copy(path, backup_path)
  }

  e <- new.env()
  load(backup_path, envir = e)

  missing <- setdiff(keep, ls(e))
  if (length(missing) > 0) {
    stop("Objects requested but not found in ", path, ": ", paste(missing, collapse = ", "))
  }

  dropped <- setdiff(ls(e), keep)
  message(sprintf("%s: keeping %s; dropping %s", path, paste(keep, collapse = ", "), paste(dropped, collapse = ", ")))

  save(list = keep, envir = e, file = path)
  invisible(NULL)
}

clean_one(
  "data-raw/train_test_counts_and_clin.RData",
  keep = c("a_counts", "b_counts", "gene_anns", "pheno_data")
)

clean_one(
  "data-raw/unseen_counts_for_survival_modeling.RData",
  keep = c("unseen_a_counts", "unseen_b_counts", "all_clin")
)

# independent_cohorts_counts.Rdata needs no changes: bryan_counts, meeks_counts,
# uromol_counts, and independent_clin are all used by validateExternalCohort().
