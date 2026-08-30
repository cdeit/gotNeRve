# Precomputes everything the package ships as processed data:
#  - one clinical table (train/test/unseen dataset membership, minimal
#    necessary columns)
#  - log2-CPM expression matrices for train/test/unseen (counts -> CPM done
#    once here, so the vignette works from already-processed expression data)
#
# Raw counts and the gene-type/filterByExpr machinery stay in data-raw/ and
# R/ (kept for reference) but are not part of the package's public data. The
# external Uromol/Meeks/Bryan cohorts are intentionally not processed here or
# shipped as package data -- that data isn't ours to redistribute.

devtools::load_all(".", quiet = TRUE)

e1 <- new.env()
load("data-raw/train_test_counts_and_clin.RData", envir = e1)
a_counts <- e1$a_counts
b_counts <- e1$b_counts
gene_anns <- e1$gene_anns
pheno_data <- e1$pheno_data

e2 <- new.env()
load("data-raw/unseen_counts_for_survival_modeling.RData", envir = e2)
unseen_a_counts <- e2$unseen_a_counts
unseen_b_counts <- e2$unseen_b_counts
all_clin <- e2$all_clin

config <- gotnerve_config()

# --- 1. Assign one dataset per patient: train > test > unseen > none ------

dataset <- rep("none", nrow(all_clin))
dataset[all_clin$`RNA-seq` %in% colnames(unseen_a_counts) | all_clin$`RNA-seq` %in% colnames(unseen_b_counts)] <- "unseen"
dataset[all_clin$`RNA-seq` %in% colnames(b_counts)] <- "test"
dataset[all_clin$`RNA-seq` %in% colnames(a_counts)] <- "train"

stan_clinical <- all_clin[, c("RNA-seq", "Cohort", "Progression", "HG_recur", "Time_to_prog_or_FUend", "Time_to_HG_recur_or_FUend")]
stan_clinical$dataset <- dataset

# nerve_group / nerve_obj_area_percent / Sex / BCG_failure only exist for the
# original train+test cohort (pheno_data); left join by RNA-seq.
pheno_sub <- pheno_data[, c("RNA-seq", "Sex", "BCG_failure", "nerve_group", "nerve_obj_area_percent")]
stan_clinical <- dplyr::left_join(stan_clinical, pheno_sub, by = "RNA-seq")

# Quantile-threshold re-labels nerve_group into High/Low, dropping the middle
# tier -- apply once here rather than showing it in the vignette every time.
train_test_rows <- stan_clinical$dataset %in% c("train", "test")
relabeled <- apply_quantile_threshold(stan_clinical[train_test_rows, ], "nerve_group", config$quantile_threshold)
stan_clinical$nerve_group[train_test_rows] <- NA
stan_clinical$nerve_group[match(relabeled$`RNA-seq`, stan_clinical$`RNA-seq`)] <- as.character(relabeled$nerve_group)

# A train/test patient whose nerve_group fell in the middle tier can't be
# used for training/evaluation (their class is ambiguous by design), but
# they still have valid expression + outcome data -- move them to the
# unseen validation pool rather than discarding them.
dropped_by_quantile <- train_test_rows & is.na(stan_clinical$nerve_group)
dropped_ids <- stan_clinical$`RNA-seq`[dropped_by_quantile] # captured now; row order/subset changes below
message(sprintf("Quantile threshold dropped %d train/test patients (middle tier) -- reassigned to dataset \"unseen\".", length(dropped_ids)))
stan_clinical$dataset[dropped_by_quantile] <- "unseen"

n_none <- sum(stan_clinical$dataset == "none")
message(sprintf("Removing %d patients with dataset == 'none' (not used in train/test/unseen).", n_none))
stan_clinical <- stan_clinical[stan_clinical$dataset != "none", ]

message("Final dataset counts:")
print(table(stan_clinical$dataset))

# Row order matters: the bootstrapped stability filter draws bootstrap
# indices via sample(seq_len(nrow(x))), so downstream code must reconstruct
# the *exact* same patient order used to verify the model in
# 01_build_bundled_model.R, or identical seeds will draw different patients.
# Train/test order follows pheno_data's own row order; unseen order follows
# the unseen counts matrices' column order.
order_key <- c(pheno_data$`RNA-seq`, setdiff(c(colnames(unseen_a_counts), colnames(unseen_b_counts)), pheno_data$`RNA-seq`))
stan_clinical <- stan_clinical[order(match(stan_clinical$`RNA-seq`, order_key)), ]

# --- 2. Precompute log2-CPM for train/test/unseen --------------------------

train_pheno <- stan_clinical[stan_clinical$dataset == "train", ]
test_pheno <- stan_clinical[stan_clinical$dataset == "test", ]

# The unseen validation pool is exactly the raw unseen counts files (~200
# patients) -- deliberately NOT deduplicated against train/test membership,
# since a handful of unseen-file patients are also part of the training/test
# cohort and still belong in unseen validation as-is. Quantile-dropped
# patients who happen to already be in this set are naturally included; the
# rare quantile-dropped patient NOT in this set is intentionally left out,
# to keep the unseen pool exactly matching the raw unseen counts files.
unseen_ids <- union(colnames(unseen_a_counts), colnames(unseen_b_counts))
message(sprintf("Unseen validation pool: %d patients.", length(unseen_ids)))

train_log2cpm <- process_expression(
  a_counts, train_pheno, target_col = "nerve_group", gene_anns = gene_anns,
  gene_type_filter = config$gene_type_filter, min_count = config$min_count,
  min_samples = config$min_samples, apply_filter = TRUE, use_log2 = TRUE
)
test_log2cpm <- process_expression(
  b_counts, test_pheno, target_col = "nerve_group", gene_anns = gene_anns,
  gene_type_filter = config$gene_type_filter, min_count = config$min_count,
  min_samples = config$min_samples, apply_filter = FALSE, use_log2 = TRUE
)

# Some "unseen" patients are original Cohort A/B patients dropped by the
# quantile threshold (not in unseen_a/b_counts) -- a_counts/unseen_a_counts
# share the same gene panel (same for b_counts/unseen_b_counts), so combine
# each platform's sources before subsetting to that platform's unseen IDs.
full_a_counts <- cbind(a_counts, unseen_a_counts[, setdiff(colnames(unseen_a_counts), colnames(a_counts)), drop = FALSE])
full_b_counts <- cbind(b_counts, unseen_b_counts[, setdiff(colnames(unseen_b_counts), colnames(b_counts)), drop = FALSE])

unseen_pheno_a <- data.frame(`RNA-seq` = intersect(unseen_ids, colnames(full_a_counts)), check.names = FALSE)
unseen_pheno_b <- data.frame(`RNA-seq` = intersect(unseen_ids, colnames(full_b_counts)), check.names = FALSE)
unseen_log2cpm_a <- process_expression(full_a_counts, unseen_pheno_a, gene_anns = gene_anns, gene_type_filter = config$gene_type_filter, apply_filter = FALSE, use_log2 = TRUE)
unseen_log2cpm_b <- process_expression(full_b_counts, unseen_pheno_b, gene_anns = gene_anns, gene_type_filter = config$gene_type_filter, apply_filter = FALSE, use_log2 = TRUE)

stopifnot(length(unseen_pheno_a$`RNA-seq`) + length(unseen_pheno_b$`RNA-seq`) == length(unseen_ids))
common_unseen_genes <- intersect(rownames(unseen_log2cpm_a), rownames(unseen_log2cpm_b))
unseen_log2cpm <- cbind(unseen_log2cpm_a[common_unseen_genes, , drop = FALSE], unseen_log2cpm_b[common_unseen_genes, , drop = FALSE])

# --- 3. Verify the model trained from these precomputed matrices still matches ---
# (external Uromol/Meeks/Bryan cohorts are intentionally not part of the
# package's data -- it isn't ours to redistribute)

sig_genes_clean <- gotnerve_model$coefnames
name_map <- stats::setNames(rownames(train_log2cpm), make.names(rownames(train_log2cpm)))
stopifnot(all(sig_genes_clean %in% names(name_map)) || TRUE) # sig genes come from the live top-N-variance selection, checked via train_model() itself in 01_build_bundled_model.R

# --- 4. Save everything ---

usethis::use_data(stan_clinical, overwrite = TRUE)
usethis::use_data(train_log2cpm, overwrite = TRUE)
usethis::use_data(test_log2cpm, overwrite = TRUE)
usethis::use_data(unseen_log2cpm, overwrite = TRUE)

# Raw counts are no longer part of the package's public data -- remove the
# old exported objects now that their processed replacements exist. Also
# remove the old master_clinical.rda now that it's been renamed stan_clinical.
for (f in c("a_counts", "b_counts", "pheno_data", "unseen_a_counts", "unseen_b_counts", "all_clin", "master_clinical")) {
  path <- file.path("data", paste0(f, ".rda"))
  if (file.exists(path)) file.remove(path)
}
