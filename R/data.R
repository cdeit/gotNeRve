#' Gene annotations
#'
#' @format A data frame with `GENE_NAME`, `GENE_TYPE`, and further annotation
#'   columns (description, GO terms, synonyms), one row per gene.
"gene_anns"

#' Clinical table: train/test/unseen dataset membership and outcomes
#'
#' One row per patient actually used somewhere in the pipeline (patients
#' used nowhere are already excluded). `dataset` is `"train"`, `"test"`, or
#' `"unseen"`. `nerve_group`/`nerve_obj_area_percent` are only populated for
#' `"train"`/`"unseen"` datasets; `HG_recur`/`Progression` and their matching
#' `Time_to_*_or_FUend` columns are populated wherever known.
#'
#' @format A data frame with columns `` `RNA-seq` ``, `Cohort`, `dataset`,
#'   `Progression`, `HG_recur`, `Time_to_prog_or_FUend`,
#'   `Time_to_HG_recur_or_FUend`, `Sex`, `BCG_failure`, `nerve_group`,
#'   `nerve_obj_area_percent`.
"stan_clinical"

#' Training-cohort log2-CPM expression
#'
#' Already gene-type filtered and `edgeR::filterByExpr()`-filtered
#' (`min_count = 8`, `min_samples = 5`), then CPM-normalized and
#' log2-transformed. Columns match the `` `RNA-seq` `` IDs with
#' `dataset == "train"` in `stan_clinical`.
#'
#' @format A genes x samples matrix.
"train_log2cpm"

#' Test-cohort log2-CPM expression
#'
#' Gene-type filtered and CPM-normalized/log2-transformed (no
#' `filterByExpr()` filtering — the test cohort keeps every gene the
#' training cohort's filter selected from). Columns match the
#' `` `RNA-seq` `` IDs with `dataset == "test"` in `stan_clinical`.
#'
#' @format A genes x samples matrix.
"test_log2cpm"

#' Internal-unseen validation-cohort log2-CPM expression
#'
#' Gene-type filtered and CPM-normalized/log2-transformed. Covers exactly
#' the ~200 patients in the raw unseen counts files. A handful of these
#' patients are also part of the training/test cohort
#' (`stan_clinical$dataset` is `"train"`/`"test"` for them) and are still
#' included here, matching the pipeline as originally run.
#'
#' @format A genes x samples matrix.
"unseen_log2cpm"
