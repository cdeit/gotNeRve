#' Impute missing signature genes using training-cohort medians
#'
#' Used by [gotNeRve()] when `new_data` is missing one of the model's
#' signature genes: fills it in with that gene's median log2-CPM value
#' across the training cohort, rather than assuming zero expression (which,
#' after centering/scaling on the training mean, would look like unusually
#' low expression instead of a neutral, typical value).
#'
#' @param expr Genes (rows) x samples (cols) log2-CPM matrix, gene names
#'   already run through `make.names()`.
#' @param sig_genes `make.names()`-ed signature gene names required by the
#'   model.
#' @param reference Genes x samples log2-CPM matrix (original, non-`make.names()`
#'   gene names) used to compute the per-gene medians. Defaults to the
#'   package's bundled `train_log2cpm`.
#'
#' @return `expr` subset and reordered to exactly `sig_genes`, with any
#'   missing rows filled in with that gene's median value in `reference`.
impute_missing_genes <- function(expr, sig_genes, reference = NULL) {
  if (is.null(reference)) reference <- gotNeRve::train_log2cpm

  missing_genes <- setdiff(sig_genes, rownames(expr))
  if (length(missing_genes) > 0) {
    ref_name_map <- stats::setNames(rownames(reference), make.names(rownames(reference)))
    orig_names <- ref_name_map[missing_genes]
    if (any(is.na(orig_names))) {
      stop(
        "Missing signature gene(s) not found in the reference dataset: ",
        paste(missing_genes[is.na(orig_names)], collapse = ", ")
      )
    }
    medians <- apply(reference[orig_names, , drop = FALSE], 1, stats::median, na.rm = TRUE)
    filler <- matrix(
      medians,
      nrow = length(missing_genes), ncol = ncol(expr),
      dimnames = list(missing_genes, colnames(expr))
    )
    expr <- rbind(expr, filler)
  }
  expr[sig_genes, , drop = FALSE]
}

#' Score new expression data against a fitted model (internal)
#'
#' The shared scoring logic behind [gotNeRve()] and the package's own
#' training pipeline (which scores the held-out test cohort against a
#' freshly fit model, not the bundled one — see [train_model()]). Aligns
#' `x` to the model's signature genes (`model$coefnames`), imputing any
#' that are missing with that gene's median log2-CPM value in the training
#' cohort (see [impute_missing_genes()]), and returns the predicted class
#' and probability. Kept separate from [gotNeRve()] so that function can
#' stay a plain "score against the bundled model" entry point with no
#' `model` argument, while training/evaluation code can still score
#' against an arbitrary in-progress model.
#'
#' @param x Genes (rows) x samples (cols) matrix or data frame.
#' @param model A fitted `caret::train` model.
#' @param input_type One of `"counts"`, `"cpm"`, or `"log2cpm"` — see
#'   [process_gene_expression()].
#' @param gene_anns Gene annotation data frame with `GENE_NAME`/`GENE_TYPE`.
#'   Only used when `input_type = "counts"`; defaults to the package's
#'   bundled `gene_anns`.
#' @param gene_type_filter `GENE_TYPE` values to keep when `input_type =
#'   "counts"`.
#'
#' @return A data frame with one row per sample: `` `RNA-seq` `` (from
#'   `colnames(x)`), `Pred_Class`, and one probability column per outcome
#'   level (e.g. `High`, `Low`).
#' @noRd
score_nerve_class <- function(x, model,
                               input_type = c("counts", "cpm", "log2cpm"),
                               gene_anns = NULL,
                               gene_type_filter = c("protein_coding", "miRNA", "lncRNA")) {
  input_type <- match.arg(input_type)
  if (is.null(gene_anns) && input_type == "counts") gene_anns <- gotNeRve::gene_anns

  sample_ids <- colnames(x)
  if (length(sample_ids) == 0) {
    return(data.frame(`RNA-seq` = character(0), Pred_Class = factor(character(0), levels = model$obsLevels), check.names = FALSE))
  }

  expr <- process_gene_expression(x, input_type, gene_anns, gene_type_filter)
  rownames(expr) <- make.names(rownames(expr))

  sig_genes <- model$coefnames
  expr <- impute_missing_genes(expr, sig_genes)

  newdata <- as.data.frame(t(expr))
  colnames(newdata) <- make.names(colnames(newdata))

  pred_class <- stats::predict(model, newdata = newdata, type = "raw")
  pred_prob <- stats::predict(model, newdata = newdata, type = "prob")

  cbind(
    data.frame(`RNA-seq` = sample_ids, Pred_Class = pred_class, check.names = FALSE),
    pred_prob
  )
}

#' Predict STaN class from gene expression using the bundled gotNeRve model
#'
#' The package's single prediction entry point for new data: loads the
#' bundled, locked model automatically and scores `x` against it — no
#' setup or configuration required. Aligns `x` to the model's signature
#' genes, imputing any that are missing with that gene's median log2-CPM
#' value in the training cohort (see [impute_missing_genes()]).
#'
#' @param x Genes (rows) x samples (cols) matrix or data frame.
#' @param input_type One of `"counts"`, `"cpm"`, or `"log2cpm"` — see
#'   [process_gene_expression()]. Use `"log2cpm"` if `x` is already
#'   fully processed and should be used as-is.
#'
#' @return A data frame with one row per sample: `sample_id` (from
#'   `colnames(x)`), `predicted_STaN_class`, and `probability_<level>` for
#'   each outcome level (e.g. `probability_high`, `probability_low`).
#' @export
gotNeRve <- function(x, input_type = c("counts", "cpm", "log2cpm")) {
  input_type <- match.arg(input_type)
  result <- score_nerve_class(
    x, gotnerve_model, input_type,
    gene_anns = gotNeRve::gene_anns,
    gene_type_filter = c("protein_coding", "miRNA", "lncRNA")
  )

  prob_cols <- setdiff(colnames(result), c("RNA-seq", "Pred_Class"))
  colnames(result)[colnames(result) %in% prob_cols] <- paste0("probability_", tolower(prob_cols))

  names(result)[names(result) == "RNA-seq"] <- "sample_id"
  names(result)[names(result) == "Pred_Class"] <- "predicted_STaN_class"
  result
}
