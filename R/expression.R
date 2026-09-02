#' Restrict gene names to a set of GENE_TYPE annotations
#'
#' Exact match against `gene_anns$GENE_TYPE`: a compound type string like
#' `"lncRNA; protein_coding"` is treated as its own distinct type and will
#' NOT match a filter of `"protein_coding"` alone.
#'
#' @param gene_names Character vector of gene names to restrict.
#' @param gene_anns Data frame with `GENE_NAME` and `GENE_TYPE` columns.
#' @param gene_type_filter Character vector of `GENE_TYPE` values to keep.
#' @return The subset of `gene_names` whose annotated type is in
#'   `gene_type_filter`.
#' @noRd
filter_genes_by_type <- function(gene_names, gene_anns, gene_type_filter) {
  valid_genes <- gene_anns$GENE_NAME[gene_anns$GENE_TYPE %in% gene_type_filter]
  intersect(gene_names, valid_genes)
}

#' Counts to log2-CPM, the training-pipeline version
#'
#' Restricts to a gene-type filter, matches samples against a phenotype
#' table, optionally applies `edgeR::filterByExpr()`, then computes (log2)
#' CPM with
#' `edgeR::cpm(normalized.lib.sizes = FALSE)`. This is the training-time
#' function — it needs `pheno_sub` to know which samples to keep and (when
#' `apply_filter = TRUE`) `target_col` to group by for `filterByExpr()`. For
#' scoring already-collected new data (no cohort/phenotype bookkeeping
#' needed), use [process_gene_expression()] instead.
#'
#' @param counts_mat Genes (rows) x samples (cols) raw counts matrix.
#' @param pheno_sub Data frame with an `` `RNA-seq` `` column identifying
#'   which columns of `counts_mat` to keep.
#' @param target_col Column of `pheno_sub` used as the `filterByExpr()`
#'   grouping factor. Only required when `apply_filter = TRUE`.
#' @param gene_anns Data frame with `GENE_NAME`/`GENE_TYPE` columns.
#' @param gene_type_filter `GENE_TYPE` values to keep before normalization.
#' @param min_count,min_samples Passed to `edgeR::filterByExpr()` as
#'   `min.count` and `min.prop = min_samples / ncol(sub_counts)`.
#' @param apply_filter Whether to apply `edgeR::filterByExpr()`. `TRUE` for
#'   the training cohort, `FALSE` for a held-out cohort that must keep every
#'   gene the training cohort kept.
#' @param use_log2 Return log2-CPM (`edgeR::cpm(log = TRUE, prior.count = 1)`)
#'   instead of linear CPM.
#'
#' @return A genes x samples (log2) CPM matrix.
process_expression <- function(counts_mat, pheno_sub, target_col = NULL, gene_anns,
                               gene_type_filter = c("protein_coding", "miRNA", "lncRNA"),
                               min_count = 8, # slightly lower than edgeR default of 10 to include potentially low expression axonal transcripts
                               min_samples = 5,
                               apply_filter = TRUE, use_log2 = TRUE) {
  valid_genes <- filter_genes_by_type(rownames(counts_mat), gene_anns, gene_type_filter)
  common_pts <- intersect(pheno_sub$`RNA-seq`, colnames(counts_mat))

  sub_counts <- counts_mat[valid_genes, common_pts, drop = FALSE]
  sub_pheno <- pheno_sub[match(common_pts, pheno_sub$`RNA-seq`), , drop = FALSE]

  if (apply_filter) {
    keep <- edgeR::filterByExpr(
      sub_counts,
      group = sub_pheno[[target_col]],
      min.count = min_count,
      min.prop = min_samples / ncol(sub_counts)
    )
    sub_counts <- sub_counts[keep, , drop = FALSE]
  }

  dge <- edgeR::DGEList(sub_counts)
  edgeR::cpm(dge, normalized.lib.sizes = FALSE, log = use_log2, prior.count = 1)
}

#' Process user-supplied expression data of a declared type
#'
#' The general, prediction-time counterpart to [process_expression()]: takes a
#' plain genes x samples matrix (or data frame) plus an explicit declaration
#' of what it already is, and returns log2-CPM without re-deriving anything
#' the caller has already done. No phenotype table, target column, or
#' `filterByExpr()` step is involved (those only make sense for the training
#' cohort's own gene selection).
#'
#' @param expr_data Genes (rows) x samples (cols) matrix or data frame.
#' @param input_type One of `"counts"` (raw counts; gene-type filtered and
#'   converted with `edgeR::DGEList()`/`edgeR::cpm()`), `"cpm"` (already
#'   library-size normalized on a linear scale; only log2-transformed, as
#'   `log2(cpm + 1)`), or `"log2cpm"` (already fully processed; returned
#'   unchanged).
#' @param gene_anns Data frame with `GENE_NAME`/`GENE_TYPE` columns. Only used
#'   when `input_type = "counts"`; ignored (with a warning if genes can't be
#'   filtered) for `"cpm"`/`"log2cpm"` since those are already normalized.
#' @param gene_type_filter `GENE_TYPE` values to keep when `input_type =
#'   "counts"`.
#'
#' @return A genes x samples log2-CPM matrix.
process_gene_expression <- function(expr_data,
                                   input_type = c("counts", "cpm", "log2cpm"),
                                   gene_anns = NULL,
                                   gene_type_filter = c("protein_coding", "miRNA", "lncRNA")) {
  input_type <- match.arg(input_type)
  mat <- as.matrix(expr_data)

  if (input_type == "log2cpm") {
    return(mat)
  }

  if (!is.null(gene_anns)) {
    valid_genes <- filter_genes_by_type(rownames(mat), gene_anns, gene_type_filter)
    mat <- mat[valid_genes, , drop = FALSE]
  }

  if (input_type == "counts") {
    dge <- edgeR::DGEList(mat)
    return(edgeR::cpm(dge, normalized.lib.sizes = FALSE, log = TRUE, prior.count = 1))
  }

  # input_type == "cpm": already library-size normalized on a linear scale.
  log2(mat + 1)
}
