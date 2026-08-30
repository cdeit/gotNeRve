#' Re-label a `nerve_*` target into extreme-quantile High/Low groups
#'
#' For a target column starting with `"nerve_group"`, looks up its corresponding
#' continuous column (`"nerve_group_x"` -> `"nerve_obj_area_percent_x"`),
#' computes the `q_val`/`1 - q_val` quantiles *within each cohort
#' independently*, relabels patients in the extremes as `"Low"`/`"High"`, and
#' drops everyone in between. Non-`nerve_group*` targets, or a `NULL` q_val,
#' pass through unchanged.
#'
#' @param pheno Phenotype data frame.
#' @param target Target column name.
#' @param q_val Quantile cutoff (e.g. `0.33` for tertiles), or `NULL` to skip.
#' @return `pheno`, re-labeled and with intermediate patients removed (for
#'   applicable targets), otherwise unchanged.
apply_quantile_threshold <- function(pheno, target, q_val) {
  if (is.null(q_val) || !grepl("^nerve_group", target)) {
    return(pheno)
  }
  cont_col <- gsub("^nerve_group", "nerve_obj_area_percent", target)
  if (!cont_col %in% colnames(pheno)) {
    return(pheno)
  }

  if ("Cohort" %in% colnames(pheno)) {
    cohorts <- unique(pheno$Cohort)
  } else {
    cohorts <- "All"
    pheno$Cohort_Temp <- "All"
  }

  res_list <- lapply(cohorts, function(cohort) {
    p_sub <- if ("Cohort" %in% colnames(pheno)) pheno[pheno$Cohort == cohort, ] else pheno[pheno$Cohort_Temp == cohort, ]
    quantiles <- stats::quantile(p_sub[[cont_col]], probs = c(q_val, 1 - q_val), na.rm = TRUE)

    new_class <- rep(NA, nrow(p_sub))
    new_class[p_sub[[cont_col]] <= quantiles[1]] <- "Low"
    new_class[p_sub[[cont_col]] >= quantiles[2]] <- "High"

    p_sub[[target]] <- factor(new_class, levels = c("Low", "High"))
    p_sub[!is.na(p_sub[[target]]), ]
  })

  res <- do.call(rbind, res_list)
  if ("Cohort_Temp" %in% colnames(res)) res$Cohort_Temp <- NULL
  res
}

#' Select the highest-variance genes in a training-cohort expression matrix
#'
#' @param train_expr Genes x samples (log2-CPM) matrix for the training
#'   cohort.
#' @param test_expr Genes x samples matrix for the held-out cohort; only
#'   genes present in both are returned.
#' @param top_n_var Number of top-variance training-cohort genes to consider.
#' @return Character vector of gene names, ordered by decreasing
#'   training-cohort variance.
top_variance_genes <- function(train_expr, test_expr, top_n_var = 1000) {
  gene_vars <- apply(train_expr, 1, stats::var, na.rm = TRUE)
  top_train_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(top_n_var, length(gene_vars)))]
  intersect(top_train_genes, rownames(test_expr))
}
