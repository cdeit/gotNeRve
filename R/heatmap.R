#' Signature heatmap: final model genes across train + test patients
#'
#' Row-scaled (z-score) log2-CPM for the genes with a non-zero coefficient
#' in the final model, with a left row annotation of each gene's final
#' glmnet coefficient and top annotations for cohort, clinical covariates,
#' true/predicted class, predicted probability, and total nerve area
#' percent — no other nerve-quantification annotation tracks are shown.
#'
#' @param fit A `gotnerve_fit` from [train_model()].
#' @param order_by Column of `fit$pheno` to sort columns by (ascending),
#'   or `NULL` to hierarchically cluster columns instead. Defaults to
#'   `"nerve_obj_area_percent"`.
#' @param split_cohort If `TRUE`, split columns into separate Cohort A/B
#'   blocks. Defaults to `FALSE` (one continuous block).
#' @return A `ComplexHeatmap::Heatmap` object.
plot_signature_heatmap <- function(fit, order_by = "nerve_obj_area_percent", split_cohort = FALSE) {
  raw_coef <- stats::coef(fit$model$finalModel, s = fit$model$bestTune$lambda)
  all_coefs <- as.vector(raw_coef)
  names(all_coefs) <- rownames(raw_coef)
  # Sign-flip: raw glmnet coefficients target the 2nd factor level ("Low"),
  # so negate to show each gene's association with "High" instead.
  all_coefs <- -all_coefs[fit$model$coefnames]

  sig_genes_clean <- names(all_coefs)[all_coefs != 0]
  name_map <- stats::setNames(rownames(fit$train_cpm), make.names(rownames(fit$train_cpm)))
  sig_genes <- name_map[sig_genes_clean]

  x_mat <- cbind(fit$train_cpm[sig_genes, , drop = FALSE], fit$test_cpm[sig_genes, , drop = FALSE])
  x_mat <- t(scale(t(x_mat)))
  x_mat[is.na(x_mat)] <- 0

  all_ids <- c(fit$train_ids, fit$test_ids)
  clin <- fit$pheno[match(all_ids, fit$pheno$`RNA-seq`), ]

  true_class <- as.character(clin[[fit$target]])
  cv_preds <- fit$cv_model$pred
  pos_class <- levels(fit$train_data[[fit$target]])[1]

  pred_class <- c(as.character(cv_preds$pred), as.character(fit$test_predictions$Pred_Class))
  pred_prob <- c(cv_preds[[pos_class]], fit$test_predictions[[pos_class]])

  top_df <- data.frame(
    Cohort = factor(clin$Cohort),
    True_Class = factor(true_class, levels = c("Low", "High")),
    Predicted_Class = factor(pred_class, levels = c("Low", "High")),
    Prediction_Prob_of_High = pred_prob
  )
  top_cols <- list(
    Cohort = c("A" = "#1D3557", "B" = "#00A9A5"),
    True_Class = c("Low" = "#ffb347", "High" = "#d94801"),
    Predicted_Class = c("Low" = "#ffb347", "High" = "#d94801"),
    Prediction_Prob_of_High = circlize::colorRamp2(
      seq(0, 1, length.out = 9),
      c("#ffffff", "#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000")
    )
  )
  if ("BCG_failure" %in% colnames(clin)) {
    top_df$BCG_failure <- factor(clin$BCG_failure, levels = c("No", "Yes"))
    top_cols$BCG_failure <- c("No" = "gray", "Yes" = "black")
  }
  if ("Progression" %in% colnames(clin)) {
    top_df$Progression <- factor(clin$Progression, levels = c("No Progression", "Progression"))
    top_cols$Progression <- c("No Progression" = "gray", "Progression" = "black")
  }
  if ("Sex" %in% colnames(clin)) {
    top_df$Sex <- factor(clin$Sex)
    top_cols$Sex <- c("Female" = "#A6BE54", "F" = "#A6BE54", "Male" = "#006400", "M" = "#006400")
  }
  top_df <- top_df[, c(intersect(c("Cohort", "BCG_failure", "Progression", "Sex"), colnames(top_df)), "True_Class", "Predicted_Class", "Prediction_Prob_of_High")]

  nerve_col <- "nerve_obj_area_percent"
  top_ann_args <- c(as.list(top_df), list(col = top_cols, gap = grid::unit(0.5, "mm")))
  if (nerve_col %in% colnames(clin)) {
    top_ann_args$`STaN+\narea (%)` <- ComplexHeatmap::anno_barplot(clin[[nerve_col]], height = grid::unit(1.5, "cm"), bar_width = 0.9, gp = grid::gpar(fill = "#fdb462", col = "white"))
  }
  top_ann <- do.call(ComplexHeatmap::HeatmapAnnotation, top_ann_args)

  gene_coefs <- all_coefs[sig_genes_clean]
  names(gene_coefs) <- sig_genes
  abs_max <- max(abs(gene_coefs))
  coef_col_fun <- circlize::colorRamp2(c(-abs_max, 0, abs_max), c("#4A7BE8", "#F8F5F0", "#E65A5A"))
  row_ann <- ComplexHeatmap::rowAnnotation(
    Final_Coef = gene_coefs[rownames(x_mat)],
    col = list(Final_Coef = coef_col_fun),
    annotation_name_rot = 0, annotation_name_side = "top"
  )

  cluster_cols <- is.null(order_by) || !(order_by %in% colnames(clin))
  col_order <- if (!cluster_cols) order(clin[[order_by]], na.last = TRUE) else NULL
  column_split <- if (isTRUE(split_cohort)) factor(clin$Cohort) else NULL

  z_pal <- c("#4A2FA8", "#7A62E8", "#C4B9FF", "#FAF5EE", "#F4B0B0", "#E04B4B", "#A1122A")
  z_breaks <- seq(-max(abs(x_mat)), max(abs(x_mat)), length.out = length(z_pal))

  ComplexHeatmap::Heatmap(
    x_mat,
    col = circlize::colorRamp2(z_breaks, z_pal),
    name = "Gene expression z-score",
    top_annotation = top_ann,
    left_annotation = row_ann,
    cluster_columns = cluster_cols,
    cluster_column_slices = FALSE,
    column_order = if (!cluster_cols) col_order else NULL,
    cluster_rows = TRUE,
    column_split = column_split,
    show_column_names = FALSE,
    row_names_gp = grid::gpar(fontsize = 10),
    row_title = "Genes in the Final Predictive Model",
    row_title_side = "left",
    column_title = "Patients in Train and Test Cohorts"
  )
}
