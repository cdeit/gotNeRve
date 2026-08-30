#' ROC curve for the training (out-of-fold) or test cohort
#'
#' Train uses the glmnet cross-validated out-of-fold predictions
#' (`fit$cv_model$pred`). Test uses `fit$test_predictions`, obtained via
#' [gotNeRve()] on the raw test-cohort counts.
#'
#' @param fit A `gotnerve_fit` from [train_model()].
#' @param dataset `"train"` or `"test"`.
#' @return A ggplot object.
plot_roc <- function(fit, dataset = c("train", "test")) {
  dataset <- match.arg(dataset)

  if (dataset == "train") {
    cv_preds <- fit$cv_model$pred
    class_lvls <- levels(cv_preds$obs)
    actual <- cv_preds$obs
    probs_pos <- cv_preds[[class_lvls[1]]]
    title_suffix <- paste0("Cohort ", fit$train_cohort, " (Training Out-of-Fold CV)")
  } else {
    actual <- fit$test_data[[fit$target]]
    class_lvls <- levels(actual)
    probs_pos <- fit$test_predictions[[class_lvls[1]]]
    title_suffix <- paste0("Cohort ", fit$test_cohort, " (Test)")
  }

  roc_obj <- pROC::roc(actual, probs_pos, levels = rev(class_lvls), direction = "<", quiet = TRUE)
  auc_val <- pROC::auc(roc_obj)
  n_pos <- sum(actual == class_lvls[1])
  n_neg <- sum(actual == class_lvls[2])

  pROC::ggroc(roc_obj, colour = "#e74c3c", linewidth = 1.2, legacy.axes = TRUE) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    ggplot2::labs(
      title = paste("ROC Curve -", title_suffix),
      subtitle = sprintf(
        "AUC = %.3f | Target: %s\nn = %d (%d %s, %d %s)",
        auc_val, class_lvls[1], n_pos + n_neg, n_pos, class_lvls[1], n_neg, class_lvls[2]
      ),
      x = "1 - Specificity", y = "Sensitivity"
    ) +
    cowplot::theme_cowplot()
}

#' Confusion matrix for the training (out-of-fold) or test cohort
#'
#' @inheritParams plot_roc
#' @return A ggplot object.
plot_confusion_matrix <- function(fit, dataset = c("train", "test")) {
  dataset <- match.arg(dataset)

  if (dataset == "train") {
    cv_preds <- fit$cv_model$pred
    actual <- cv_preds$obs
    preds <- cv_preds$pred
    title_suffix <- paste0("Cohort ", fit$train_cohort, " (Training Out-of-Fold CV)")
  } else {
    actual <- fit$test_data[[fit$target]]
    preds <- fit$test_predictions$Pred_Class
    title_suffix <- paste0("Cohort ", fit$test_cohort, " (Test)")
  }

  cm <- caret::confusionMatrix(preds, actual)
  acc <- cm$overall[["Accuracy"]]

  cm_df <- as.data.frame(cm$table)
  cm_df$Prediction <- factor(cm_df$Prediction, levels = rev(levels(cm_df$Prediction)))
  cm_df$Reference <- factor(cm_df$Reference, levels = rev(levels(cm_df$Reference)))

  ggplot2::ggplot(cm_df, ggplot2::aes(x = .data$Reference, y = .data$Prediction, fill = .data$Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$Freq), vjust = 1, size = 6) +
    ggplot2::scale_fill_gradient(low = "#ffffff", high = "#2c3e50") +
    ggplot2::labs(
      title = paste0("Confusion Matrix - ", title_suffix, "\nAcc: ", round(acc, 3)),
      x = "True Class", y = "Predicted Class"
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "none")
}

#' Predicted-probability boxplot by true class
#'
#' @inheritParams plot_roc
#' @return A ggplot object.
plot_risk_distribution <- function(fit, dataset = c("train", "test")) {
  dataset <- match.arg(dataset)

  if (dataset == "train") {
    cv_preds <- fit$cv_model$pred
    class_lvls <- levels(cv_preds$obs)
    actual <- cv_preds$obs
    probs_pos <- cv_preds[[class_lvls[1]]]
    title_suffix <- paste0("Cohort ", fit$train_cohort, " (Training Out-of-Fold CV)")
  } else {
    actual <- fit$test_data[[fit$target]]
    class_lvls <- levels(actual)
    probs_pos <- fit$test_predictions[[class_lvls[1]]]
    title_suffix <- paste0("Cohort ", fit$test_cohort, " (Test)")
  }

  df <- data.frame(Actual = actual, Prob = probs_pos)
  pal <- c("Low" = "#ffb347", "High" = "#d94801")
  for (l in levels(actual)) if (!l %in% names(pal)) pal[[l]] <- "gray50"

  ggplot2::ggplot(df, ggplot2::aes(x = factor(.data$Actual, levels = rev(levels(.data$Actual))), y = .data$Prob, fill = .data$Actual)) +
    ggplot2::geom_boxplot(alpha = 0.7, outliers = FALSE) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::ylim(0, 1) +
    ggpubr::stat_compare_means(method = "wilcox.test") +
    ggplot2::labs(title = paste("Predicted Risk Distribution -", title_suffix), y = paste("Probability of", class_lvls[1]), x = "True Class") +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "none")
}

#' Overlayed ROC: training CV-out-of-fold vs. held-out test cohort
#'
#' Sensitivity/specificity/accuracy are always reported at a fixed 0.5
#' probability cutoff (no optimized-threshold option).
#'
#' @param fit A `gotnerve_fit` from [train_model()].
#' @return A ggplot object.
plot_roc_overlay <- function(fit) {
  cv_preds <- fit$cv_model$pred
  class_lvls <- levels(cv_preds$obs)

  train_roc <- pROC::roc(cv_preds$obs, cv_preds[[class_lvls[1]]], levels = rev(class_lvls), direction = "<", quiet = TRUE)
  test_actual <- fit$test_data[[fit$target]]
  test_roc <- pROC::roc(test_actual, fit$test_predictions[[class_lvls[1]]], levels = rev(class_lvls), direction = "<", quiet = TRUE)

  # Fixed 0.5 probability cutoff only -- no Youden's-J-optimized threshold.
  train_coords <- as.list(suppressWarnings(pROC::coords(train_roc, x = 0.5, input = "threshold", ret = c("sensitivity", "specificity", "accuracy")))[1, ])
  test_coords <- as.list(suppressWarnings(pROC::coords(test_roc, x = 0.5, input = "threshold", ret = c("sensitivity", "specificity", "accuracy")))[1, ])

  roc_df <- rbind(
    data.frame(fpr = 1 - test_roc$specificities, tpr = test_roc$sensitivities, Cohort = "Test"),
    data.frame(fpr = 1 - train_roc$specificities, tpr = train_roc$sensitivities, Cohort = "CV OOF (Train)")
  )

  ggplot2::ggplot(roc_df, ggplot2::aes(x = .data$fpr, y = .data$tpr, color = .data$Cohort, linetype = .data$Cohort, alpha = .data$Cohort)) +
    ggplot2::geom_path(linewidth = 1.2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    ggplot2::scale_linetype_manual(values = c("Test" = "solid", "CV OOF (Train)" = "dashed")) +
    ggplot2::scale_color_manual(values = c("Test" = "#00B8A9", "CV OOF (Train)" = "#1D3557")) +
    ggplot2::scale_alpha_manual(values = c("Test" = 1.0, "CV OOF (Train)" = 0.75)) +
    ggplot2::labs(
      title = "Model Generalization: Train CV OOF vs. Test Cohort",
      subtitle = sprintf(
        "CV OOF - AUC: %.3f | Sens: %.3f | Spec: %.3f | Acc: %.3f\nTest - AUC: %.3f | Sens: %.3f | Spec: %.3f | Acc: %.3f",
        pROC::auc(train_roc), train_coords$sensitivity, train_coords$specificity, train_coords$accuracy,
        pROC::auc(test_roc), test_coords$sensitivity, test_coords$specificity, test_coords$accuracy
      ),
      x = "1 - Specificity (False Positive Rate)", y = "Sensitivity (True Positive Rate)"
    ) +
    cowplot::theme_cowplot()
}

#' GLMnet coefficient shrinkage path for the final model's active genes
#'
#' @param fit A `gotnerve_fit` from [train_model()].
#' @param label_n Number of genes (ranked by |coefficient| at the fitted
#'   lambda) to label directly on the plot.
#' @return A ggplot object.
plot_coefficient_path <- function(fit, label_n = 14) {
  fm <- fit$model$finalModel
  # Negated: raw glmnet coefficients target the 2nd factor level ("Low"),
  # so this shows each gene's association with "High" instead, consistent
  # with the heatmap/boxplot coefficient convention.
  beta <- -as.matrix(fm$beta)
  lambdas <- fm$lambda

  active_genes <- rownames(beta)[rowSums(beta != 0) > 0]
  beta_active <- beta[active_genes, , drop = FALSE]

  df_lines <- expand.grid(Gene = active_genes, LambdaIndex = seq_along(lambdas), stringsAsFactors = FALSE)
  df_lines$LogLambda <- log(lambdas[df_lines$LambdaIndex])
  df_lines$Estimate <- as.vector(beta_active)

  min_ll <- min(df_lines$LogLambda)
  final_coefs <- beta_active[, length(lambdas)]
  top_genes <- head(names(sort(abs(final_coefs), decreasing = TRUE)), label_n)
  df_labels <- data.frame(Gene = top_genes, LogLambda = min_ll, Estimate = final_coefs[top_genes], stringsAsFactors = FALSE)

  ggplot2::ggplot(df_lines, ggplot2::aes(x = .data$LogLambda, y = .data$Estimate, group = .data$Gene, color = .data$Gene)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.8) +
    ggplot2::geom_line(alpha = 0.8, linewidth = 0.8) +
    ggrepel::geom_text_repel(
      data = df_labels, ggplot2::aes(label = .data$Gene, x = .data$LogLambda, y = .data$Estimate),
      hjust = 1, direction = "y", nudge_x = -0.2, segment.color = "gray70",
      size = 4, fontface = "bold", show.legend = FALSE
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(20, 30, 20, 10), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = "GLMnet Coefficient Shrinkage Path", x = "Log Lambda", y = "Coefficients") +
    ggplot2::coord_cartesian(clip = "off")
}

#' Top 50 bootstrap-stable features
#'
#' Bar chart of the bootstrapped feature-stability filter's inclusion
#' frequencies, with a dashed line at the inclusion-frequency threshold used
#' to decide which genes survived into the tuning sweep.
#'
#' @param fit A `gotnerve_fit` from [train_model()] (must have been trained
#'   with `glm_stability = TRUE`).
#' @param threshold Inclusion-frequency cutoff to draw as a vertical line
#'   (defaults to `fit`'s own `gotnerve_config()` value if available, else
#'   `0.25`).
#' @return A ggplot object.
plot_stability_frequencies <- function(fit, threshold = 0.25) {
  df <- fit$boot_summary[order(-fit$boot_summary$Inclusion_Freq), ][seq_len(min(50, nrow(fit$boot_summary))), ]

  ggplot2::ggplot(df, ggplot2::aes(x = .data$Inclusion_Freq, y = stats::reorder(.data$Gene, .data$Inclusion_Freq))) +
    ggplot2::geom_bar(stat = "identity", fill = "#3498db", color = "black") +
    ggplot2::geom_vline(xintercept = threshold, color = "#e74c3c", linetype = "dashed", linewidth = 1.2) +
    cowplot::theme_cowplot() +
    ggplot2::labs(
      title = "Top 50 bootstrap-stable features",
      subtitle = "Dashed red line indicates threshold cutoff for model inclusion",
      x = "Bootstrap inclusion frequency", y = "Feature"
    ) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9, face = "bold"))
}

#' Classic `cv.glmnet` cross-validation error curve at the fitted alpha
#'
#' Refits `glmnet::cv.glmnet()` at the model's tuned alpha purely for this
#' base-R diagnostic plot (`set.seed(42)` for reproducibility), independent
#' of the tuning sweep used to select alpha/lambda originally.
#'
#' @param fit A `gotnerve_fit` from [train_model()].
#' @return Invisibly, the `cv.glmnet` object; called for its plot side effect.
plot_cv_error <- function(fit) {
  y <- fit$train_data[[fit$target]]
  x <- as.matrix(fit$train_data[, setdiff(colnames(fit$train_data), fit$target), drop = FALSE])
  set.seed(42)
  cv_fit <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = fit$model$bestTune$alpha)
  graphics::plot(cv_fit)
  graphics::title(paste("cv.glmnet Error Curve (Fixed Alpha =", fit$model$bestTune$alpha, ")"), line = 3)
  invisible(cv_fit)
}

#' Boxplots of the final model's signature genes by true class
#'
#' @param fit A `gotnerve_fit` from [train_model()].
#' @param n_genes Number of genes (ranked by |final coefficient|) to plot, or
#'   `NULL` for all genes with a non-zero final coefficient.
#' @return A ggplot object (faceted by gene).
plot_gene_boxplots <- function(fit, n_genes = NULL) {
  final_coefs <- as.vector(stats::coef(fit$model$finalModel, s = fit$model$bestTune$lambda))
  names(final_coefs) <- rownames(stats::coef(fit$model$finalModel, s = fit$model$bestTune$lambda))
  final_coefs <- final_coefs[setdiff(names(final_coefs), "(Intercept)")]
  final_coefs <- final_coefs[final_coefs != 0]
  ordered_genes <- names(sort(abs(final_coefs), decreasing = TRUE))
  top_genes <- if (is.null(n_genes)) ordered_genes else head(ordered_genes, n_genes)

  name_map <- stats::setNames(rownames(fit$train_cpm), make.names(rownames(fit$train_cpm)))
  orig_genes <- name_map[top_genes]

  all_cpm <- cbind(fit$train_cpm[orig_genes, , drop = FALSE], fit$test_cpm[orig_genes, , drop = FALSE])
  pheno_sub <- fit$pheno[, c("RNA-seq", fit$target)]
  colnames(pheno_sub) <- c("Sample", "Group")

  plot_df <- as.data.frame(all_cpm)
  plot_df$Gene <- rownames(plot_df)
  plot_df <- tidyr::pivot_longer(plot_df, cols = -"Gene", names_to = "Sample", values_to = "Expression")
  plot_df <- dplyr::left_join(plot_df, pheno_sub, by = "Sample")
  plot_df$Gene <- factor(plot_df$Gene, levels = orig_genes)
  plot_df$Group <- factor(plot_df$Group, levels = c("Low", "High"))

  pal <- c("Low" = "#ffb347", "High" = "#d94801")

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Group, y = .data$Expression, fill = .data$Group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.6, size = 1.2, color = "black") +
    ggplot2::facet_wrap(~Gene, scales = "free_y", ncol = 4) +
    ggplot2::scale_fill_manual(values = pal) +
    ggpubr::stat_compare_means(method = "wilcox.test", size = 3.5) +
    cowplot::theme_cowplot() +
    ggplot2::labs(x = fit$target, y = "log2(CPM + 1)") +
    ggplot2::theme(legend.position = "bottom")
}
