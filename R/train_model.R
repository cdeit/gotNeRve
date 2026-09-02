#' Train the default gotNeRve elastic-net signature model
#'
#' Cohort split -> gene-type filter + CPM (via [process_expression()]) ->
#' top-variance gene selection -> bootstrapped feature-stability filter ->
#' alpha/lambda tuning sweep -> final fixed-tune refit. Only the elastic-net
#' classifier path is implemented.
#'
#' @section RNG fidelity:
#' To reproduce the exact default model, the
#' random-draw order below matters precisely, including one
#' call that looks dead but isn't: a `foldid` vector is drawn via `sample()`
#' before the stability filter's pre-CV `glmnet::cv.glmnet()` call but is
#' *not* passed to it — it still consumes the RNG stream, so removing it
#' would change every downstream random draw.
#'
#' @param pheno_data Phenotype data frame (e.g. the bundled `pheno_data`).
#' @param a_counts,b_counts Cohort A / Cohort B genes x samples raw counts
#'   matrices.
#' @param gene_anns Gene annotation data frame with `GENE_NAME`/`GENE_TYPE`.
#' @param config A config list from [gotnerve_config()].
#'
#' @return A `gotnerve_fit` object (a list) with elements: `model` (the final
#'   `caret::train` glmnet fit), `cv_model` (the out-of-fold tuning
#'   predictions/results used for train-side ROC and KM plots),
#'   `train_data`/`test_data` (model-ready data frames), `train_ids`/
#'   `test_ids`, `target` (cleaned target name), `train_cohort`/`test_cohort`,
#'   `boot_summary` (stability-filter bootstrap audit trail), `train_cpm`/
#'   `test_cpm` (log2-CPM matrices restricted to the top-variance gene set,
#'   original gene names), `pheno` (combined, filtered train+test phenotype
#'   data), and `test_predictions` (the held-out cohort's predicted class/
#'   probability, obtained via the same scoring logic [gotNeRve()] uses
#'   internally, for full consistency with how any other new cohort is
#'   scored).
train_model <- function(pheno_data, a_counts, b_counts, gene_anns,
                        config = gotnerve_config()) {
  target <- config$target_var

  pheno_filtered <- apply_quantile_threshold(pheno_data, target, config$quantile_threshold)
  pheno_filtered <- pheno_filtered[!is.na(pheno_filtered[[target]]), ]

  train_cohort <- config$train_cohort
  test_cohort <- if (train_cohort == "A") "B" else "A"

  train_pheno <- pheno_filtered[pheno_filtered$Cohort == train_cohort, ]
  test_pheno <- pheno_filtered[pheno_filtered$Cohort == test_cohort, ]

  train_counts_raw <- if (train_cohort == "A") a_counts else b_counts
  test_counts_raw <- if (train_cohort == "A") b_counts else a_counts

  train_cpm_mat <- process_expression(
    train_counts_raw, train_pheno, target_col = target, gene_anns = gene_anns,
    gene_type_filter = config$gene_type_filter,
    min_count = config$min_count, min_samples = config$min_samples,
    apply_filter = TRUE, use_log2 = TRUE
  )
  test_cpm_mat <- process_expression(
    test_counts_raw, test_pheno, target_col = target, gene_anns = gene_anns,
    gene_type_filter = config$gene_type_filter,
    min_count = config$min_count, min_samples = config$min_samples,
    apply_filter = FALSE, use_log2 = TRUE
  )

  top_genes <- top_variance_genes(train_cpm_mat, test_cpm_mat, config$top_n_var)

  train_cpm_final <- as.data.frame(t(train_cpm_mat[top_genes, , drop = FALSE]))
  train_cpm_final$`RNA-seq` <- rownames(train_cpm_final)
  test_cpm_final <- as.data.frame(t(test_cpm_mat[top_genes, , drop = FALSE]))
  test_cpm_final$`RNA-seq` <- rownames(test_cpm_final)

  train_p <- train_pheno[, c("RNA-seq", target)]
  test_p <- test_pheno[, c("RNA-seq", target)]

  train_data <- dplyr::inner_join(train_p, train_cpm_final, by = "RNA-seq")
  test_data <- dplyr::inner_join(test_p, test_cpm_final, by = "RNA-seq")

  train_ids <- train_data$`RNA-seq`
  test_ids <- test_data$`RNA-seq`
  train_data$`RNA-seq` <- NULL
  test_data$`RNA-seq` <- NULL

  train_vec <- as.factor(make.names(train_data[[target]]))
  test_vec <- as.factor(make.names(test_data[[target]]))

  lvls <- levels(train_vec)
  pos_match <- grep("yes|progression|high|positive|case", lvls, ignore.case = TRUE, value = TRUE)
  if (length(pos_match) == 1 && length(lvls) == 2) {
    train_vec <- stats::relevel(train_vec, ref = pos_match[1])
    test_vec <- stats::relevel(test_vec, ref = pos_match[1])
  }
  train_data[[target]] <- train_vec
  test_data[[target]] <- test_vec

  clean_names <- make.names(colnames(train_data))
  colnames(train_data) <- clean_names
  colnames(test_data) <- clean_names
  target_clean <- make.names(target)

  fam <- "binomial"
  boot_summary_res <- NULL

  set.seed(42) # reproducibility seed for the stability filter and tuning sweep

  if (isTRUE(config$glm_stability)) {
    x_mat <- as.matrix(train_data[, setdiff(colnames(train_data), target_clean), drop = FALSE])
    y_vec <- train_data[[target_clean]]

    nfolds <- max(3, config$cv_folds)
    foldid <- sample(rep(seq_len(nfolds), length.out = nrow(x_mat))) # consumed by RNG only; intentionally not passed to cv.glmnet below
    pre_cv <- suppressWarnings(glmnet::cv.glmnet(
      x_mat, y_vec,
      family = fam, alpha = config$glm_stability_alpha, standardize = TRUE, nfolds = nfolds
    ))
    stab_rule <- if (identical(config$glm_stability_rule, "oneSE")) "lambda.1se" else "lambda.min"
    stable_lambda <- pre_cv[[stab_rule]]

    n_iters <- config$glm_bootstraps
    kept_genes <- vector("list", n_iters)
    pos_counts <- stats::setNames(numeric(ncol(x_mat)), colnames(x_mat))
    neg_counts <- stats::setNames(numeric(ncol(x_mat)), colnames(x_mat))

    for (b in seq_len(n_iters)) {
      boot_idx <- sample(seq_len(nrow(x_mat)), replace = TRUE)
      boot_fit <- suppressWarnings(glmnet::glmnet(
        x_mat[boot_idx, , drop = FALSE], y_vec[boot_idx],
        family = fam, alpha = config$glm_stability_alpha, standardize = TRUE, lambda = stable_lambda
      ))
      coefs <- stats::coef(boot_fit, s = stable_lambda)
      c_vec <- as.vector(coefs)
      names(c_vec) <- rownames(coefs)
      c_vec <- c_vec * -1 # caret targets the 1st factor level; glmnet targets the 2nd

      pos_vars <- setdiff(names(c_vec)[c_vec > 0], "(Intercept)")
      neg_vars <- setdiff(names(c_vec)[c_vec < 0], "(Intercept)")
      if (length(pos_vars) > 0) pos_counts[pos_vars] <- pos_counts[pos_vars] + 1
      if (length(neg_vars) > 0) neg_counts[neg_vars] <- neg_counts[neg_vars] + 1

      kept_genes[[b]] <- setdiff(names(c_vec)[c_vec != 0], "(Intercept)")
    }

    gene_freq <- table(unlist(kept_genes)) / n_iters

    tot_nonzero <- pos_counts + neg_counts
    valid_idx <- which(tot_nonzero > 0)
    pct_pos <- stats::setNames(rep(0, length(tot_nonzero)), names(tot_nonzero))
    pct_neg <- stats::setNames(rep(0, length(tot_nonzero)), names(tot_nonzero))
    pct_pos[valid_idx] <- pos_counts[valid_idx] / tot_nonzero[valid_idx]
    pct_neg[valid_idx] <- neg_counts[valid_idx] / tot_nonzero[valid_idx]
    max_sign_cons <- pmax(pct_pos, pct_neg)

    boot_summary_res <- data.frame(
      Gene = names(tot_nonzero[valid_idx]),
      Inclusion_Freq = as.numeric(tot_nonzero[valid_idx] / n_iters),
      Pos_Coef_Freq = round(pct_pos[valid_idx] * 100, 2),
      Neg_Coef_Freq = round(pct_neg[valid_idx] * 100, 2),
      Sign_Consistency = round(max_sign_cons[valid_idx] * 100, 2),
      stringsAsFactors = FALSE
    )

    if (config$glm_stability_sign_threshold > 0) {
      valid_sign_genes <- names(max_sign_cons)[max_sign_cons >= config$glm_stability_sign_threshold]
      stable_genes <- intersect(names(gene_freq)[gene_freq >= config$glm_stability_threshold], valid_sign_genes)
    } else {
      stable_genes <- names(gene_freq)[gene_freq >= config$glm_stability_threshold]
    }

    if (length(stable_genes) >= 2) {
      train_data <- train_data[, c(target_clean, stable_genes)]
    } else {
      warning("Bootstrapped stability filter yielded fewer than 2 features; falling back to full dataset.")
    }
  }

  x_tune <- as.matrix(train_data[, setdiff(colnames(train_data), target_clean), drop = FALSE])
  y_tune <- train_data[[target_clean]]
  alpha_seq <- seq(0, 1, by = config$glm_alpha_step)
  nfolds <- max(3, config$cv_folds)
  foldid <- sample(rep(seq_len(nfolds), length.out = nrow(x_tune)))

  # For each alpha, take its own lambda.min CV error (and SE at that point).
  sweep_fits <- vector("list", length(alpha_seq))
  sweep_results <- data.frame(alpha = alpha_seq, cvm_min = NA_real_, cvsd_min = NA_real_, lambda_min = NA_real_)
  for (i in seq_along(alpha_seq)) {
    cv_fit_tune <- suppressWarnings(glmnet::cv.glmnet(
      x_tune, y_tune,
      family = fam, alpha = alpha_seq[i], standardize = TRUE, keep = TRUE, foldid = foldid
    ))
    idx_min <- which.min(cv_fit_tune$cvm)
    sweep_fits[[i]] <- cv_fit_tune
    sweep_results$cvm_min[i] <- cv_fit_tune$cvm[idx_min]
    sweep_results$cvsd_min[i] <- cv_fit_tune$cvsd[idx_min]
    sweep_results$lambda_min[i] <- cv_fit_tune$lambda[idx_min]
  }

  # 1-SE rule applied across alpha: find the global minimum, add its SE to
  # get a threshold, then keep the largest (simplest/most-regularized) alpha
  # at or below that threshold, using its own lambda.min. (The bootstrapped
  # stability filter above is unaffected by this and still uses lambda.1se
  # via config$glm_stability_rule.)
  global_idx <- which.min(sweep_results$cvm_min)
  threshold <- sweep_results$cvm_min[global_idx] + sweep_results$cvsd_min[global_idx]
  within_1se <- which(sweep_results$cvm_min <= threshold)
  best_i <- within_1se[which.max(sweep_results$alpha[within_1se])]

  best_alpha <- sweep_results$alpha[best_i]
  best_lambda <- sweep_results$lambda_min[best_i]
  best_cv_fit <- sweep_fits[[best_i]]

  lambda_idx <- which.min(abs(best_cv_fit$lambda - best_lambda))
  oof_raw <- best_cv_fit$fit.preval[, lambda_idx]
  probs_oof <- 1 / (1 + exp(-oof_raw))

  cv_preds <- data.frame(rowIndex = seq_len(nrow(x_tune)), obs = y_tune, Resample = "OOF")
  cv_preds$pred <- factor(ifelse(probs_oof > 0.5, levels(y_tune)[2], levels(y_tune)[1]), levels = levels(y_tune))
  cv_preds[[levels(y_tune)[1]]] <- 1 - probs_oof
  cv_preds[[levels(y_tune)[2]]] <- probs_oof

  cv_model <- list(pred = cv_preds, results = data.frame(alpha = best_alpha, lambda = best_lambda, ROC = 1))

  best_tune <- data.frame(alpha = best_alpha, lambda = best_lambda)
  final_trControl <- caret::trainControl(method = "none", classProbs = TRUE, summaryFunction = caret::twoClassSummary)
  model_fit <- caret::train(
    stats::as.formula(paste(target_clean, "~ .")),
    data = train_data, method = "glmnet", trControl = final_trControl,
    tuneGrid = best_tune, preProcess = config$ml_preproc, metric = "ROC"
  )

  test_predictions <- score_nerve_class(test_counts_raw, model_fit, input_type = "counts", gene_anns = gene_anns, gene_type_filter = config$gene_type_filter)
  test_predictions <- test_predictions[match(test_ids, test_predictions$`RNA-seq`), ] # align rows to test_data/test_ids order

  structure(
    list(
      model = model_fit,
      cv_model = cv_model,
      train_data = train_data,
      test_data = test_data,
      train_ids = train_ids,
      test_ids = test_ids,
      target = target_clean,
      train_cohort = train_cohort,
      test_cohort = test_cohort,
      boot_summary = boot_summary_res,
      train_cpm = train_cpm_mat[top_genes, , drop = FALSE],
      test_cpm = test_cpm_mat[top_genes, , drop = FALSE],
      pheno = rbind(train_pheno, test_pheno),
      test_predictions = test_predictions
    ),
    class = "gotnerve_fit"
  )
}
