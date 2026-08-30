#' Configuration for the gotNeRve default pipeline
#'
#' Returns a named list of settings consumed by [train_model()]. Defaults
#' match the pipeline as originally run, with one intentional change:
#' `min_count = 8` (originally 10).
#'
#' @param target_var Phenotype column in `pheno_data` to predict.
#' @param quantile_threshold Extreme-quantile cutoff used to re-label
#'   `nerve_*` targets into High/Low, dropping the middle tier. Set to `NULL`
#'   to skip this step for non-`nerve_*` targets.
#' @param train_cohort Which cohort ("A" or "B") trains the model; the other
#'   cohort becomes the held-out test set.
#' @param gene_type_filter `GENE_TYPE` values (exact match against
#'   `gene_anns$GENE_TYPE`) to keep before normalization.
#' @param min_count,min_samples Passed to `edgeR::filterByExpr()`.
#' @param top_n_var Number of highest-variance training-cohort genes kept.
#' @param ml_preproc `caret::train()` `preProcess` methods for the final fit.
#' @param cv_folds Number of folds used by `glmnet::cv.glmnet()`.
#' @param glm_alpha_step Increment for the elastic-net alpha sweep, from 0 to
#'   1 inclusive. The alpha/lambda combination is chosen by a 1-SE rule
#'   applied across alpha: at each alpha, take its own `lambda.min` CV error;
#'   find the global minimum of those; add its standard error to get a
#'   threshold; among alphas at or below that threshold, keep the largest
#'   (most-regularized/simplest) one, using its own `lambda.min`. (The
#'   bootstrapped stability filter is separate and always uses
#'   `lambda.1se`, via `glm_stability_rule`.)
#' @param glm_stability Whether to apply the bootstrapped feature-stability
#'   filter before the alpha/lambda tuning sweep.
#' @param glm_bootstraps Number of stability-filter bootstrap resamples.
#' @param glm_stability_threshold Minimum bootstrap inclusion frequency for a
#'   gene to survive the stability filter.
#' @param glm_stability_sign_threshold Minimum coefficient-sign consistency
#'   (across bootstraps in which a gene was non-zero) required to survive.
#' @param glm_stability_alpha Fixed alpha used only during the stability
#'   bootstrap (independent of `glm_alpha_range`).
#' @param glm_stability_rule `"oneSE"` or `"best"`, used to pick the fixed
#'   lambda for the stability bootstrap.
#'
#' @return A named list.
gotnerve_config <- function(target_var = "nerve_group",
                            quantile_threshold = 0.33,
                            train_cohort = "A",
                            gene_type_filter = c("protein_coding", "miRNA", "lncRNA"),
                            min_count = 8,
                            min_samples = 5,
                            top_n_var = 1000,
                            ml_preproc = c("center", "scale"),
                            cv_folds = 10,
                            glm_alpha_step = 0.1,
                            glm_stability = TRUE,
                            glm_bootstraps = 1000,
                            glm_stability_threshold = 0.25,
                            glm_stability_sign_threshold = 0.50,
                            glm_stability_alpha = 0.5,
                            glm_stability_rule = "oneSE") {
  list(
    target_var = target_var,
    quantile_threshold = quantile_threshold,
    train_cohort = train_cohort,
    gene_type_filter = gene_type_filter,
    min_count = min_count,
    min_samples = min_samples,
    top_n_var = top_n_var,
    ml_preproc = ml_preproc,
    cv_folds = cv_folds,
    glm_alpha_step = glm_alpha_step,
    glm_stability = glm_stability,
    glm_bootstraps = glm_bootstraps,
    glm_stability_threshold = glm_stability_threshold,
    glm_stability_sign_threshold = glm_stability_sign_threshold,
    glm_stability_alpha = glm_stability_alpha,
    glm_stability_rule = glm_stability_rule
  )
}
