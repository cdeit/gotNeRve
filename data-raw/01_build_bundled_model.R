# Trains the model from raw counts and ships it as internal (non-exported)
# package data in R/sysdata.rda, used as gotNeRve()'s default `model` so
# gotNeRve(new_data) works with zero setup.
#
# Run this whenever train_model()'s algorithm changes and the bundled model
# needs to be refreshed.

devtools::load_all(".", quiet = TRUE)

e1 <- new.env()
load("data-raw/train_test_counts_and_clin.RData", envir = e1)

fit <- train_model(e1$pheno_data, e1$a_counts, e1$b_counts, e1$gene_anns, gotnerve_config())

message(sprintf(
  "Trained: alpha = %s, lambda = %.6g, %d signature genes.",
  fit$model$bestTune$alpha, fit$model$bestTune$lambda, length(fit$model$coefnames)
))

gotnerve_model <- fit$model
usethis::use_data(gotnerve_model, internal = TRUE, overwrite = TRUE)
