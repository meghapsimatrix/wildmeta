run_cwb.rma.mv <- function(full_model,
                           C_mat,
                           R,
                           auxiliary_dist = "Rademacher",
                           boot_adjust = "CR0",
                           test_adjust = "CR0",
                           test_type = "chi-sq") {

  # assembling data ---------------------------------------------------------

  X_mat <- full_model$X
  effect_size <- full_model$yi
  v <- full_model$vi

  study <- clubSandwich:::findCluster.rma.mv(full_model)


  # null model --------------------------------------------------------------

  Xnull <- constrain_predictors(X_mat, C_mat)

  Xnull_f <- matrix(NA, nrow = nrow(full_model$X.f), ncol = ncol(Xnull))
  Xnull_f[full_model$not.na,] <- Xnull

  null_model <- update(full_model, yi = full_model$yi.f,  mods = ~ 0 + Xnull_f)



  # residuals and predicted values ------------------------------------------

  res <- stats::residuals(null_model)
  pred <- effect_size - res


  # Adjust ------------------------------------------------------------------

  if (adjust == TRUE) {
    split_res <- split(res, study)
    # JEP: Is the as.matrix() coercion necessary?
    e_tilde_j <- purrr::map(split_res, as.matrix)
    B_j <- attr(clubSandwich::vcovCR(null_model,
                                     cluster = study,
                                     type = "CR2",
                                     inverse_var = TRUE), "adjustments")

    # JEP: What if the cluster ID is not in sort order?
    # I think might need to use unsplit() here.
    res <- unlist(purrr::pmap(list(B_j, e_tilde_j), function(x, y) as.vector(x %*% y)))

  }

  # bootstrap ---------------------------------------------------------------
  num_cluster <- unique(study)
  k_j <- as.numeric(table(study))

  bootstraps <- replicate(n = R, {

    wts <- return_wts(auxiliary_dist = auxiliary_dist, cluster_var = num_cluster)
    # JEP: What if the clusters are not in sort order?
    eta <- rep(wts, k_j)
    y_boot <- pred + res * eta

    # JAMES for missing data
    y_new <- rep(NA, length = nrow(full_model$X.f))
    y_new[full_model$not.na] <- y_boot

    boot_mod <- update(full_model, yi = y_new)

    cov_mat <- clubSandwich::vcovCR(boot_mod, type = "CR1")

    res <- clubSandwich::Wald_test(boot_mod,
                                   constraints = C_mat,
                                   vcov = cov_mat,
                                   test = "Naive-F") # test-type?

  }, simplify = FALSE) %>%
    dplyr::bind_rows()

  booties <- bootstraps$Fstat


  class(booties) <- "wildmeta"

  return(booties)

}
