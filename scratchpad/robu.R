run_CWB.robu <- function(full_model,
                         C_mat,
                         R,
                         auxiliary_dist = "Rademacher",
                         adjust = FALSE) {

  dep <- full_model$modelweights

  full_dat <- full_model$data.full %>%
    dplyr::mutate(id = rownames(.)) %>%
    dplyr::rename(effect_size = 1,
                  v = 2)

  X_mat <- full_model$X.full %>%
    select(-1) %>%
    as.matrix()

  full_formula <- paste(full_model$reg_table[, 1], collapse = " + ")
  full_formula <- stringr::str_replace(full_formula, "X.Intercept.", "1")

  Xnull <- constrain_predictors(Xmat = X_mat, Cmat = C_mat)


  null_model <- robumeta::robu(stats::as.formula(paste("effect_size ~ ", null_formula)),
                               studynum = study,
                               var.eff.size = v,
                               small = FALSE,
                               modelweights = dep,
                               data = dat)

  dat$res <- clubSandwich:::residuals_CS.robu(null_model)

  # residuals and transformed residuals -------------------------------------
  dat$pred <- with(dat, effect_size - res)
  split_res <- split(dat$res, dat$study)

  # Adjust ------------------------------------------------------------------

  if(adjust == TRUE){
    e_tilde_j <- purrr::map(split_res, as.matrix)
    B_j <- attr(clubSandwich::vcovCR(null_model,
                                     cluster = dat$study,
                                     type = "CR2",
                                     inverse_var = TRUE), "adjustments")
    dat$res <- unlist(purrr::pmap(list(B_j, e_tilde_j), function(x, y) as.vector(x %*% y)))

  }

  # bootstrap ---------------------------------------------------------------
  num_cluster <- unique(dat$study)
  k_j <- as.numeric(table(dat$study))

  bootstraps <- replicate(n = R, {

    wts <- return_wts(auxiliary_dist = auxiliary_dist, cluster_var = num_cluster)
    dat$eta <- rep(wts, k_j)
    dat$new_y <- with(dat, pred + res * eta)

    boot_mod <- robumeta::robu(stats::as.formula(paste("new_y ~ ", full_formula)),
                               studynum = study,
                               var.eff.size = v,
                               small = FALSE,
                               modelweights = dep,
                               data = dat)

    cov_mat <- clubSandwich::vcovCR(boot_mod, type = "CR1")

    res <- clubSandwich::Wald_test(boot_mod,
                                   constraints = clubSandwich::constrain_zero(indices),
                                   vcov = cov_mat,
                                   test = "Naive-F")

  }, simplify = FALSE) %>%
    dplyr::bind_rows()

  class(bootstraps) <- "wildmeta"

  return(bootstraps)

}
