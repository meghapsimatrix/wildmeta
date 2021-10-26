run_CWB.robu <- function(full_model,
                         C_mat,
                         R,
                         auxiliary_dist = "Rademacher",
                         adjust = FALSE) {


  # info about model --------------------------------------------------------

  dep <- full_model$modelweights
  intercept <- sum(str_detect(full_model$reg_table[, 1], "X.Intercept."))

  # assembling data ---------------------------------------------------------

  full_dat <- full_model$data.full %>%
    dplyr::mutate(id = rownames(.)) %>%
    dplyr::rename(effect_size = 1,
                  v = 2)

  x_dat <- full_model$X.full %>%
    dplyr::mutate(id = rownames(.)) %>%
    dplyr::select(-1)

  dat <- full_dat %>%
    dplyr::left_join(x_dat, by = "id")


  # full formula ------------------------------------------------------------

  full_formula <- paste(full_model$reg_table[, 1], collapse = " + ")
  full_formula <- stringr::str_replace(full_formula, "X.Intercept.", "1")


  # null_model --------------------------------------------------------------

  X_mat <- full_model$X.full %>%
    select(-1) %>%
    as.matrix()

  Xnull <- constrain_predictors(Xmat = X_mat, Cmat = C_mat)

  dat <- bind_cols(dat, as.data.frame(Xnull))

  null_formula <- paste("effect_size ~ 0 + ", paste(colnames(as.data.frame(Xnull)), collapse = " + "))

  null_model <- robumeta::robu(stats::as.formula(null_formula),
                               studynum = study,
                               var.eff.size = v,
                               small = FALSE,
                               modelweights = dep,
                               data = dat)

  # residuals and predicted values ------------------------------------------

  dat$res <- clubSandwich:::residuals_CS.robu(null_model)
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

    # JAMES intercept sitation here :)
    if(intercept == 0){

      boot_formula <- stats::as.formula(paste("new_y ~ 0 + ", full_formula))

    } else {

      boot_formula <- stats::as.formula(paste("new_y ~", full_formula))

    }

    boot_mod <- robumeta::robu(boot_formula,
                               studynum = study,
                               var.eff.size = v,
                               small = FALSE,
                               modelweights = dep,
                               data = dat)

    cov_mat <- clubSandwich::vcovCR(boot_mod, type = "CR1")

    res <- clubSandwich::Wald_test(boot_mod,
                                   constraints = C_mat,
                                   vcov = cov_mat,
                                   test = "Naive-F")

  }, simplify = FALSE) %>%
    dplyr::bind_rows()

  booties <- bootstraps$Fstat


  class(booties) <- "wildmeta"

  return(booties)

}
