run_CWB.rma.mv <- function(full_model,
                           C_mat,
                           R,
                           auxiliary_dist = "Rademacher",
                           adjust = FALSE) {

  # JAMES
  # need something here to make the constrain matrix if users add indices?
  # like 1:3 in C_mat instead of constrain_equal(...)?

  X_mat <- full_model$X
  effect_size <- full_model$yi
  v <- full_model$vi

  study <- clubSandwich:::findCluster.rma.mv(full_model)

  intercept <- sum(str_detect(rownames(full_model$beta), "intrcpt"))

  # Null model --------------------------------------------------------------

  Xnull <- constrain_predictors(X_mat, C_mat)

  Xnull_f <- matrix(NA, nrow = nrow(full_model$X.f), ncol = ncol(Xnull))
  Xnull_f[full_model$not.na,] <- Xnull

  # JAMES - does it matter if we do 0 here for all kinds of models or no?
  # likes if original model has an intercept
  # I think it doesn't matter bc we are not doing constraint test but just wanted to check
  null_model <- update(full_model, yi = full_model$yi.f,  mods = ~ 0 + Xnull_f)


  res <- stats::residuals(null_model)


  # residuals and transformed residuals -------------------------------------
  pred <- effect_size - res
  split_res <- split(res, study)


  # Adjust ------------------------------------------------------------------

  if(adjust == TRUE){
    e_tilde_j <- purrr::map(split_res, as.matrix)
    B_j <- attr(clubSandwich::vcovCR(null_model,
                                     cluster = study,
                                     type = "CR2",
                                     inverse_var = TRUE), "adjustments")
    res <- unlist(purrr::pmap(list(B_j, e_tilde_j), function(x, y) as.vector(x %*% y)))

  }

  # bootstrap ---------------------------------------------------------------
  num_cluster <- unique(study)
  k_j <- as.numeric(table(study))

  bootstraps <- replicate(n = R, {

    wts <- return_wts(auxiliary_dist = auxiliary_dist, cluster_var = num_cluster)
    eta <- rep(wts, k_j)
    y_boot <- pred + res * eta

    # JAMES for missing data
    y_new <- rep(NA, length = nrow(full_model$X.f))
    y_new[full_model$not.na] <- y_boot

    # JAMES CHECK THIS:
    # some way to do this so it's not dropping intercept
    # so it can use the original constrain matrix for F tests below?

    if(intercept == 0){

    boot_mod <- update(full_model, yi = y_new, mods = ~ 0 + full_model$X.f)

    } else{

    boot_mod <- update(full_model, yi = y_new, mods = ~ full_model$X.f)

    }


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
