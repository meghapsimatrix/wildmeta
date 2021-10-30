
# Constrain Predictors ----------------------------------------------------

constrain_predictors <- function(Xmat, Cmat) {

  q <- nrow(Cmat)
  p <- ncol(Cmat)
  if (ncol(Xmat) != ncol(Cmat)) stop("Constraint matrix must have same number of columns as predictor matrix.")

  XtX_inv <- chol2inv(chol(crossprod(Xmat)))
  Cnull <- diag(nrow = p) - XtX_inv %*% t(Cmat) %*% chol2inv(chol(Cmat %*% XtX_inv %*% t(Cmat))) %*% Cmat
  Xnull <- qr.X(qr(Xmat %*% Cnull), ncol = p - q)

  return(Xnull)

}



# Auxiliary distribution --------------------------------------------------

return_wts <- function(auxiliary_dist, cluster_var) {

  if(auxiliary_dist == "Rademacher"){

    wts <- sample(c(-1, 1),
                  size = length(cluster_var),
                  replace = TRUE,
                  prob = rep(1/2, 2))

  } else if(auxiliary_dist == "Mammen"){

    wts <- sample(c(- (sqrt(5) - 1)/2,
                    (sqrt(5) + 1)/2),
                  size = length(cluster_var),
                  replace = TRUE,
                  prob = c((sqrt(5) + 1) /(2 * sqrt(5)), (sqrt(5) - 1) /(2 * sqrt(5)) ))

  } else if(auxiliary_dist == "Webb six"){

    wts <- sample(c(-sqrt(3/2),
                    -sqrt(2/2),
                    -sqrt(1/2),
                    sqrt(1/2),
                    sqrt(2/2),
                    sqrt(3/2)),
                  size = length(cluster_var),
                  replace = TRUE,
                  prob = rep(1/6, 6))

  } else if(auxiliary_dist == "uniform"){

    wts <- runif(n = length(cluster_var),
                 min = -sqrt(3),
                 max = sqrt(3))

  } else if(auxiliary_dist == "standard normal"){

    wts <- rnorm(n = length(cluster_var))

  }

  return(wts)

}




# estimate null model -----------------------------------------------------

estimate_null.rma.mv <- function(full_model,
                                 C_mat,
                                 R) {

  # info -------------------------------------------------------------------

  X_mat <- full_model$X
  cluster <- clubSandwich:::findCluster.rma.mv(full_model)


  # null model --------------------------------------------------------------

  Xnull <- constrain_predictors(X_mat, C_mat)

  Xnull_f <- matrix(NA, nrow = nrow(full_model$X.f), ncol = ncol(Xnull))
  Xnull_f[full_model$not.na,] <- Xnull

  null_model <- update(full_model, yi = full_model$yi.f,  mods = ~ 0 + Xnull_f)

  return(null_model)

}


estimate_null.robu <- function(full_model,
                               C_mat,
                               R) {


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

  cluster <- full_model$data.full$study


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

  return(null_model)


}
