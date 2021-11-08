
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

  dep <- full_model$modelweights

  # assembling data ---------------------------------------------------------

  es_dat <- full_model$data.full[, c("effect.size", "var.eff.size", "study")]

  # null_model --------------------------------------------------------------

  X_mat <- as.matrix(full_model$X.full[, -1])
  Xnull <- constrain_predictors(Xmat = X_mat, Cmat = C_mat)

  null_dat <- cbind(es_dat, as.data.frame(Xnull))

  null_formula <- paste("effect.size ~ 0 + ", paste(colnames(as.data.frame(Xnull)), collapse = " + "))

  null_model <- robumeta::robu(stats::as.formula(null_formula),
                               studynum = study,
                               var.eff.size = var.eff.size,
                               small = FALSE,
                               modelweights = dep,
                               data = null_dat)

  return(null_model)


}


# get the cluster ---------------------------------------------------------

get_cluster.robu <- function(full_model){

  cluster <- full_model$data.full$study

  return(cluster)
}


get_cluster.rma.mv <- function(full_model){

  cluster <- clubSandwich:::findCluster.rma.mv(full_model)

  return(cluster)

}


# get the F  --------------------------------------------------------------

get_boot_F.robu <- function(y_boot,
                            full_model,
                            C_mat){

  # info about model --------------------------------------------------------

  dep <- full_model$modelweights
  intercept <- sum(stringr::str_detect(full_model$reg_table[, 1], "X.Intercept."))

  dat <- full_model$data.full
  dat$new_y <- y_boot

  # full formula ------------------------------------------------------------

  full_formula <- paste(full_model$reg_table[, 1], collapse = " + ")
  full_formula <- stringr::str_replace(full_formula, "X.Intercept.", "1")



  # restimate the full model ------------------------------------------------

  if(intercept == 0){

    boot_formula <- stats::as.formula(paste("new_y ~ 0 + ", full_formula))

  } else {

    boot_formula <- stats::as.formula(paste("new_y ~", full_formula))

  }

  boot_mod <- robumeta::robu(boot_formula,
                             studynum = study,
                             var.eff.size = var.eff.size,
                             small = FALSE,
                             modelweights = dep,
                             data = dat)

  cov_mat <- clubSandwich::vcovCR(boot_mod, type = "CR1")

  res <- clubSandwich::Wald_test(boot_mod,
                                 constraints = C_mat,
                                 vcov = cov_mat,
                                 test = "Naive-F")

  res <- res$Fstat

  return(res)



}


get_boot_F.rma.mv <- function(y_boot,
                              full_model,
                              C_mat){


  y_new <- rep(NA, length = nrow(full_model$X.f))
  y_new[full_model$not.na] <- y_boot


  boot_mod <- tryCatch(update(full_model, formula = y_new ~ .),
                       error = function(e) NA)

  if(!is.na(boot_mod)){

  cov_mat <- clubSandwich::vcovCR(boot_mod, type = "CR1")

  res <- clubSandwich::Wald_test(boot_mod,
                                 constraints = C_mat,
                                 vcov = cov_mat,
                                 test = "Naive-F")

  res <- res$Fstat


  } else{

    res <- NA
  }

  return(res)



}



# get fitted values -------------------------------------------------------


get_fitted.robu <- function(model){

  fits <- fitted.robu(null_model)  # where is this from?

  return(fits)
}

get_fitted.rma.mv <- function(model){

  fits <- stats::fitted.values(model)

  return(fits)
}


get_res.robu <- function(model){

  res <- residuals.robu(full_model)
  return(res)
}

get_res.rma.mv <- function(model){


  res <- stats::residuals(model)
  return(res)
}
