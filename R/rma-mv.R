
# estimate null model -----------------------------------------------------
#' @importFrom metafor update.rma
#' @importFrom clubSandwich findCluster.rma.mv
#' @export

estimate_null.rma.mv <- function(full_model,
                                 C_mat) {

  # info -------------------------------------------------------------------

  X_mat <- full_model$X
  cluster <- clubSandwich::findCluster.rma.mv(full_model)


  # null model --------------------------------------------------------------

  Xnull <- constrain_predictors(X_mat, C_mat)

  Xnull_f <- matrix(NA, nrow = nrow(full_model$X.f), ncol = ncol(Xnull))
  Xnull_f[full_model$not.na,] <- Xnull

  null_model <- metafor::update.rma(full_model, yi = full_model$yi.f,  mods = ~ 0 + Xnull_f)

  return(null_model)

}


# get the cluster ---------------------------------------------------------
#' @export

get_cluster.rma.mv <- function(full_model) {

  cluster <- clubSandwich::findCluster.rma.mv(full_model)

  return(cluster)

}


# get the F  --------------------------------------------------------------
#' @export

get_boot_F.rma.mv <- function(full_model,
                              y_boot,
                              C_mat,
                              cluster,
                              type = "CR0",
                              test = "Naive-F") {

  new_dat <- eval(full_model$call$data)

  if (is.null(full_model$formula.yi)) {
    yi_name <- as.character(full_model$call$yi)
  } else {
    yi_name <- as.character(full_model$formula.yi[[2]])
  }

  y_new <- rep(NA, length = nrow(full_model$X.f))
  y_new[full_model$not.na] <- y_boot

  new_dat[[yi_name]] <- y_new


  boot_mod <- tryCatch(metafor::update.rma(full_model, data = new_dat),
                       error = function(e) NA)

  if (inherits(boot_mod, "rma.mv")) {

    cov_mat <- clubSandwich::vcovCR(boot_mod, cluster = cluster, type = type)

    res <- clubSandwich::Wald_test(boot_mod,
                                   constraints = C_mat,
                                   vcov = cov_mat,
                                   test = test)

    res <- res$Fstat


  } else {

    res <- NA
  }

  return(res)

}

# get fitted values -------------------------------------------------------
#' @importFrom metafor fitted.rma
#' @export

get_fitted.rma.mv <- function(model){

  fits <- metafor::fitted.rma(model)

  return(fits)
}

# get residuals -------------------------------------------------------
#' @importFrom metafor residuals.rma
#' @export

get_res.rma.mv <- function(model) {

  res <- metafor::residuals.rma(model)
  return(res)
}
