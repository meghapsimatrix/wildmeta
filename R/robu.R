
# estimate null model -----------------------------------------------------
#' @importFrom robumeta robu
#' @export

estimate_null.robu <- function(full_model,
                               C_mat) {

  ord <- order(order(full_model$study_orig_id))

  dep <- full_model$modelweights

  # assembling data ---------------------------------------------------------

  es_dat <- full_model$data.full[ord, c("effect.size", "var.eff.size", "study")]

  # null_model --------------------------------------------------------------

  X_mat <- as.matrix(full_model$X.full[ord, -1])
  Xnull <- as.data.frame(constrain_predictors(Xmat = X_mat, Cmat = C_mat))

  null_dat <- cbind(es_dat, Xnull)

  null_formula <- paste("effect.size ~ 0 + ", paste(colnames(Xnull), collapse = " + "))

  null_model <- robumeta::robu(stats::as.formula(null_formula),
                               studynum = study,
                               var.eff.size = var.eff.size,
                               small = FALSE,
                               modelweights = dep,
                               data = null_dat)

  return(null_model)


}

# get the cluster ---------------------------------------------------------
#' @export

get_cluster.robu <- function(full_model) {

  ord <- order(order(full_model$study_orig_id))
  cluster <- full_model$data.full$study[ord]

  return(cluster)
}

# get the F  --------------------------------------------------------------
#' @export

get_boot_F.robu <- function(full_model,
                            y_boot,
                            C_mat,
                            cluster,
                            type = "CR0",
                            test = "Naive-F") {

  # use update robu to fit bootstrapped model

  boot_mod <- update_robu(full_model,
                          y = y_boot)

  cov_mat <- clubSandwich::vcovCR(boot_mod, cluster = cluster, type = type)

  res <- clubSandwich::Wald_test(boot_mod,
                                 constraints = C_mat,
                                 vcov = cov_mat,
                                 test = test)

  res <- res$Fstat

  return(res)

}

#' @export

get_boot_F_f.robu <- function(full_model,
                              C_mat,
                              cluster,
                              type = "CR0",
                              test = "Naive-F") {

  function(y_boot) {

    # use update robu to fit bootstrapped model

    boot_mod <- update_robu(full_model,
                            y = y_boot)

    cov_mat <- clubSandwich::vcovCR(boot_mod, cluster = cluster, type = type)

    res <- clubSandwich::Wald_test(boot_mod,
                                   constraints = C_mat,
                                   vcov = cov_mat,
                                   test = test)

    res <- res$Fstat

    return(res)

  }

}

# get fitted values -------------------------------------------------------
#' @export

get_fitted.robu <- function(model) {

  ord <- order(order(model$study_orig_id))
  fits <- as.numeric(model$data.full$pred[ord])

  return(fits)
}

# get residuals -------------------------------------------------------
#' @export

get_res.robu <- function(model) {

  ord <- order(order(model$study_orig_id))
  res <- model$data.full$e.r[ord]

  return(res)
}

# get model coefficients ---------------------------------------------
#' @export

coef.robu <- function(object, ...) {
  cf <- object$reg_table$b.r
  names(cf) <- object$reg_table$labels
  cf[!is.na(cf)]
}
