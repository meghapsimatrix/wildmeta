
# estimate null model -----------------------------------------------------
#' @importFrom metafor update.rma
#' @importFrom clubSandwich findCluster.rma.mv
#' @export

estimate_null.rma.mv <- function(full_model,
                                 C_mat) {

  # set up child environment ---------------------------------------------------
  null_env <- new.env(parent = attr(full_model$random[[1]], ".Environment"))

  # handle formulas in yi call
  yi <- full_model$call$yi
  if (length(yi) > 1) full_model$call$yi <- yi[[2]]

  # compute null matrix --------------------------------------------------------

  Xnull <- constrain_predictors(full_model$X, C_mat)


  if (is.null(full_model$subset)) {
    obs_rows <- full_model$not.na
  } else {
    obs_rows <- full_model$subset
    obs_rows[full_model$subset] <- full_model$not.na
  }

  Xnull_f <- matrix(NA, nrow = length(obs_rows), ncol = ncol(Xnull))
  Xnull_f[obs_rows,] <- Xnull
  null_env$Xnull_f <- Xnull_f

  # estimate null model --------------------------------------------------------

  null_model_call <- metafor::update.rma(full_model, mods = ~ 0 + Xnull_f, evaluate = FALSE)

  eval(null_model_call, envir = null_env)

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

  # set up child environment ---------------------------------------------------
  boot_env <- new.env(parent = attr(full_model$random[[1]], ".Environment"))

  # handle formulas in yi call
  yi <- full_model$call$yi
  if (length(yi) > 1) {
    y_name <- paste(as.character(yi[[2]]), "boot", sep = "_")
    yi[[2]] <- NULL
    full_model$call$mods <- yi
  } else {
    y_name <- paste(as.character(yi), "boot", sep = "_")
  }

  if (is.null(full_model$subset)) {
    obs_rows <- full_model$not.na
  } else {
    obs_rows <- full_model$subset
    obs_rows[full_model$subset] <- full_model$not.na
  }

  y_new <- rep(NA, length = length(obs_rows))
  y_new[obs_rows] <- y_boot
  assign(y_name, y_new, envir = boot_env)

  arg_list <- list(object = full_model, yi = as.symbol(y_name), evaluate = FALSE)
  boot_model_call <- do.call(metafor::update.rma, args = arg_list)

  boot_mod <- tryCatch(eval(boot_model_call, envir = boot_env),
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
