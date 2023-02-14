#' @export

find_env.rma.uni <- function(mod) {
  if (inherits(mod$formula.yi, "formula")) {
    environment(mod$formula.yi)
  } else if (inherits(mod$formula.mods, "formula")) {
    environment(mod$formula.mods)
  } else {
    parent.frame()
  }
}

#' @export

find_env.rma.mv <- function(mod) {
  if (inherits(mod$random[[1]], "formula")) {
    environment(mod$random[[1]])
  } else if (inherits(mod$formula.yi, "formula")) {
    environment(mod$formula.yi)
  } else if (inherits(mod$formula.mods, "formula")) {
    environment(mod$formula.mods)
  } else {
    parent.frame()
  }
}

# estimate null model -----------------------------------------------------
#' @importFrom stats reformulate
#' @importFrom metafor update.rma
#' @importFrom clubSandwich findCluster.rma.mv
#' @export

estimate_null.rma <- function(full_model, C_mat) {

  # set up child environment ---------------------------------------------------
  eval_env <- find_env(full_model)
  null_env <- new.env(parent = eval_env)

  # handle formulas in yi call
  yi <- full_model$call$yi
  if (inherits(yi, "call")) full_model$call$yi <- yi[[2]]

  # Find name for null predictor matrix
  data_names <- names(eval(full_model$call$data, envir = eval_env))
  Xnull_name <- "X_null"
  while (Xnull_name %in% data_names) Xnull_name <- paste(Xnull_name, "null", sep = "_")
  mod_formula <- stats::reformulate(Xnull_name, intercept = FALSE, env = null_env)

  # compute null matrix --------------------------------------------------------

  if (is.null(full_model$subset)) {
    obs_rows <- full_model$not.na
  } else {
    obs_rows <- full_model$subset
    obs_rows[full_model$subset] <- full_model$not.na
  }

  Xnull <- matrix(NA, nrow = length(obs_rows), ncol = ncol(C_mat) - nrow(C_mat))
  Xnull[obs_rows,] <- constrain_predictors(full_model$X, C_mat)
  assign(Xnull_name, Xnull, envir = null_env)

  # estimate null model --------------------------------------------------------
  arg_list <- list(object = full_model, mods = mod_formula, evaluate = FALSE)
  null_model_call <- do.call(metafor::update.rma, args = arg_list)

  eval(null_model_call, envir = null_env)

}


# get the cluster ---------------------------------------------------------
#' @export

get_cluster.rma.mv <- function(full_model) {

  cluster <- clubSandwich::findCluster.rma.mv(full_model)

  return(cluster)

}

#' @export

get_cluster.rma.uni <- function(full_model) {

  cluster <- factor(1:full_model$k)

  return(cluster)

}



# get fitted values -------------------------------------------------------
#' @importFrom metafor fitted.rma
#' @export

get_fitted.rma <- function(model) {

  fits <- fitted.rma(model)

  return(fits)
}

# get residuals -------------------------------------------------------
#' @importFrom metafor residuals.rma
#' @export

get_res.rma <- function(model) {

  res <- residuals.rma(model)
  return(res)
}

# get indicators for complete observations----------------------------
#' @export

get_obs_rows.rma <- function(model) {

  if (is.null(model$subset)) {
    obs_rows <- model$not.na
  } else {
    obs_rows <- model$subset
    obs_rows[model$subset] <- model$not.na
  }

  return(obs_rows)
}

# get the F  --------------------------------------------------------------
#' @export

get_boot_F.rma <- function(full_model,
                              y_boot,
                              C_mat,
                              cluster,
                              type = "CR0",
                              test = "Naive-F") {

  # set up child environment ---------------------------------------------------
  eval_env <- find_env(full_model)
  boot_env <- new.env(parent = eval_env)

  # handle formulas in yi call
  yi <- full_model$call$yi
  if (inherits(yi, "call")) {
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

  if (inherits(boot_mod, "rma")) {

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

#' @export

get_boot_F_f.rma <- function(full_model,
                             C_mat,
                             cluster,
                             type = "CR0",
                             test = "Naive-F") {

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
    full_model$call$subset <- NULL
  }

  boot_env <- new.env()
  boot_env$dat <- full_model$data[obs_rows,,drop=FALSE]

  arg_list <- list(object = full_model, yi = as.symbol(y_name), data = as.symbol('dat'), evaluate = FALSE)

  if (inherits(full_model, "rma.mv")) {
    if (!is.null(full_model$call$V)) {
      boot_env$Vmat <- full_model$V
      arg_list$V <- as.symbol('Vmat')
    }
    if (!is.null(full_model$call$W)) {
      boot_env$Wmat <- full_model$W
      arg_list$W <- as.symbol('Wmat')
    }
  }

  function(y_boot, cluster = cluster) {

    boot_env$dat[[y_name]] <- y_boot

    boot_model_call <- do.call(metafor::update.rma, args = arg_list)
    boot_mod <- tryCatch(
      eval(boot_model_call, envir = boot_env),
      error = function(e) NA
    )

    if (inherits(boot_mod, "rma")) {

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

}
