#' @title Calculate bootstrap outcomes or test statistics using cluster wild bootstrapping
#'
#' @description Calculate bootstrap outcomes or test statistics using cluster wild bootstrapping for meta-analytic models fit using robu() and rma.mv()
#'
#' @param full_model Model fit using `robu()` or `rma.mv()`. For cluster wild bootstrapping, a null model is recommended, with null model indicating a model containing all variables except the ones being tested.
#' @param cluster Vector indicating the clustering variable.
#' @param f Optional name (unquoted) of a function to be used to calculate bootstrap test statistics based on the bootstrapped outcomes. Default value is NULL. If f is NULL, this function returns a list with length of the number of bootstraps containing bootstrapped outcomes.
#' @param ... Optional arguments to be passed to the function specified in the f argument above.
#' @param auxiliary_distribution Character string indicating the auxiliary distribution to be used for cluster wild bootstrapping, with available options: "Rademacher", "Mammen", "Webb six", "uniform", "standard normal". The default is set to "Rademacher." We recommend the Rademacher distribution for models that have at least 10 clusters. For models with less than 10 clusters, we recommend the use of "Webb six" distribution.
#' @param adjust 	Character string specifying which small-sample adjustment should be used to multiply the residuals by, with available options "CR0", "CR1", "CR2", "CR3", or "CR4". The default is set to CRO, which will multiply the residuals by identity matrices and therefore, wil not add any adjustments to the boostrapping algorithm.
#' @param simplify Logical, with TRUE indicating the bootstrapped outcomes or F statistics will be simplified to a vector or matrix and FALSE indicating the results will be returned as a list.
#'
#'
#' @return A list or matrix containing either the bootstrapped outcomes or bootstrapped test statistics.
#'
#' @export
#'
#' @examples
#' library(clubSandwich)
#' library(robumeta)
#'
#' model <- robu(d ~ 0 + study_type + hrs + test,
#'              studynum = study,
#'               var.eff.size = V,
#'               small = FALSE,
#'               data = SATcoaching)
#'
#'
#' bootstraps <- run_cwb(
#'   model = full_model,
#'   cluster =  full_model$data.full$study,
#'   R = 12,
#'   adjust = "CR2",
#'   simplify = FALSE
#' )
#'
#'
#'


run_cwb <- function(model,
                    cluster,
                    f = NULL,
                    ...,
                    R,
                    auxiliary_dist = "Rademacher",
                    adjust = "CR0",
                    simplify = FALSE) {


  # coerce cluster variable to factor
  if (!is.factor(cluster)) cluster <- as.factor(cluster)

  # # residuals and predicted values ------------------------------------------

  # need to figure out the robu situation
  # model$fitted.values <- fitted.robu(model)
  # model$residuals <- residuals.robu(model)

  res <- get_res(model)
  pred <- get_fitted(model)

  # Adjust ------------------------------------------------------------------

  if (adjust %in% c("CR1", "CR2", "CR3", "CR4")) {
    split_res <- split(res, cluster)
    B_j <- attr(clubSandwich::vcovCR(model,
                                     cluster = cluster,
                                     type = adjust), "adjustments")
    res_list <- purrr::map2(B_j, split_res, ~ as.vector(.x %*% .y))
    res <- unsplit(res_list, cluster)
  }

  # bootstrap ---------------------------------------------------------------
  num_cluster <- unique(cluster)

  bootstraps <- replicate(n = R, {

    wts <- return_wts(auxiliary_dist = auxiliary_dist, cluster_var = num_cluster)
    eta <- wts[cluster]
    y_boot <- pred + res * eta

    #y_new <- rep(NA, length = nrow(model$X.f))  # how is this working with robu?
    #y_new[model$not.na] <- y_boot

  }, simplify = simplify & is.null(f))

  if (is.null(f)) {
    return(bootstraps)
  }

  boot_stats <- sapply(bootstraps, f, ..., simplify = simplify) # lapply doesn;t have simplify as argument

  return(boot_stats)
}