#' @title Calculate bootstrap outcomes or test statistics using cluster wild
#'   bootstrapping
#'
#' @description Calculate bootstrap outcomes or test statistics using cluster
#'   wild bootstrapping for meta-analytic models fit using
#'   \code{robumeta::robu()}, \code{metafor::rma.mv()}, and \code{metafor::rma.uni()}.
#'
#' @param model Fitted \code{robumeta::robu()},
#'   \code{metafor::rma.mv()}, or \code{metafor::rma.uni()}
#'   model. For cluster wild bootstrapping, a null model is recommended, with
#'   null model indicating a model containing all variables except the ones
#'   being tested.
#' @param cluster Vector indicating which observations belong to the same
#'   cluster.
#' @param R Number of bootstrap replications.
#' @param f Optional function to be used to calculate bootstrap test statistics
#'   based on the bootstrapped outcomes. If f is \code{NULL} (the default), this
#'   function returns a list containing bootstrapped outcomes.
#' @param ... Optional arguments to be passed to the function specified in
#'   \code{f}.
#' @param auxiliary_dist Character string indicating the auxiliary distribution
#'   to be used for cluster wild bootstrapping, with available options:
#'   "Rademacher", "Mammen", "Webb six", "uniform", "standard normal". The
#'   default is set to "Rademacher." We recommend the Rademacher distribution
#'   for models that have at least 10 clusters. For models with less than 10
#'   clusters, we recommend the use of "Webb six" distribution.
#' @param adjust Character string specifying which small-sample adjustment
#'   should be used to multiply the residuals by. The available options are
#'   \code{"CRO"}, \code{"CR1"}, \code{"CR2"}, \code{"CR3"}, or \code{"CR4"},
#'   with a default of \code{"CRO"}.
#' @param simplify Logical, with \code{TRUE} indicating the bootstrapped
#'   outcomes or F statistics will be simplified to a vector or matrix and
#'   \code{FALSE} (the default) indicating the results will be returned as a
#'   list.
#' @param seed Optional seed value to ensure reproducibility.
#' @param future_args Optional list of additional arguments passed to the
#'   \code{future_*()} functions used in calculating results across bootstrap
#'   replications. Ignored if the future.apply package is not available.
#' @param future_f_args Optional list of additional arguments passed to the
#'   \code{future_*()} function used in calculating \code{f} results (but not
#'   bootstrap outcome vectors) across bootstrap replications. Ignored if the
#'   future.apply package is not available.
#'
#' @return A list or matrix containing either the bootstrapped outcomes or
#'   bootstrapped test statistics.
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
#'   model = model,
#'   cluster =  model$data.full$study,
#'   R = 12,
#'   adjust = "CR2",
#'   simplify = FALSE
#' )
#'
#' bootstraps


run_cwb <- function(model,
                    cluster,
                    R,
                    f = NULL,
                    ...,
                    auxiliary_dist = "Rademacher",
                    adjust = "CR0",
                    simplify = FALSE,
                    seed = NULL,
                    future_args = NULL,
                    future_f_args = NULL) {

  future_available <- requireNamespace("future.apply", quietly = TRUE)

  if (!is.null(seed)) set.seed(seed)

  # coerce cluster variable to factor
  if (!is.factor(cluster)) cluster <- as.factor(cluster)

  # # residuals and predicted values ------------------------------------------

  res <- get_res(model)
  pred <- get_fitted(model)

  # Adjust ------------------------------------------------------------------

  if (adjust %in% c("CR1", "CR2", "CR3", "CR4")) {
    split_res <- split(res, cluster)
    B_j <- attr(clubSandwich::vcovCR(model,
                                     cluster = cluster,
                                     type = adjust), "adjustments")
    res_list <- Map(function(x, y) as.vector(x %*% y), x = B_j, y = split_res)
    res <- unsplit(res_list, cluster)
  }

  # bootstrap ---------------------------------------------------------------
  n_clusters <- length(unique(cluster))

  replicate_args <- list(
    n = R,
    expr = quote({
      wts <- wild_wts(auxiliary_dist = auxiliary_dist, n_clusters = n_clusters)
      eta <- wts[cluster]
      y_boot <- pred + res * eta
    }),
    simplify = simplify & is.null(f)
  )

  bootstraps <- if (future_available) {
    do.call(future.apply::future_replicate,
            args = c(replicate_args, future_args))
  } else {
    do.call(replicate, args = replicate_args)
  }

  if (is.null(f)) {
    return(bootstraps)
  }

  sapply_args <- list(
    X = bootstraps,
    FUN = f,
    cluster = cluster,
    ...,
    simplify = simplify
  )

  boot_stats <- if (future_available) {
    do.call(future.apply::future_sapply,
            args = c(sapply_args, future_args, future_f_args))
  } else {
    do.call(sapply, args = sapply_args)
  }

  return(boot_stats)
}
