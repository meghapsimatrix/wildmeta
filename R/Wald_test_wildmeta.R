#' @title Calculate p-values with cluster wild bootstrapping for meta-regression
#'   models.
#'
#' @description Calculate p-values for single coefficient and multiple contrast
#'   hypothesis tests using cluster wild bootstrapping.
#'
#' @param full_model Model fit using \code{robumeta::robu()},
#'   \code{metafor::rma.mv()}, or \code{metafor::rma.uni()} that includes the full set of moderators in the
#'   meta-regression model.
#' @param constraints A q X p constraint matrix be tested. Alternately, a
#'   function to create such a matrix, specified using
#'   \code{clubSandwich::constrain_equal()} or
#'   \code{clubSandwich::constrain_zero()}.
#' @param cluster Vector of identifiers indicating which observations
#'   belong to the same cluster. If \code{NULL} (the default), then the
#'   clustering variable will be inferred based on the structure of
#'   \code{full_mod}.
#' @param type Character string specifying which small-sample adjustment is used
#'   to calculate the Wald test statistic. The available options are
#'   \code{"CRO"}, \code{"CR1"}, \code{"CR2"}, \code{"CR3"}, or \code{"CR4"},
#'   with a default of \code{"CRO"}.
#' @param test Character string specifying which (if any) small-sample
#'   adjustment is used in calculating the test statistic. Default is
#'   \code{"Naive-F"}, which does not make any small-sample adjustment.
#' @inheritParams run_cwb
#'
#' @return A \code{data.frame} containing the name of the test, the adjustment
#'   used for the bootstrap process, the type of variance-covariance matrix
#'   used, the type of test statistic, the number of bootstrap replicates, and
#'   the bootstrapped p-value.
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
#' C_mat <- constrain_equal(1:3, coefs = coef(model))
#'
#' Wald_test_cwb(full_model = model,
#'               constraints = C_mat,
#'               R = 12)
#'
#' # Equivalent, using constrain_equal()
#' Wald_test_cwb(full_model = model,
#'               constraints = constrain_equal(1:3),
#'               R = 12)
#'
#' @importFrom clubSandwich Wald_test
#' @importFrom stats coef

Wald_test_cwb <- function(full_model,
                          constraints,
                          R,
                          cluster = NULL,
                          auxiliary_dist = "Rademacher",
                          adjust = "CR0",
                          type = "CR0",
                          test = "Naive-F",
                          seed = NULL,
                          future_args = NULL) {

  if (inherits(constraints, "function")) {
    constraints <- constraints(stats::coef(full_model))
  }

  # compute the null model
  null_model <- estimate_null(full_model,
                              C_mat = constraints)

  # detect clusters if not specified
  if (is.null(cluster)) cluster <- get_cluster(null_model)

  # evaluate f on each bootstrap
  future_f_args <- if (inherits(full_model,"rma")) {
    list(
      future.packages = c("clubSandwich","metafor"),
      future.envir = find_env(full_model)
    )
  } else {
    NULL
  }

  boots <- run_cwb(null_model,
                   cluster = cluster,
                   R = R,
                   f = get_boot_F,
                   full_model = full_model,
                   C_mat = constraints, type = type, test = test,
                   auxiliary_dist = auxiliary_dist,
                   adjust = adjust,
                   simplify = TRUE,
                   seed = seed,
                   future_args = future_args,
                   future_f_args = future_f_args)

  full_vcov <- clubSandwich::vcovCR(full_model, type = type, cluster = cluster)
  org_F <- clubSandwich::Wald_test(full_model,
                                   constraints = constraints,
                                   vcov = full_vcov,
                                   test = test)

  org_F <- org_F$Fstat

  p_val <- mean(boots > org_F, na.rm = TRUE)
  boot_test <- if (adjust != "CR0") "CWB Adjusted" else "CWB"

  p_boot <- data.frame(
    Test = boot_test,
    Adjustment = adjust,
    CR_type = type,
    Statistic = test,
    R = R,
    p_val = p_val
  )

  class(p_boot) <- c("Wald_test_wildmeta", class(p_boot))
  attr(p_boot, "bootstraps") <- boots
  attr(p_boot, "original") <- org_F

  return(p_boot)

}
