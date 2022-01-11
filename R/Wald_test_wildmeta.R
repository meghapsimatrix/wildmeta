#' @title Calculate p-values with cluster wild bootstrapping for meta-regression
#'   models.
#'
#' @description Calculate p-values for single coefficient and multiple contrast
#'   hypothesis tests using cluster wild bootstrapping.
#'
#' @param full_model Model fit using \code{robumeta::robu()} and
#'   \code{metafor::rma.mv()} that includes the full set of moderators in the
#'   meta-regression model.
#' @param constraint_matrix A q X p constraint matrix be tested. Can be
#'   specified using \code{clubSandwich::constrain_equal()},
#'   \code{clubSandwich::constrain_zero()}, or
#'   \code{clubSandwich::constrain_pairwise()}.
#' @param R Number of bootstrap replications.
#' @param auxiliary_distribution Character string indicating the auxiliary
#'   distribution to be used for cluster wild bootstrapping, with available
#'   options: "Rademacher", "Mammen", "Webb six", "uniform", "standard normal".
#'   The default is set to "Rademacher." We recommend the Rademacher
#'   distribution for models that have at least 10 clusters. For models with
#'   less than 10 clusters, we recommend the use of "Webb six" distribution.
#' @param adjust 	Character string specifying which small-sample adjustment is
#'   used to multiply the residuals in the bootstrap process, with available
#'   options \code{"CRO"}, \code{"CR1"}, \code{"CR2"}, \code{"CR3"}, or
#'   \code{"CR4"}. The default is set to CRO, which will multiply the residuals
#'   by identity matrices and therefore, will not make any adjustments to the
#'   bootstrapping algorithm.
#' @param type Character string specifying which small-sample adjustment is used
#'   to calculate the Wald test statistic. The available options are
#'   \code{"CRO"}, \code{"CR1"}, \code{"CR2"}, \code{"CR3"}, or \code{"CR4"},
#'   with a default of \code{"CRO"}.
#' @param test Character string specifying which (if any) small-sample
#'   adjustment is used in calculating the test statistic. Default is
#'   \code{"Naive-F"}, which does not make any small-sample adjustment.
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
#' Wald_test_cwb(full_model = full_model,
#'               constraint_matrix = C_mat,
#'               R = 12)
#'
#' @importFrom clubSandwich Wald_test

Wald_test_cwb <- function(full_model,
                          constraint_matrix,
                          cluster = NULL,
                          R,
                          auxiliary_dist = "Rademacher",
                          adjust = "CR0",
                          type = "CR0",
                          test = "Naive-F") {

  # added the null model
  null_model <- estimate_null(full_model,
                              C_mat = constraint_matrix)

  # for run_cwb_new need to pull out the clusters
  # will there be an issue with missing data in clusters for rma.mv?
  if (is.null(cluster)) cluster <- get_cluster(null_model)


  boots <- run_cwb(null_model,
                   cluster = cluster,
                   R = R,
                   adjust = adjust,
                   auxiliary_dist = auxiliary_dist,
                   f = get_boot_F,  # this goes to sapply
                   full_model = full_model,
                   C_mat = constraint_matrix,
                   type = type,
                   test = test,
                   simplify = TRUE)

  full_vcov <- clubSandwich::vcovCR(full_model, type = type, cluster = cluster)
  org_F <- clubSandwich::Wald_test(full_model,
                                   constraints = C_mat,
                                   vcov = full_vcov,
                                   test = test)

  org_F <- org_F$Fstat

  p_val <- mean(boots > org_F, na.rm = TRUE)
  boot_test <- if (adjust != "CR0") "CWB Adjusted" else "CWB"

  p_boot <- data.frame(
    Test = boot_test,
    Adjustment = adjust,
    `CR type` = type,
    `Statistic` = test,
    R = R,
    p_val = p_val
  )

  class(p_boot) <- c("Wald_test_wildmeta", class(p_boot))
  attr(p_boot, "bootstraps") <- boots
  attr(p_boot, "original") <- org_F

  return(p_boot)

}
