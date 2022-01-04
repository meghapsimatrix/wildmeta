#' @title Calculate p-values with cluster wild bootstrapping for meta-regression models.
#'
#' @description Calculates p-values for single coefficient and multiple contrast hypothesis tests using cluster wild bootstrapping.
#'
#' @param full_model Model fit using `robu()` or `rma.mv()` that includes all moderators of interest.
#' @param constraint_matrix A q X p constraint matrix be tested. Can be specified using constrain_equal, constrain_zero, or constrain_pairwise from the clubSandwich package.
#' @param R number of bootstrap replications.
#' @param auxiliary_distribution Character string indicating the auxiliary distribution to be used for cluster wild bootstrapping, with available options: "Rademacher", "Mammen", "Webb six", "uniform", "standard normal". The default is set to "Rademacher." We recommend the Rademacher distribution for models that have at least 10 clusters. For models with less than 10 clusters, we recommend the use of "Webb six" distribution.
#' @param adjust 	Character string specifying which small-sample adjustment should be used to multiply the residuals by, with available options "CR0", "CR1", "CR2", "CR3", or "CR4". The default is set to CRO, which will multiply the residuals by identity matrices and therefore, wil not add any adjustments to the boostrapping algorithm.
#'
#'
#' @return A tibble containing the name of the test and the p-value.
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

Wald_test_cwb <- function(full_model,
                          constraint_matrix,
                          R,
                          auxiliary_dist = "Rademacher",
                          adjust = "CR0"){

  # added the null model
  null_model <- estimate_null(full_model,
                              C_mat = constraint_matrix)

  # for run_cwb_new need to pull out the clusters
  # will there be an issue with missing data in clusters for rma.mv?
  cluster <- get_cluster(null_model)


  boots <- run_cwb(null_model,
                   cluster = cluster,
                   R = R,
                   adjust = adjust,
                   auxiliary_dist = auxiliary_dist,
                   f = get_boot_F,  # this goes to sapply
                   full_model = full_model,
                   C_mat = constraint_matrix, # this is additional argument for sapply
                   simplify = TRUE)

  org_F <- clubSandwich::Wald_test(full_model,
                                   constraints = C_mat,
                                   vcov = clubSandwich::vcovCR(full_model, type = "CR1"),
                                   test = "Naive-F")

  org_F <- org_F$Fstat

  p_val <- mean(boots > org_F, na.rm = TRUE)
  test <- "CWB"


  if (adjust != "CR0") {
    test <- "CWB Adjusted"
  }

  p_boot <- data.frame(test = test, p_val = p_val)

  class(p_boot) <- c("Wald_test_wildmeta", class(p_boot))
  attr(p_boot, "bootstraps") <- boots
  attr(p_boot, "original") <- org_F



  return(p_boot)

}
