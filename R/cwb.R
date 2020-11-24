#' @title Calculate p-values with cluster wild bootstrapping for meta-regressions
#'
#' @description Calculates p-values for single coefficient and multiple contrast hypothesis tests using cluster wild bootstrapping.
#'
#' @param data data frame or tibble containing effect sizes, variance of effect sizes, study identification number, and moderating variables.
#' @param smd name of the column containing the standardized mean differences.
#' @param var name of the column containing the variances of the standardized mean differences.
#' @param cluster name of the column containing study identifiers.
#' @param full_model full model fit using robumeta or metafor.
#' @param null_model null model fit using robumeta or metafor. The null model includes all covariates from the full model that are not to be tested.
#' @param indices indices of the variables to be tested in single coefficient test or multiple contrast hypothesis test.
#' @param R number of bootstrap replications.
#'
#'
#'
#' @return A tibble containing the p-value.
#'
#' @export
#'
#' @examples
#'
#'
#'


cwb <- function(data,
                smd,
                var,
                cluster,
                full_model,
                null_model,
                indices,
                R = 999) {

  data$smd <- data %>%
    dplyr::pull({{smd}})

  data$var <- data %>%
    dplyr::pull({{var}})

  data$cluster <- data %>%
    dplyr::pull({{cluster}})


  # residuals and transformed residuals -------------------------------------

  data$res <- clubSandwich:::residuals_CS.robu(null_model)
  data$pred <- with(data, smd - res)
  split_res <- split(data$res, data$cluster)


  # Rademacher weights ------------------------------------------------------

  num_cluster <- unique(data$cluster)
  k_j <- as.numeric(table(data$cluster))

  system.time(

    bootstraps <- purrr::rerun(.n = R, {

      wts <- sample(c(-1, 1), size = length(num_cluster), replace = TRUE)
      data$eta <- rep(wts, k_j)
      data$new_y <- with(data, pred + res * eta)

      boot_mod <- robumeta::robu(stats::as.formula(paste("new_y ~ ",
                                        paste(full_model$reg_table$labels[!stringr::str_detect(full_model$reg_table$labels, "Intercept")],
                                              collapse = "+"))),
                       studynum = cluster,
                       var.eff.size = var,
                       small = FALSE,
                       data = data)

      cov_mat <- clubSandwich::vcovCR(boot_mod, type = "CR1")

      res <- clubSandwich::Wald_test(boot_mod,
                       constraints = clubSandwich::constrain_zero(indices),
                       vcov = cov_mat,
                       test = "Naive-F")

    }) %>%
      dplyr::bind_rows()

  )

  org_F <- clubSandwich::Wald_test(full_model,
                     constraints = clubSandwich::constrain_zero(indices),
                     vcov = clubSandwich::covCR(full_model, type = "CR1"),
                     test = "Naive-F") %>%
    dplyr::pull(Fstat)


  p_boot <- bootstraps %>%
    dplyr::summarize(p_val = mean(Fstat > org_F)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(test = "CWB") %>%
    dplyr::select(test, p_val)


  return(p_boot)

}
