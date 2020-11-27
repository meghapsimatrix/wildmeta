#' @title Calculate p-values with cluster wild bootstrapping for meta-regressions
#'
#' @description Calculates p-values for single coefficient and multiple contrast hypothesis tests using cluster wild bootstrapping.
#'
#' @param dat data frame or tibble containing effect sizes, variance of effect sizes, study identification number, and moderating variables.
#' @param smd name of the column containing the standardized mean differences.
#' @param var name of the column containing the variances of the standardized mean differences.
#' @param cluster name of the column containing study identifiers.
#' @param covs_full_form a string containing the right hand side of the meta-regression model (e.g., "X1 + X2 + X3 + X4 + X5").
#' @param indices indices of the variables to be tested in single coefficient test or multiple contrast hypothesis test.
#' @param R number of bootstrap replications.
#' @param adjust logical indicating whether or not to multiply residuals by CR2 adjustment matrices when bootstrapping.
#'
#'
#' @return A tibble containing the name of the test and a p-value.
#'
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' library(clubSandwich)
#' cwb(dat = SATcoaching,
#'     smd = d,
#'     var = V,
#'     cluster = study,
#'     covs_full_form = "study_type",
#'     indices = 2:3)
#'


cwb <- function(dat,
                smd,
                var,
                cluster,
                covs_full_form,
                indices,
                R = 999,
                adjust = FALSE) {

  dat$smd <- dat %>%
    dplyr::pull({{smd}})

  dat$var <- dat %>%
    dplyr::pull({{var}})

  dat$cluster <- dat %>%
    dplyr::pull({{cluster}})


  full_formula <- paste("smd ~ ", covs_full_form)

  full_model <- robumeta::robu(stats::as.formula(full_formula),
                               studynum = cluster,
                               var.eff.size = var,
                               small = FALSE,
                               data = dat)


  # Null model --------------------------------------------------------------
  terms <- c("1", unlist(stringr::str_split(covs_full_form, " \\+ ")))
  null_cov <- terms[which(!(seq_along(terms) %in% indices))]
  null_formula <- paste("smd ~ ", null_cov)

  null_model <- robumeta::robu(stats::as.formula(null_formula),
                               studynum = cluster,
                               var.eff.size = var,
                               small = FALSE,
                               data = dat)

  # residuals and transformed residuals -------------------------------------

  dat$res <- clubSandwich:::residuals_CS.robu(null_model)
  dat$pred <- with(dat, smd - res)
  split_res <- split(dat$res, dat$cluster)

  if(adjust == TRUE){
    e_tilde_j <- purrr::map(split_res, as.matrix)
    B_j <- attr(clubSandwich::vcovCR(null_model,
                                     cluster = dat$cluster,
                                     type = "CR2",
                                     inverse_var = TRUE), "adjustments")
    dat$res <- unlist(purrr::pmap(list(B_j, e_tilde_j), function(x, y) as.vector(x %*% y)))

  }


  # Bootstrapping -----------------------------------------------------------

  num_cluster <- unique(dat$cluster)
  k_j <- as.numeric(table(dat$cluster))

  system.time(

    bootstraps <- purrr::rerun(.n = R, {

      wts <- sample(c(-1, 1), size = length(num_cluster), replace = TRUE)
      dat$eta <- rep(wts, k_j)
      dat$new_y <- with(dat, pred + res * eta)

      boot_mod <- robumeta::robu(stats::as.formula(paste("new_y ~ ", covs_full_form)),
                                 studynum = cluster,
                                 var.eff.size = var,
                                 small = FALSE,
                                 dat = dat)

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
                                   vcov = clubSandwich::vcovCR(full_model, type = "CR1"),
                                   test = "Naive-F") %>%
    dplyr::pull(Fstat)


  p_boot <- bootstraps %>%
    dplyr::summarize(p_val = mean(Fstat > org_F)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(test = "CWB") %>%
    dplyr::select(test, p_val)

  if(adjust == TRUE){
    p_boot$test <- "CWB Adjusted"
  }


  return(p_boot)
}
