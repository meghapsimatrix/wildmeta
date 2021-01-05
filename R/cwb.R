#' @title Calculate p-values with cluster wild bootstrapping for meta-regressions
#'
#' @description Calculates p-values for single coefficient and multiple contrast hypothesis tests using cluster wild bootstrapping.
#'
#' @param full_model a model fit using `robu()` that includes all moderators of interest.
#' @param indices indices of the variables to be tested in single coefficient test or multiple contrast hypothesis test.
#' @param R number of bootstrap replications.
#' @param adjust logical indicating whether or not to multiply residuals by CR2 adjustment matrices when bootstrapping.
#'
#'
#' @return A tibble containing the name of the test, p-value, and the name of the working model. CE denotes correlated effects and HE denotes hierarchical effects).
#'
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' library(clubSandwich)
#' library(robumeta)
#'
#' full <- robu(d ~ study_type,
#'              studynum = study,
#'              var.eff.size = V,
#'              small = FALSE,
#'              data = SATcoaching)
#'
#' cwb(full_model = full,
#'     indices = 2:3,
#'     R = 99)
#'


cwb <- function(full_model,
                indices,
                R = 999,
                adjust = FALSE) {

  dep <- full_model$modelweights

  full_dat <- full_model$data.full %>%
    dplyr::mutate(id = rownames(.))

  x_dat <- full_model$X.full %>%
    dplyr::mutate(id = rownames(.)) %>%
    dplyr::select(-1)

  dat <- full_dat %>%
    dplyr::left_join(x_dat, by = "id")

  null_formula <- paste(full_model$reg_table[, 1][ - indices], collapse = " + ")
  null_formula <- stringr::str_replace(null_formula, "X.Intercept.", "1")

  full_formula <- paste(full_model$reg_table[, 1], collapse = " + ")
  full_formula <- stringr::str_replace(full_formula, "X.Intercept.", "1")


  null_model <- robumeta::robu(stats::as.formula(paste("effect.size ~ ", null_formula)),
                               studynum = study,
                               var.eff.size = var.eff.size,
                               small = FALSE,
                               modelweights = dep,
                               data = dat)

  # residuals and transformed residuals -------------------------------------

  dat$res <- clubSandwich:::residuals_CS.robu(null_model)
  dat$pred <- with(dat, effect.size - res)
  split_res <- split(dat$res, dat$study)


  # Adjust ------------------------------------------------------------------

  if(adjust == TRUE){
    e_tilde_j <- purrr::map(split_res, as.matrix)
    B_j <- attr(clubSandwich::vcovCR(null_model,
                                     cluster = dat$study,
                                     type = "CR2",
                                     inverse_var = TRUE), "adjustments")
    dat$res <- unlist(purrr::pmap(list(B_j, e_tilde_j), function(x, y) as.vector(x %*% y)))

  }



  # bootstrap ---------------------------------------------------------------


  num_cluster <- unique(dat$study)
  k_j <- as.numeric(table(dat$study))

  system.time(

    bootstraps <- purrr::rerun(.n = R, {

      wts <- sample(c(-1, 1), size = length(num_cluster), replace = TRUE)
      dat$eta <- rep(wts, k_j)
      dat$new_y <- with(dat, pred + res * eta)

      boot_mod <- robumeta::robu(stats::as.formula(paste("new_y ~ ", full_formula)),
                                 studynum = study,
                                 var.eff.size = var.eff.size,
                                 small = FALSE,
                                 modelweights = dep,
                                 data = dat)

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

  p_boot$working_model <- dep
  p_boot$working_model <- ifelse(p_boot$working_model == "CORR", "CE",
                                 "HE")


  return(p_boot)
}
