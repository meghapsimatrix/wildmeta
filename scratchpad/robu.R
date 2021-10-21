run_CWB.robu <- function(full_model,
                         C_mat,
                         R,
                         auxiliary_dist = "Rademacher",
                         adjust = FALSE) {

  dep <- full_model$modelweights

  full_dat <- full_model$data.full %>%
    dplyr::mutate(id = rownames(.)) %>%
    dplyr::rename(effect_size = 1,
                  v = 2)

  x_dat <- full_model$X.full %>%
    dplyr::mutate(id = rownames(.)) %>%
    dplyr::select(-1)

  dat <- full_dat %>%
    dplyr::left_join(x_dat, by = "id")

  null_formula <- paste(full_model$reg_table[, 1][ - indices], collapse = " + ")
  null_formula <- stringr::str_replace(null_formula, "X.Intercept.", "1")

  full_formula <- paste(full_model$reg_table[, 1], collapse = " + ")
  full_formula <- stringr::str_replace(full_formula, "X.Intercept.", "1")


  null_model <- robumeta::robu(stats::as.formula(paste("effect_size ~ ", null_formula)),
                               studynum = study,
                               var.eff.size = v,
                               small = FALSE,
                               modelweights = dep,
                               data = dat)

  dat$res <- clubSandwich:::residuals_CS.robu(null_model)


}
