suppressPackageStartupMessages(library(robumeta))

data("SATcoaching", package = "clubSandwich")

# not dealing with missing data for the moment
SATcoaching <- subset(SATcoaching, !is.na(hrs))

# make random weights
SATcoaching$wt <- 1 + rpois(nrow(SATcoaching), lambda = 5)

mod <- robu(d ~ 1, studynum = study,
                  var.eff.size = V,
                  userweights = wt,
                  small = TRUE,
                  data = SATcoaching)
check_dfs <- TRUE
tol <- testthat_tolerance()

check_update <- function(mod, check_dfs = TRUE, tol = testthat_tolerance()) {

  # compare update results to original model results

  y <- mod$data.full$effect.size[order(order(as.factor(mod$study_orig_id)))]

  if (mod$small) {
    vcov_type <- "CR2"
    test_type <- "Satterthwaite"
  }  else {
    vcov_type <- "CR0"
    test_type <- "naive-tp"
  }

  update_mod <- update_robu(mod, y = y, vcov = vcov_type)

  # Check bread calculations
  bread_update <- bread(update_mod)
  X_update <- model.matrix(update_mod)
  w <- update_mod$w_ij
  XwX_inv <- chol2inv(chol(crossprod(X_update, w * X_update)))
  expect_equal(bread_update, update_mod$nobs * XwX_inv)

  update_tests <- clubSandwich::conf_int(update_mod,
                                         vcov = update_mod$vcov,
                                         test = test_type)

  expect_equal(
    as.vector(mod$b.r),
    as.vector(update_mod$coefficients)
  )

  expect_equal(
    mod$VR.r, as.matrix(update_mod$vcov),
    ignore_attr = TRUE, tolerance = tol
  )

  if (mod$small) {
    expect_equal(
      mod$reg_table$SE, update_tests$SE, tolerance = tol
    )
  }

  if (check_dfs) {
    expect_equal(
      mod$reg_table$dfs, update_tests$df
    )
  }

  upup_mod <- update_robu(update_mod, y = y, vcov = vcov_type)
  upup_tests <- clubSandwich::conf_int(upup_mod,
                                         vcov = upup_mod$vcov_type,
                                         test = test_type)

  expect_equal(
    as.vector(update_mod$coefficients),
    as.vector(upup_mod$coefficients)
  )

  expect_equal(
    as.matrix(update_mod$vcov), as.matrix(upup_mod$vcov),
    ignore_attr = TRUE
  )

  expect_equal(
    update_tests$SE, upup_tests$SE
  )

  expect_equal(
    update_tests$df, upup_tests$df
  )

  # compare update to original with a random outcome

  y_rand <- rnorm(length(y))
  update_rand <- update_robu(mod, y_rand, vcov = vcov_type)
  rand_tests <- clubSandwich::conf_int(update_rand,
                                       vcov = update_rand$vcov,
                                       test = test_type)

  rand_data <- eval(mod$cl$data)
  rand_data$y_rand <- y_rand
  mod_call <- mod$cl
  mod_call$formula <- update.formula(mod_call$formula, y_rand ~ . )
  mod_call$data <- as.symbol("rand_data")
  mod_rand <- eval(mod_call)

  expect_equal(
    as.vector(mod_rand$b.r),
    as.vector(update_rand$coefficients)
  )

  expect_equal(
    mod_rand$VR.r, as.matrix(update_rand$vcov),
    ignore_attr = TRUE, tolerance = tol
  )

  if (mod_rand$small) {
    expect_equal(
      mod_rand$reg_table$SE, rand_tests$SE, tolerance = tol
    )
  }

  if (check_dfs) {
    expect_equal(
      mod_rand$reg_table$dfs, rand_tests$df
    )
  }

  upup_rand <- update_robu(update_rand, y = y_rand, vcov = vcov_type)
  upup_rand_tests <- clubSandwich::conf_int(upup_rand,
                                       vcov = upup_rand$vcov_type,
                                       test = test_type)

  expect_equal(
    as.vector(update_rand$coefficients),
    as.vector(upup_rand$coefficients)
  )

  expect_equal(
    as.matrix(update_rand$vcov), as.matrix(upup_rand$vcov),
    ignore_attr = TRUE
  )

  expect_equal(
    rand_tests$SE, upup_rand_tests$SE
  )

  expect_equal(
    rand_tests$df, upup_rand_tests$df
  )

}

test_that("update_robu.default works for CE models", {

  meta_CE <- robu(d ~ 1, studynum = study,
                  var.eff.size = V,
                  small = FALSE,
                  data = SATcoaching)

  check_update(meta_CE)

  meta_CE <- robu(d ~ 1, studynum = study,
                  var.eff.size = V,
                  small = TRUE,
                  data = SATcoaching)

  check_update(meta_CE)


  reg_CE <- robu(d ~ 0 + study_type + hrs + test,
                 studynum = study,
                 var.eff.size = V,
                 small = FALSE, modelweights = "CORR",
                 data = SATcoaching)

  check_update(reg_CE)

  reg_CE <- robu(d ~ 0 + study_type + hrs + test,
                 studynum = study,
                 var.eff.size = V,
                 small = TRUE, modelweights = "CORR",
                 data = SATcoaching)

  check_update(reg_CE)

})

test_that("update_robu.default works for HE models", {

  meta_HE <- robu(d ~ 1, studynum = study,
                  var.eff.size = V,
                  small = FALSE, modelweights = "HIER",
                  data = SATcoaching)

  check_update(meta_HE)

  meta_HE <- robu(d ~ 1, studynum = study,
                  var.eff.size = V,
                  small = TRUE, modelweights = "HIER",
                  data = SATcoaching)

  check_update(meta_HE)

  reg_HE <- robu(d ~ 0 + study_type + hrs + test,
                  studynum = study,
                  var.eff.size = V,
                  small = FALSE, modelweights = "HIER",
                  data = SATcoaching)

  check_update(reg_HE)

  reg_HE <- robu(d ~ 0 + study_type + hrs + test,
                 studynum = study,
                 var.eff.size = V,
                 small = TRUE, modelweights = "HIER",
                 data = SATcoaching)

  check_update(reg_HE)

})

test_that("update_robu.default works for user-weighted models",{

  meta_user <- robu(d ~ 1, studynum = study,
                    var.eff.size = V,
                    userweights = wt,
                    small = FALSE,
                    data = SATcoaching)

  check_update(meta_user)

  meta_user <- robu(d ~ 1, studynum = study,
                    var.eff.size = V,
                    userweights = wt,
                    small = TRUE,
                    data = SATcoaching)

  check_update(meta_user, check_dfs = FALSE)

  reg_user <- robu(d ~ 0 + study_type + hrs + test,
                    studynum = study,
                    var.eff.size = V,
                    userweights = wt,
                    small = FALSE,
                    data = SATcoaching)

  check_update(reg_user)

  reg_user <- robu(d ~ 0 + study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   userweights = wt,
                   small = TRUE,
                   data = SATcoaching)

  check_update(reg_user, check_dfs = FALSE)

})

test_that("update_robu.default doesn't works for lm objects",{


  reg_user <- lm(d ~ 0 + study_type + hrs + test,
                 weights = wt,
                 data = SATcoaching)

  expect_error(update_robu(reg_user, y = rnorm(nrow(SATcoaching))))

  reg_user$mod_label <- "linear regression"

  expect_error(update_robu(reg_user, y = rnorm(nrow(SATcoaching))))

})
