library(robumeta)

data("SATcoaching", package = "clubSandwich")

# not dealing with missing data for the moment
SATcoaching <- subset(SATcoaching, !is.na(hrs))

# make random weights
SATcoaching$wt <- 1 + rpois(nrow(SATcoaching), lambda = 5)

mod <- robu(d ~ 1, studynum = study,
                var.eff.size = V,
                small = FALSE,
                data = SATcoaching)
y <- SATcoaching$d

check_update <- function(mod, y) {

  if (mod$small) {
    vcov_type <- "CR2"
    test_type <- "Satterthwaite"
  }  else {
    vcov_type <- "CR0"
    test_type <- "naive-tp"
  }

  update_mod <- update_robu(mod, y = y, vcov = vcov_type)
  update_tests <- clubSandwich::conf_int(update_mod,
                                         vcov = update_mod$vcov,
                                         test = test_type)

  expect_equal(
    as.vector(mod$b.r),
    as.vector(update_mod$coefficients)
  )

  expect_equal(
    mod$VR.r, as.matrix(update_mod$vcov),
    ignore_attr = TRUE
  )

  if (mod$small) {
    expect_equal(
      mod$reg_table$SE, update_tests$SE
    )
  }

  expect_equal(
    mod$reg_table$dfs, update_tests$df
  )

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
    ignore_attr = TRUE
  )

  if (mod_rand$small) {
    expect_equal(
      mod_rand$reg_table$SE, rand_tests$SE
    )
  }

  expect_equal(
    mod_rand$reg_table$dfs, rand_tests$df
  )
}

test_that("update_robu.default works for CE models",{
  meta_CE <- robu(d ~ 1, studynum = study,
                  var.eff.size = V,
                  small = FALSE,
                  data = SATcoaching)

  check_update(meta_CE, y = SATcoaching$d)

  meta_CE <- robu(d ~ 1, studynum = study,
                  var.eff.size = V,
                  small = TRUE,
                  data = SATcoaching)

  check_update(meta_CE, y = SATcoaching$d)


  reg_CE <- robu(d ~ 0 + study_type + hrs + test,
                 studynum = study,
                 var.eff.size = V,
                 small = FALSE, modelweights = "CORR",
                 data = SATcoaching)

  check_update(reg_CE, y = SATcoaching$d)

  reg_CE <- robu(d ~ 0 + study_type + hrs + test,
                 studynum = study,
                 var.eff.size = V,
                 small = TRUE, modelweights = "CORR",
                 data = SATcoaching)

  check_update(reg_CE, y = SATcoaching$d)

})

test_that("update_robu.default works for HE models",{

  meta_HE <- robu(d ~ 1, studynum = study,
                  var.eff.size = V,
                  small = FALSE, modelweights = "HIER",
                  data = SATcoaching)

  check_update(meta_HE, y = SATcoaching$d)

  meta_HE <- robu(d ~ 1, studynum = study,
                  var.eff.size = V,
                  small = TRUE, modelweights = "HIER",
                  data = SATcoaching)

  check_update(meta_HE, y = SATcoaching$d)

  reg_HE <- robu(d ~ 0 + study_type + hrs + test,
                  studynum = study,
                  var.eff.size = V,
                  small = FALSE, modelweights = "HIER",
                  data = SATcoaching)

  check_update(reg_HE, y = SATcoaching$d)

  reg_HE <- robu(d ~ 0 + study_type + hrs + test,
                 studynum = study,
                 var.eff.size = V,
                 small = TRUE, modelweights = "HIER",
                 data = SATcoaching)

  check_update(reg_HE, y = SATcoaching$d)

})

test_that("update_robu.default works for user-weighted models",{

  meta_user <- robu(d ~ 1, studynum = study,
                    var.eff.size = V,
                    userweights = wt,
                    small = FALSE,
                    data = SATcoaching)

  check_update(meta_user, y = SATcoaching$d)

  meta_user <- robu(d ~ 1, studynum = study,
                    var.eff.size = V,
                    userweights = wt,
                    small = TRUE,
                    data = SATcoaching)

  check_update(meta_user, y = SATcoaching$d)

  reg_user <- robu(d ~ 0 + study_type + hrs + test,
                    studynum = study,
                    var.eff.size = V,
                    userweights = wt,
                    small = FALSE,
                    data = SATcoaching)

  check_update(reg_user, y = SATcoaching$d)

  reg_user <- robu(d ~ 0 + study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   userweights = wt,
                   small = TRUE,
                   data = SATcoaching)

  check_update(reg_user, y = SATcoaching$d)

})
