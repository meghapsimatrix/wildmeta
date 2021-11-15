library(robumeta)

data("SATcoaching", package = "clubSandwich")

# not dealing with missing data for the moment
SATcoaching <- subset(SATcoaching, !is.na(hrs))

# make random weights
SATcoaching$wt <- 1 + rpois(nrow(SATcoaching), lambda = 5)

check_update <- function(mod, y) {

  update_mod <- update_robu(mod, y = y)

  expect_equal(
    as.vector(mod$b.r),
    as.vector(update_mod$coefficients)
  )

  y_rand <- rnorm(length(y))
  update_rand <- update_robu(mod, y_rand)

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
}

test_that("update_robu.default works for CE models",{
  meta_CE <- robu(d ~ 1, studynum = study,
                  var.eff.size = V,
                  small = FALSE,
                  data = SATcoaching)
  mod <- meta_CE
  y <- SATcoaching$d
  res <- update_robu(mod, y = y)
  check_update(meta_CE, y = SATcoaching$d)

  reg_CE <- robu(d ~ 0 + study_type + hrs + test,
                 studynum = study,
                 var.eff.size = V,
                 small = FALSE, modelweights = "CORR",
                 data = SATcoaching)

  check_update(reg_CE, y = SATcoaching$d)
})

test_that("update_robu.default works for HE models",{

  meta_HE <- robu(d ~ 1, studynum = study,
                  var.eff.size = V,
                  small = FALSE, modelweights = "HIER",
                  data = SATcoaching)

  check_update(meta_HE, y = SATcoaching$d)

  reg_HE <- robu(d ~ 0 + study_type + hrs + test,
                  studynum = study,
                  var.eff.size = V,
                  small = FALSE, modelweights = "HIER",
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

  reg_user <- robu(d ~ 0 + study_type + hrs + test,
                    studynum = study,
                    var.eff.size = V,
                    userweights = wt,
                    small = FALSE,
                    data = SATcoaching)

  check_update(reg_user, y = SATcoaching$d)

})
