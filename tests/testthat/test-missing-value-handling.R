suppressPackageStartupMessages(library(robumeta))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(clubSandwich))

data("corrdat")

# create missingness in outcomes
corrdat_miss_y <- corrdat
missing_y <- as.logical(rbinom(nrow(corrdat), size = 1L, prob = 0.1))
corrdat_miss_y$effectsize[missing_y] <- NA
corrdat_full_y <- subset(corrdat_miss_y, !missing_y)

# create missingness in predictors
corrdat_miss_x <- corrdat
missing_x <- as.logical(rbinom(nrow(corrdat), size = 1L, prob = 0.1))
corrdat_miss_x$followup[missing_x] <- NA
corrdat_full_x <- subset(corrdat_miss_x, !missing_x)

# create missingness in clusters
corrdat_miss_cl <- corrdat
missing_cl <- as.logical(rbinom(nrow(corrdat), size = 1L, prob = 0.1))
corrdat_miss_cl$studyid[missing_cl] <- NA
corrdat_full_cl <- subset(corrdat_miss_cl, !missing_cl)


# missingness in outcomes and predictors
corrdat_miss_yx <- corrdat_miss_y
corrdat_miss_yx$followup[missing_x] <- NA
corrdat_full_yx <- subset(corrdat_miss_x, !missing_y & !missing_x)


# missingness in outcomes and clusters
corrdat_miss_yc <- corrdat_miss_y
corrdat_miss_yc$studyid[missing_cl] <- NA
corrdat_full_yc <- subset(corrdat_miss_yc, !missing_y & !missing_cl)


# missingness in predictors and clusters
corrdat_miss_xc <- corrdat_miss_x
corrdat_miss_xc$studyid[missing_cl] <- NA
corrdat_full_xc <- subset(corrdat_miss_xc, !missing_x & !missing_cl)


# missingness everywhere
corrdat$effectsize[missing_y] <- NA
corrdat$followup[missing_x] <- NA
corrdat$studyid[missing_cl] <- NA
corrdat_full <- subset(corrdat, !missing_y & !missing_x & !missing_cl)
corrdat <- corrdat



compare_robus <- function(dat_miss, dat_full, ...) {

  mod_miss <- robu(effectsize ~ binge + followup + males + college,
                   var.eff.size = var, studynum = studyid,
                   data = dat_miss,
                   modelweights = "CORR")

  mod_full <- robu(effectsize ~ binge + followup + males + college,
                   var.eff.size = var, studynum = studyid,
                   data = dat_full,
                   modelweights = "CORR")


  test_miss <- Wald_test_cwb(mod_miss, constraints = constrain_zero(2:4), ...)
  test_full <- Wald_test_cwb(mod_full, constraints = constrain_zero(2:4), ...)

  expect_equal(coef(mod_miss), coef(mod_full))
  expect_equal(attr(test_miss, "original"), attr(test_full, "original"))
  expect_equal(attr(test_miss, "bootstraps"), attr(test_full, "bootstraps"))

}


test_that("Wald_test_cwb() works with robu objects that have missing values.", {

  compare_robus(corrdat_miss_y, corrdat_full_y,
                R = 12, auxiliary_dist = "Rademacher",
                adjust = "CR0", type = "CR0",
                test = "Naive-F", seed = 11)

  compare_robus(corrdat_miss_x, corrdat_full_x,
                R = 12, auxiliary_dist = "Mammen",
                adjust = "CR2", type = "CR0",
                test = "EDT", seed = 12)

  compare_robus(corrdat_miss_cl, corrdat_full_cl,
                R = 12, auxiliary_dist = "Rademacher",
                adjust = "CR0", type = "CR1",
                test = "HTZ", seed = 13)

  compare_robus(corrdat_miss_yx, corrdat_full_yx,
                R = 12, auxiliary_dist = "Rademacher",
                adjust = "CR0", type = "CR0",
                test = "Naive-F", seed = 14)

  compare_robus(corrdat_miss_yc, corrdat_full_yc,
                R = 12, auxiliary_dist = "Rademacher",
                adjust = "CR0", type = "CR0",
                test = "EDF", seed = 15)

  compare_robus(corrdat_miss_xc, corrdat_full_xc,
                R = 12, auxiliary_dist = "Rademacher",
                adjust = "CR0", type = "CR0",
                test = "HTA", seed = 16)

  compare_robus(corrdat, corrdat_full,
                R = 12, auxiliary_dist = "Rademacher",
                adjust = "CR0", type = "CR0",
                test = "Naive-F", seed = 17)
})


compare_rmas <- function(dat_miss, dat_full, ...) {

  dat_miss <- dat_miss
  dat_full <- dat_full

  V_miss <- impute_covariance_matrix(vi = dat_miss$var,
                                     cluster = dat_miss$studyid,
                                     r = 0.7)

  suppressWarnings(
    mod_miss <- rma.mv(effectsize ~ binge + followup + males + college,
                        V = V_miss,
                        random = ~ 1 | studyid,
                        data = dat_miss)
  )

  V_full <- impute_covariance_matrix(vi = dat_full$var,
                                     cluster = dat_full$studyid,
                                     r = 0.7)

  mod_full <- rma.mv(effectsize ~ binge + followup + males + college,
                     V = V_full,
                     random = ~ 1 | studyid,
                     data = dat_full)

  suppressWarnings(
    test_miss <- Wald_test_cwb(mod_miss, constraints = constrain_zero(2:4), ...)
  )
  test_full <- Wald_test_cwb(mod_full, constraints = constrain_zero(2:4), ...)

  expect_equal(coef(mod_miss), coef(mod_full))
  expect_equal(attr(test_miss, "original"), attr(test_full, "original"))
  expect_equal(attr(test_miss, "bootstraps"), attr(test_full, "bootstraps"))

}

test_that("Wald_test_cwb() works with rma.mv objects that have missing values.", {

  compare_rmas(corrdat_miss_y, corrdat_full_y,
               R = 3, auxiliary_dist = "Rademacher",
               adjust = "CR0", type = "CR0",
               test = "Naive-F", seed = 11)

  compare_rmas(corrdat_miss_x, corrdat_full_x,
               R = 3, auxiliary_dist = "Mammen",
               adjust = "CR2", type = "CR0",
               test = "EDT", seed = 12)

  compare_rmas(corrdat_miss_yx, corrdat_full_yx,
               R = 3, auxiliary_dist = "Rademacher",
               adjust = "CR0", type = "CR0",
               test = "Naive-F", seed = 14)

})
