suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(clubSandwich))


data("oswald2013", package = "robumeta")
oswald2013$yi <- atanh(oswald2013$R)
oswald2013$vi <- 1 / (oswald2013$N - 3)
oswald2013$esID <- 1:nrow(oswald2013)
oswald2013$wt <- 1 + rpois(nrow(oswald2013), lambda = 1)
table(oswald2013$wt)


V <- impute_covariance_matrix(vi = oswald2013$vi, cluster = oswald2013$Study, r = 0.4)

mod_A <- rma.mv(yi = yi, V = V,
                random = ~ 1 | Study,
                data = oswald2013,
                sparse = TRUE)

mod_B <- rma.mv(yi ~ 0 + Crit.Cat, V = V,
                random = ~ 1 | Study,
                data = oswald2013,
                sparse = TRUE)


Cmat_A <- constrain_equal("Crit.Cat", reg_ex = TRUE, coef(mod_B))

test_that("estimate_null() works for rma.mv objects.", {

  mod_A_con <- estimate_null(mod_B, C_mat = Cmat_A)
  expect_equal(get_fitted(mod_A), get_fitted(mod_A_con))
  expect_equal(get_res(mod_A), get_res(mod_A_con))

})

test_that("get_cluster() works for rma.mv objects.", {
  study_fac <- as.factor(oswald2013$Study)
  expect_equal(get_cluster(mod_A), study_fac)
  expect_equal(get_cluster(mod_B), study_fac)
})

test_that("get_boot_F() works for rma.mv objects.", {

  y_boot <- oswald2013$yi

  Cmat_int <- constrain_zero(1, coefs = coef(mod_A))
  expect_equal(
    get_boot_F(mod_A, y_boot = y_boot, C_mat = Cmat_int, cluster = get_cluster(mod_A)),
    Wald_test(mod_A, constraints = Cmat_int, vcov = "CR0", test = "Naive-F")$Fstat
  )

  expect_equal(
    get_boot_F(mod_B, y_boot = y_boot, C_mat = Cmat_A, cluster = get_cluster(mod_B),
               type = "CR1", test = "HTZ"),
    Wald_test(mod_B, constraints = Cmat_A, vcov = "CR1", test = "HTZ")$Fstat
  )


})

test_that("run_cwb options work for rma.mv objects.", {

  boot_yi_A <- run_cwb(mod_A, cluster = get_cluster(mod_A), R = 9,
                       auxiliary_dist = "Rademacher")
  fit_A <- get_fitted(mod_A)
  abs_res_A <- abs(get_res(mod_A))
  A_check <- sapply(boot_yi_A, function(x) all.equal(abs(x - fit_A), abs_res_A))
  expect_true(all(A_check))

  boot_yi_B <- run_cwb(mod_B, cluster = get_cluster(mod_B), R = 9,
                       auxiliary_dist = "Rademacher")
  fit_B <- get_fitted(mod_B)
  abs_res_B <- abs(get_res(mod_B))
  B_check <- sapply(boot_yi_B, function(x) all.equal(abs(x - fit_B), abs_res_B))
  expect_true(all(B_check))

})


test_that("Wald_test_cwb() results do not depend on sort order.", {

  skip_on_cran()

  ord <- sample(1:nrow(oswald2013))
  oswald_scram <- oswald2013[ord,]

  V_scram <- impute_covariance_matrix(vi = oswald_scram$vi, cluster = oswald_scram$Study, r = 0.4)

  scram_A <- rma.mv(yi = yi, V = V_scram,
                    random = ~ 1 | Study,
                    data = oswald_scram,
                    sparse = TRUE)
  expect_equal(coef(mod_A), coef(scram_A))
  expect_equal(as.numeric(get_res(mod_A)[ord]), as.numeric(get_res(scram_A)))
  expect_equal(as.numeric(get_fitted(mod_A)[ord]), as.numeric(get_fitted(scram_A)))
  expect_equal(get_cluster(mod_A)[ord], get_cluster(scram_A))

  scram_B <- rma.mv(yi ~ 0 + Crit.Cat, V = V_scram,
                    random = ~ 1 | Study,
                    data = oswald_scram,
                    sparse = TRUE)
  expect_equal(coef(mod_B), coef(scram_B))
  expect_equal(as.numeric(get_res(mod_B)[ord]), as.numeric(get_res(scram_B)))
  expect_equal(as.numeric(get_fitted(mod_B)[ord]), as.numeric(get_fitted(scram_B)))
  expect_equal(get_cluster(mod_B)[ord], get_cluster(scram_B))


  orig_A <- Wald_test_cwb(mod_B, constraints = Cmat_A,
                          R = 4,
                          auxiliary_dist = "Rademacher",
                          adjust = "CR0",
                          type = "CR0",
                          test = "Naive-F",
                          seed = 1)

  scram_A <- Wald_test_cwb(scram_B, constraints = Cmat_A,
                           R = 4,
                           auxiliary_dist = "Rademacher",
                           adjust = "CR0",
                           type = "CR0",
                           test = "Naive-F",
                           seed = 1)

  expect_equal(attr(orig_A, "original"), attr(scram_A, "original"))
  expect_equal(attr(orig_A, "bootstraps"), attr(scram_A, "bootstraps"))


})



test_that("Wald_test_cwb() works with user-weighted rma.mv models.", {

  mod_wt <- rma.mv(yi ~ 0 + Crit.Cat + Crit.Domain + IAT.Focus + Scoring,
                   V = V, W = wt,
                   random = ~ 1 | Study,
                   data = oswald2013,
                   sparse = TRUE)

  test_wt <- Wald_test_cwb(mod_wt,
                           constraints = constrain_equal("Crit.Cat", reg_ex = TRUE),
                           R = 3,
                           seed = 19)

  expect_s3_class(test_wt, "Wald_test_wildmeta")

})

