library(metafor)
library(clubSandwich)


data("oswald2013", package = "robumeta")
oswald2013$yi <- atanh(oswald2013$R)
oswald2013$vi <- 1 / (oswald2013$N - 3)
oswald2013$esID <- 1:nrow(oswald2013)

V <- impute_covariance_matrix(vi = oswald2013$vi, cluster = oswald2013$Study, r = 0.4)

mod_A <- rma.mv(yi = yi, V = V,
                random = ~ 1 | Study / esID,
                data = oswald2013,
                sparse = TRUE)

mod_B <- rma.mv(yi ~ 0 + Crit.Cat, V = V,
                random = ~ 1 | Study / esID,
                data = oswald2013,
                sparse = TRUE)

mod_C1 <- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring, V = V,
                 random = ~ 1 | Study / esID,
                 data = oswald2013,
                 sparse = TRUE)

mod_C2 <- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring, V = V,
                 random = ~ 1 | Study,
                 data = oswald2013,
                 sparse = TRUE)

mod_D <- rma.mv(yi ~ Crit.Domain + IAT.Focus + Scoring, V = V,
                random = ~ 1 | Study / esID,
                data = oswald2013,
                sparse = TRUE)

mod_E <- rma.mv(yi ~ Crit.Cat + IAT.Focus + Scoring, V = V,
                random = ~ 1 | Study,
                data = oswald2013,
                sparse = TRUE)

mod_F <- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus, V = V,
                random = ~ 1 | Study,
                data = oswald2013,
                sparse = TRUE)

Cmat_A <- constrain_equal("Crit.Cat", reg_ex = TRUE, coef(mod_B))
Cmat_B <- constrain_zero(7:10, coef(mod_C1))
Cmat_D <- constrain_zero("Crit.Cat", reg_ex = TRUE, coef(mod_C1))
Cmat_E <- constrain_zero("Crit.Domain", reg_ex = TRUE, coef(mod_C2))
Cmat_F <- constrain_zero("Scoring", reg_ex = TRUE, coef(mod_C2))

test_that("estimate_null() works for rma.mv objects.", {

  mod_A_con <- estimate_null(mod_B, C_mat = Cmat_A)
  expect_equal(get_fitted(mod_A), get_fitted(mod_A_con))
  expect_equal(get_res(mod_A), get_res(mod_A_con))

  mod_D_con <- estimate_null(mod_C1, Cmat_D)
  expect_equal(get_fitted(mod_D), get_fitted(mod_D_con))
  expect_equal(get_res(mod_D), get_res(mod_D_con))

  mod_E_con <- estimate_null(mod_C2, Cmat_E)
  expect_equal(get_fitted(mod_E), get_fitted(mod_E_con))
  expect_equal(get_res(mod_E), get_res(mod_E_con))

  mod_F_con <- estimate_null(mod_C2, Cmat_F)
  expect_equal(get_fitted(mod_F), get_fitted(mod_F_con))
  expect_equal(get_res(mod_F), get_res(mod_F_con))

  mod_B_con <- estimate_null(mod_C1, Cmat_B)
  expect_equal(get_fitted(mod_B), get_fitted(mod_B_con))
  expect_equal(get_res(mod_B), get_res(mod_B_con))

})

test_that("get_cluster() works for rma.mv objects.", {
  study_fac <- as.factor(oswald2013$Study)
  expect_equal(get_cluster(mod_A), study_fac)
  expect_equal(get_cluster(mod_B), study_fac)
  expect_equal(get_cluster(mod_C1), study_fac)
  expect_equal(get_cluster(mod_C2), study_fac)
  expect_equal(get_cluster(mod_D), study_fac)
  expect_equal(get_cluster(mod_E), study_fac)
  expect_equal(get_cluster(mod_F), study_fac)
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

  expect_equal(
    get_boot_F(mod_C1, y_boot = y_boot, C_mat = Cmat_B, cluster = get_cluster(mod_C1),
               type = "CR2", test = "EDT"),
    Wald_test(mod_C1, constraints = Cmat_B, vcov = "CR2", test = "EDT")$Fstat
  )

  expect_equal(
    get_boot_F(mod_C1, y_boot = y_boot, C_mat = Cmat_D, cluster = get_cluster(mod_C1),
               type = "CR3", test = "chi-sq"),
    Wald_test(mod_C1, constraints = Cmat_D, vcov = "CR3", test = "chi-sq")$Fstat
  )

  expect_equal(
    get_boot_F(mod_C2, y_boot = y_boot, C_mat = Cmat_E, cluster = get_cluster(mod_C2),
               type = "CR2", test = "Naive-Fp"),
    Wald_test(mod_C2, constraints = Cmat_E, vcov = "CR2", test = "Naive-Fp")$Fstat
  )

  expect_equal(
    get_boot_F(mod_C2, y_boot = y_boot, C_mat = Cmat_F, cluster = get_cluster(mod_C2),
               type = "CR0", test = "HTA"),
    Wald_test(mod_C2, constraints = Cmat_F, vcov = "CR0", test = "HTA")$Fstat
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

  boot_yi_C <- run_cwb(mod_C1, cluster = get_cluster(mod_C), R = 9,
                       auxiliary_dist = "Rademacher")
  fit_C <- get_fitted(mod_C1)
  abs_res_C <- abs(get_res(mod_C1))
  C_check <- sapply(boot_yi_C, function(x) all.equal(abs(x - fit_C), abs_res_C))
  expect_true(all(C_check))

})


test_that("Wald_test_cwb() results do not depend on sort order.", {

  ord <- sample(1:nrow(oswald2013))
  oswald_scram <- oswald2013[ord,]

})


test_that("Wald_test_cwb() works with missing values.", {

  # create missingness in predictors

  # create missingness in outcomes

  # create missingness in clusters

})


test_that("Wald_test_cwb() works with user-weighted robu models.", {

  # create missingness in predictors

  # create missingness in outcomes

  # create missingness in clusters

})

