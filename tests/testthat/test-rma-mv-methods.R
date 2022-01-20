suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(clubSandwich))


data("oswald2013", package = "robumeta")
oswald2013$yi <- atanh(oswald2013$R)
oswald2013$vi <- 1 / (oswald2013$N - 3)
oswald2013$esID <- 1:nrow(oswald2013)
oswald2013$wt <- 1 + rpois(nrow(oswald2013), lambda = 1)
table(oswald2013$wt)
oswald2013 <<- oswald2013

V <<- impute_covariance_matrix(vi = oswald2013$vi, cluster = oswald2013$Study, r = 0.4)

mod_A <<- rma.mv(yi = yi, V = V,
                random = ~ 1 | Study / esID,
                data = oswald2013,
                sparse = TRUE)

mod_B <<- rma.mv(yi ~ 0 + Crit.Cat, V = V,
                random = ~ 1 | Study / esID,
                data = oswald2013,
                sparse = TRUE)

mod_C1 <<- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring, V = V,
                 random = ~ 1 | Study / esID,
                 data = oswald2013,
                 sparse = TRUE)

mod_C2 <<- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring, V = V,
                 random = ~ 1 | Study,
                 data = oswald2013,
                 sparse = TRUE)

mod_D <<- rma.mv(yi ~ Crit.Domain + IAT.Focus + Scoring, V = V,
                random = ~ 1 | Study / esID,
                data = oswald2013,
                sparse = TRUE)

mod_E <<- rma.mv(yi ~ Crit.Cat + IAT.Focus + Scoring, V = V,
                random = ~ 1 | Study,
                data = oswald2013,
                sparse = TRUE)

mod_F <<- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus, V = V,
                random = ~ 1 | Study,
                data = oswald2013,
                sparse = TRUE)

Cmat_A <<- constrain_equal("Crit.Cat", reg_ex = TRUE, coef(mod_B))
Cmat_B <<- constrain_zero(7:10, coef(mod_C1))
Cmat_D <<- constrain_zero("Crit.Cat", reg_ex = TRUE, coef(mod_C1))
Cmat_E <<- constrain_zero("Crit.Domain", reg_ex = TRUE, coef(mod_C2))
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

  boot_yi_C <- run_cwb(mod_C1, cluster = get_cluster(mod_C1), R = 9,
                       auxiliary_dist = "Rademacher")
  fit_C <- get_fitted(mod_C1)
  abs_res_C <- abs(get_res(mod_C1))
  C_check <- sapply(boot_yi_C, function(x) all.equal(abs(x - fit_C), abs_res_C))
  expect_true(all(C_check))

})

compare_mod_results <- function(mod, scram, ord = 1:length(get_res(mod)), tol = 1e-6) {
  expect_equal(coef(mod), coef(scram), tolerance = tol)
  expect_equal(as.numeric(get_res(mod)[ord]), as.numeric(get_res(scram)), tolerance = tol)
  expect_equal(as.numeric(get_fitted(mod)[ord]), as.numeric(get_fitted(scram)), tolerance = tol)
  expect_equal(get_cluster(mod)[ord], get_cluster(scram))
}

test_that("Wald_test_cwb() results do not depend on sort order.", {

  skip_on_cran()

  ord <- sample(1:nrow(oswald2013))
  oswald_scram <<- oswald2013[ord,]

  V_scram <<- impute_covariance_matrix(vi = oswald_scram$vi, cluster = oswald_scram$Study, r = 0.4)

  scram_A <<- rma.mv(yi = yi, V = V_scram,
                    random = ~ 1 | Study / esID,
                    data = oswald_scram,
                    sparse = TRUE)
  compare_mod_results(mod_A, scram_A, ord)

  scram_B <<- rma.mv(yi ~ 0 + Crit.Cat, V = V_scram,
                    random = ~ 1 | Study / esID,
                    data = oswald_scram,
                    sparse = TRUE)
  compare_mod_results(mod_B, scram_B, ord)

  scram_C1 <<- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring, V = V_scram,
                     random = ~ 1 | Study / esID,
                     data = oswald_scram,
                     sparse = TRUE)
  compare_mod_results(mod_C1, scram_C1, ord)

  scram_C2 <<- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring, V = V_scram,
                     random = ~ 1 | Study,
                     data = oswald_scram,
                     sparse = TRUE)
  compare_mod_results(mod_C2, scram_C2, ord)

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


  orig_D <- Wald_test_cwb(mod_C1, constraints = Cmat_D,
                          R = 3,
                          auxiliary_dist = "Mammen",
                          adjust = "CR1",
                          type = "CR1",
                          test = "Naive-Fp",
                          seed = 2)

  scram_D <- Wald_test_cwb(scram_C1, constraints = Cmat_D,
                           R = 3,
                           auxiliary_dist = "Mammen",
                           adjust = "CR1",
                           type = "CR1",
                           test = "Naive-Fp",
                           seed = 2)

  expect_equal(attr(orig_D, "original"), attr(scram_D, "original"))
  expect_equal(attr(orig_D, "bootstraps"), attr(scram_D, "bootstraps"))

  orig_E <- Wald_test_cwb(mod_C2, constraints = Cmat_E,
                          R = 5,
                          auxiliary_dist = "Webb six",
                          adjust = "CR2",
                          type = "CR2",
                          test = "HTZ",
                          seed = 3)

  scram_E <- Wald_test_cwb(scram_C2, constraints = Cmat_E,
                           R = 5,
                           auxiliary_dist = "Webb six",
                           adjust = "CR2",
                           type = "CR2",
                           test = "HTZ",
                           seed = 3)

  expect_equal(attr(orig_E, "original"), attr(scram_E, "original"))
  expect_equal(attr(orig_E, "bootstraps"), attr(scram_E, "bootstraps"))

  orig_F <- Wald_test_cwb(mod_C2, constraints = Cmat_F,
                          R = 8,
                          auxiliary_dist = "uniform",
                          adjust = "CR3",
                          type = "CR3",
                          test = "EDT",
                          seed = 4)

  scram_F <- Wald_test_cwb(scram_C2, constraints = Cmat_F,
                           R = 8,
                           auxiliary_dist = "uniform",
                           adjust = "CR3",
                           type = "CR3",
                           test = "EDT",
                           seed = 4)

  expect_equal(attr(orig_F, "original"), attr(scram_F, "original"))
  expect_equal(attr(orig_F, "bootstraps"), attr(scram_F, "bootstraps"))

})


test_that("Wald_test_cwb() works when rma.mv uses subset.", {

  oswald_sub <<- subset(oswald2013, Crit.Cat == "Microbehavior")

  V_sub <<- impute_covariance_matrix(vi = oswald_sub$vi, cluster = oswald_sub$Study, r = 0.4)

  mod_full <<- rma.mv(yi = yi, V = V_sub,
                     mods = ~ 0 + IAT.Focus + Crit.Domain,
                     random = ~ 1 | Study / esID,
                     data = oswald_sub,
                     sparse = TRUE)

  mod_sub <<- rma.mv(yi = yi, V = V,
                    mods = ~ 0 + IAT.Focus + Crit.Domain,
                    random = ~ 1 | Study / esID,
                    data = oswald2013,
                    subset = Crit.Cat == "Microbehavior",
                    sparse = TRUE)

  compare_mod_results(mod_full, mod_sub)

  test_full <- Wald_test_cwb(mod_full,
                             constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                             R = 4,
                             auxiliary_dist = "Rademacher",
                             adjust = "CR0",
                             type = "CR0",
                             test = "Naive-F",
                             seed = 19)

  test_sub <- Wald_test_cwb(mod_sub,
                            constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                            R = 4,
                            auxiliary_dist = "Rademacher",
                            adjust = "CR0",
                            type = "CR0",
                            test = "Naive-F",
                            seed = 19)

  expect_equal(attr(test_full, "original"), attr(test_sub, "original"))
  expect_equal(attr(test_full, "bootstraps"), attr(test_sub, "bootstraps"))

})



test_that("Wald_test_cwb() works with user-weighted rma.mv models.", {

  mod_wt <<- rma.mv(yi ~ 0 + Crit.Cat + Crit.Domain + IAT.Focus + Scoring,
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

