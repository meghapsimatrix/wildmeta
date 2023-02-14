suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(clubSandwich))

set.seed(20220126)
data("oswald2013", package = "robumeta")

# Create variables
oswald2013$yi <- atanh(oswald2013$R)
oswald2013$vi <- 1 / (oswald2013$N - 3)
oswald2013$esID <- 1:nrow(oswald2013)
oswald2013$wt <- 1 + rpois(nrow(oswald2013), lambda = 1)
table(oswald2013$wt)
oswald2013 <- oswald2013

# Use a subset of 30 studies
oswald2013_full <- oswald2013
study_sample <- sample(unique(oswald2013$Study), size = 30)
oswald2013 <- subset(oswald2013, Study %in% study_sample)


# Fit some models
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

mod_Y <- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring, V = V,
                 data = oswald2013,
                 sparse = TRUE)

mod_Z <- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus, V = V,
                data = oswald2013,
                sparse = TRUE)

Cmat_A <- constrain_equal("Crit.Cat", reg_ex = TRUE, coef(mod_B))
Cmat_B <- constrain_zero(7:10, coef(mod_C1))
Cmat_D <- constrain_zero("Crit.Cat", reg_ex = TRUE, coef(mod_C1))
Cmat_E <- constrain_zero("Crit.Domain", reg_ex = TRUE, coef(mod_C2))
Cmat_F <- constrain_zero("Scoring", reg_ex = TRUE, coef(mod_C2))

tol <- 1e-6

test_that("estimate_null() works for rma.mv objects.", {

  mod_A_con <- estimate_null(mod_B, C_mat = Cmat_A)
  expect_equal(get_fitted(mod_A), get_fitted(mod_A_con), tolerance = tol)
  expect_equal(get_res(mod_A), get_res(mod_A_con), tolerance = tol)

  mod_D_con <- estimate_null(mod_C1, Cmat_D)
  expect_equal(get_fitted(mod_D), get_fitted(mod_D_con), tolerance = tol)
  expect_equal(get_res(mod_D), get_res(mod_D_con), tolerance = tol)

  mod_E_con <- estimate_null(mod_C2, Cmat_E)
  expect_equal(get_fitted(mod_E), get_fitted(mod_E_con), tolerance = tol)
  expect_equal(get_res(mod_E), get_res(mod_E_con), tolerance = tol)

  mod_F_con <- estimate_null(mod_C2, Cmat_F)
  expect_equal(get_fitted(mod_F), get_fitted(mod_F_con), tolerance = tol)
  expect_equal(get_res(mod_F), get_res(mod_F_con), tolerance = tol)

  mod_B_con <- estimate_null(mod_C1, Cmat_B)
  expect_equal(get_fitted(mod_B), get_fitted(mod_B_con), tolerance = tol)
  expect_equal(get_res(mod_B), get_res(mod_B_con), tolerance = tol)

  mod_Z_con <- estimate_null(mod_Y, Cmat_F)
  expect_equal(get_fitted(mod_Z), get_fitted(mod_Z_con), tolerance = tol)
  expect_equal(get_res(mod_Z), get_res(mod_Z_con), tolerance = tol)

  oswald2013_mod <- oswald2013
  oswald2013_mod$X_null <- 2

  mod_G <- update(mod_C2, data = oswald2013_mod)
  mod_F_Xnull <- estimate_null(mod_G, Cmat_F)
  expect_equal(get_fitted(mod_F), get_fitted(mod_F_Xnull), tolerance = tol)
  expect_equal(get_res(mod_F), get_res(mod_F_Xnull), tolerance = tol)

})

test_that("get_cluster() works for rma.mv objects.", {
  study_fac <- as.factor(oswald2013$Study)
  expect_equal(get_cluster(mod_A), study_fac, tolerance = tol)
  expect_equal(get_cluster(mod_B), study_fac, tolerance = tol)
  expect_equal(get_cluster(mod_C1), study_fac, tolerance = tol)
  expect_equal(get_cluster(mod_C2), study_fac, tolerance = tol)
  expect_equal(get_cluster(mod_D), study_fac, tolerance = tol)
  expect_equal(get_cluster(mod_E), study_fac, tolerance = tol)
  expect_equal(get_cluster(mod_F), study_fac, tolerance = tol)
})

test_that("get_boot_F_F() works for rma.mv objects.", {

  y_boot <- oswald2013$yi

  Cmat_int <- constrain_zero(1, coefs = coef(mod_A))

  f_A <- get_boot_F_f(mod_A, C_mat = Cmat_int)
  expect_equal(
    f_A(y_boot = y_boot, cluster = get_cluster(mod_A)),
    Wald_test(mod_A, constraints = Cmat_int, vcov = "CR0", test = "Naive-F")$Fstat, tolerance = tol
  )

  f_B <- get_boot_F_f(mod_B, C_mat = Cmat_A, type = "CR1", test = "HTZ")
  expect_equal(
    f_B(y_boot = y_boot, cluster = get_cluster(mod_B)),
    Wald_test(mod_B, constraints = Cmat_A, vcov = "CR1", test = "HTZ")$Fstat, tolerance = tol
  )

  f_C1 <- get_boot_F_f(mod_C1, C_mat = Cmat_B, type = "CR2", test = "EDT")
  expect_equal(
    f_C1(y_boot = y_boot, cluster = get_cluster(mod_C1)),
    Wald_test(mod_C1, constraints = Cmat_B, vcov = "CR2", test = "EDT")$Fstat, tolerance = tol
  )

  f_C1 <- get_boot_F_f(mod_C1, C_mat = Cmat_D, type = "CR3", test = "chi-sq")
  expect_equal(
    f_C1(y_boot = y_boot, cluster = get_cluster(mod_C1)),
    Wald_test(mod_C1, constraints = Cmat_D, vcov = "CR3", test = "chi-sq")$Fstat, tolerance = tol
  )

  f_C2 <- get_boot_F_f(mod_C2, C_mat = Cmat_E, type = "CR2", test = "Naive-Fp")
  expect_equal(
    f_C2(y_boot = y_boot, cluster = get_cluster(mod_C2)),
    Wald_test(mod_C2, constraints = Cmat_E, vcov = "CR2", test = "Naive-Fp")$Fstat, tolerance = tol
  )

  f_C2 <- get_boot_F_f(mod_C2, C_mat = Cmat_F, type = "CR0", test = "HTA")
  expect_equal(
    f_C2(y_boot = y_boot, cluster = get_cluster(mod_C2)),
    Wald_test(mod_C2, constraints = Cmat_F, vcov = "CR0", test = "HTA")$Fstat, tolerance = tol
  )

  f_Y <- get_boot_F_f(mod_Y, C_mat = Cmat_F, type = "CR0", test = "HTA")
  expect_equal(
    f_Y(y_boot = y_boot, cluster = oswald2013$Study),
    Wald_test(mod_Y, constraints = Cmat_F, cluster = oswald2013$Study,
              vcov = "CR0", test = "HTA")$Fstat, tolerance = tol
  )


  oswald2013$yi_boot <- oswald2013$yi + rnorm(nrow(oswald2013))

  mod_G <- rma.mv(yi = yi_boot, mods = ~ Crit.Cat + Crit.Domain + IAT.Focus, V = V,
                  random = ~ 1 | Study,
                  data = oswald2013,
                  sparse = TRUE)
  mod_H <- rma.mv(yi = yi_boot ~ Crit.Cat + Crit.Domain + IAT.Focus, V = V,
                  random = ~ 1 | Study,
                  data = oswald2013,
                  sparse = TRUE)
  expect_false(all(coef(mod_F) == coef(mod_G)))
  expect_equal(coef(mod_G), coef(mod_H), tolerance = tol)
  expect_false(all(y_boot == oswald2013$yi_boot))

  Cmat_G <- constrain_equal("Crit.Cat", reg_ex = TRUE, coefs = coef(mod_G))
  F_stat <- Wald_test(mod_F, constraints = Cmat_G, vcov = "CR0", test = "HTZ")$Fstat
  G_stat <- get_boot_F_f(mod_G, C_mat = Cmat_G, type = "CR0", test = "HTZ")(y_boot = y_boot, cluster = get_cluster(mod_G))
  H_stat <- get_boot_F_f(mod_H, C_mat = Cmat_G, type = "CR0", test = "HTZ")(y_boot = y_boot, cluster = get_cluster(mod_G))

  expect_equal(G_stat, F_stat, tolerance = tol)

  expect_equal(H_stat, F_stat, tolerance = tol)

  expect_false(G_stat == Wald_test(mod_G, constraints = Cmat_G, vcov = "CR0", test = "HTZ")$Fstat)

})


test_that("run_cwb options work for rma.mv objects.", {

  boot_yi_A <- run_cwb(mod_A, cluster = get_cluster(mod_A), R = 3,
                       auxiliary_dist = "Rademacher")
  fit_A <- get_fitted(mod_A)
  abs_res_A <- abs(get_res(mod_A))
  A_check <- sapply(boot_yi_A, function(x) all.equal(abs(x - fit_A), abs_res_A))
  expect_true(all(A_check))

  boot_yi_B <- run_cwb(mod_B, cluster = get_cluster(mod_B), R = 3,
                       auxiliary_dist = "Rademacher")
  fit_B <- get_fitted(mod_B)
  abs_res_B <- abs(get_res(mod_B))
  B_check <- sapply(boot_yi_B, function(x) all.equal(abs(x - fit_B), abs_res_B))
  expect_true(all(B_check))

  boot_yi_C <- run_cwb(mod_C1, cluster = get_cluster(mod_C1), R = 3,
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
  expect_equal(get_cluster(mod)[ord], get_cluster(scram), tolerance = tol)
}


test_that("Wald_test_cwb() results do not depend on sort order.", {

  skip_on_cran()

  tol <- 1e-5

  ord <- sample(1:nrow(oswald2013))
  oswald_scram <- oswald2013[ord,]

  V_scram <- impute_covariance_matrix(vi = oswald_scram$vi, cluster = oswald_scram$Study, r = 0.4)


  scram_A <- rma.mv(yi = yi, V = V_scram,
                    random = ~ 1 | Study / esID,
                    data = oswald_scram,
                    sparse = TRUE)
  compare_mod_results(mod_A, scram_A, ord)

  scram_B <- rma.mv(yi ~ 0 + Crit.Cat, V = V_scram,
                    random = ~ 1 | Study / esID,
                    data = oswald_scram,
                    sparse = TRUE)
  compare_mod_results(mod_B, scram_B, ord)

  scram_C1 <- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring, V = V_scram,
                     random = ~ 1 | Study / esID,
                     data = oswald_scram,
                     sparse = TRUE)
  compare_mod_results(mod_C1, scram_C1, ord)

  scram_C2 <- rma.mv(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring, V = V_scram,
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

  expect_equal(attr(orig_A, "original"), attr(scram_A, "original"), tolerance = tol)
  expect_equal(attr(orig_A, "bootstraps"), attr(scram_A, "bootstraps"), tolerance = tol)


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

  expect_equal(attr(orig_D, "original"), attr(scram_D, "original"), tolerance = tol)
  expect_equal(attr(orig_D, "bootstraps"), attr(scram_D, "bootstraps"), tolerance = tol)

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

  expect_equal(attr(orig_E, "original"), attr(scram_E, "original"), tolerance = tol)
  expect_equal(attr(orig_E, "bootstraps"), attr(scram_E, "bootstraps"), tolerance = tol)

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

  expect_equal(attr(orig_F, "original"), attr(scram_F, "original"), tolerance = tol)
  expect_equal(attr(orig_F, "bootstraps"), attr(scram_F, "bootstraps"), tolerance = tol)

})


test_that("Wald_test_cwb() works when rma.mv uses subset.", {

  set.seed(20230212)
  oswald_missX <- oswald2013_full
  oswald_missX$IAT.Focus[sample(1:nrow(oswald2013_full), size = 40)] <- NA
  oswald_sub <- subset(oswald_missX, Crit.Cat == "Microbehavior")
  oswald_complete <- subset(oswald_missX, Crit.Cat == "Microbehavior" & !is.na(IAT.Focus))

  V_missX <- impute_covariance_matrix(vi = oswald_missX$vi, cluster = oswald_missX$Study, r = 0.4)
  V_sub <- impute_covariance_matrix(vi = oswald_sub$vi, cluster = oswald_sub$Study, r = 0.4)
  V_complete <- impute_covariance_matrix(vi = oswald_complete$vi, cluster = oswald_complete$Study, r = 0.4)

  mod_complete <- rma.mv(yi = yi, V = V_complete,
                         mods = ~ 0 + IAT.Focus + Crit.Domain,
                         random = ~ 1 | Study / esID,
                         data = oswald_complete,
                         sparse = TRUE)

  expect_warning(
    mod_sub <- rma.mv(yi = yi, V = V_sub,
                      mods = ~ 0 + IAT.Focus + Crit.Domain,
                      random = ~ 1 | Study / esID,
                      data = oswald_sub,
                      sparse = TRUE)
  )

  expect_warning(
    mod_missX <- rma.mv(yi = yi, V = V_missX,
                        mods = ~ 0 + IAT.Focus + Crit.Domain,
                        random = ~ 1 | Study / esID,
                        data = oswald_missX,
                        subset = Crit.Cat == "Microbehavior",
                        sparse = TRUE)
  )

  compare_mod_results(mod_complete, mod_sub)
  compare_mod_results(mod_complete, mod_missX)

  test_complete <- Wald_test_cwb(mod_sub,
                                 constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                                 R = 5,
                                 auxiliary_dist = "Rademacher",
                                 adjust = "CR0",
                                 type = "CR0",
                                 test = "Naive-F",
                                 seed = 19)

  test_sub <- Wald_test_cwb(mod_sub,
                            constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                            R = 5,
                            auxiliary_dist = "Rademacher",
                            adjust = "CR0",
                            type = "CR0",
                            test = "Naive-F",
                            seed = 19)

  expect_equal(attr(test_complete, "original"), attr(test_sub, "original"), tolerance = tol)
  expect_equal(attr(test_complete, "bootstraps"), attr(test_sub, "bootstraps"), tolerance = tol)

  test_missX <- Wald_test_cwb(mod_missX,
                              constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                              R = 5,
                              auxiliary_dist = "Rademacher",
                              adjust = "CR0",
                              type = "CR0",
                              test = "Naive-F",
                              seed = 19)

  expect_equal(attr(test_complete, "original"), attr(test_missX, "original"), tolerance = tol)
  expect_equal(attr(test_complete, "bootstraps"), attr(test_missX, "bootstraps"), tolerance = tol)

  test_missX2 <- Wald_test_cwb(mod_missX,
                               constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                               cluster = oswald_missX$Study,
                               R = 5,
                               auxiliary_dist = "Rademacher",
                               adjust = "CR0",
                               type = "CR0",
                               test = "Naive-F",
                               seed = 19)

  expect_equal(attr(test_missX2, "original"), attr(test_missX, "original"), tolerance = tol)
  expect_equal(attr(test_missX2, "bootstraps"), attr(test_missX, "bootstraps"), tolerance = tol)

  test_sub2 <- Wald_test_cwb(mod_sub,
                             constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                             cluster = oswald_sub$Study,
                             R = 5,
                             auxiliary_dist = "Rademacher",
                             adjust = "CR0",
                             type = "CR0",
                             test = "Naive-F",
                             seed = 19)

  expect_equal(attr(test_sub2, "original"), attr(test_sub, "original"), tolerance = tol)
  expect_equal(attr(test_sub2, "bootstraps"), attr(test_sub, "bootstraps"), tolerance = tol)


})



test_that("Wald_test_cwb() works with user-weighted rma.mv models.", {

  skip_if(packageVersion("clubSandwich") < '0.5.5.9999')

  W_mat <- impute_covariance_matrix(vi = oswald2013$wt, cluster = oswald2013$Study, r = 0, return_list = FALSE)

  mod_wt1 <- rma.mv(yi ~ 0 + Crit.Cat + Crit.Domain + IAT.Focus + Scoring,
                     V = V, W = W_mat,
                     random = ~ 1 | Study,
                     data = oswald2013,
                     sparse = TRUE)


  test_wt1 <- Wald_test_cwb(mod_wt1,
                            constraints = constrain_equal("Crit.Cat", reg_ex = TRUE),
                            R = 3,
                            seed = 19)

  expect_s3_class(test_wt1, "Wald_test_wildmeta")
  expect_true(!is.na(test_wt1$p_val))

  mod_wt2 <- rma.mv(yi ~ 0 + Crit.Cat + Crit.Domain + IAT.Focus + Scoring,
                     V = V, W = wt,
                     random = ~ 1 | Study,
                     data = oswald2013,
                     sparse = TRUE)

  test_wt2 <- Wald_test_cwb(mod_wt2,
                            constraints = constrain_equal("Crit.Cat", reg_ex = TRUE),
                            R = 3,
                            seed = 19)

  expect_s3_class(test_wt2, "Wald_test_wildmeta")
  expect_true(!is.na(test_wt2$p_val))

  expect_equal(test_wt1, test_wt2, tolerance = tol)

})

