suppressPackageStartupMessages(library(clubSandwich))
suppressPackageStartupMessages(library(robumeta))

data("oswald2013")
oswald2013$yi <- atanh(oswald2013$R)
oswald2013$vi <- 1 / (oswald2013$N - 3)

mod_A <- robu(yi ~ 1,
              studynum = Study, var.eff.size = vi,
              modelweights = "CORR",
              data = oswald2013)

mod_B <- robu(yi ~ 0 + Crit.Cat,
              studynum = Study, var.eff.size = vi,
              modelweights = "CORR",
              data = oswald2013)

mod_C <- robu(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring,
              studynum = Study, var.eff.size = vi,
              modelweights = "CORR",
              data = oswald2013)

mod_D <- robu(yi ~ Crit.Domain + IAT.Focus + Scoring,
              studynum = Study, var.eff.size = vi,
              modelweights = "CORR",
              data = oswald2013)

mod_E <- robu(yi ~ Crit.Cat + IAT.Focus + Scoring,
              studynum = Study, var.eff.size = vi,
              modelweights = "CORR",
              data = oswald2013)

mod_F <- robu(yi ~ Crit.Cat + Crit.Domain + IAT.Focus,
              studynum = Study, var.eff.size = vi,
              modelweights = "CORR",
              data = oswald2013)

Cmat_A <- constrain_equal("Crit.Cat", reg_ex = TRUE, coef(mod_B))
Cmat_B <- constrain_zero(7:10, coef(mod_C))
Cmat_D <- constrain_zero("Crit.Cat", reg_ex = TRUE, coef(mod_C))
Cmat_E <- constrain_zero("Crit.Domain", reg_ex = TRUE, coef(mod_C))
Cmat_F <- constrain_zero("Scoring", reg_ex = TRUE, coef(mod_C))

test_that("estimate_null() works for robu objects.", {

  mod_A_con <- estimate_null(mod_B, C_mat = Cmat_A)
  expect_equal(get_fitted(mod_A), get_fitted(mod_A_con))
  expect_equal(get_res(mod_A), get_res(mod_A_con))

  mod_D_con <- estimate_null(mod_C, Cmat_D)
  expect_equal(get_fitted(mod_D), get_fitted(mod_D_con))
  expect_equal(get_res(mod_D), get_res(mod_D_con))

  mod_E_con <- estimate_null(mod_C, Cmat_E)
  expect_equal(get_fitted(mod_E), get_fitted(mod_E_con))
  expect_equal(get_res(mod_E), get_res(mod_E_con))

  mod_F_con <- estimate_null(mod_C, Cmat_F)
  expect_equal(get_fitted(mod_F), get_fitted(mod_F_con))
  expect_equal(get_res(mod_F), get_res(mod_F_con))

  mod_B_con <- estimate_null(mod_C, Cmat_B)
  expect_equal(get_fitted(mod_B), get_fitted(mod_B_con))
  expect_equal(get_res(mod_B), get_res(mod_B_con))

})

test_that("get_cluster() works for robu objects.", {
  study_fac <- as.numeric(as.factor(oswald2013$Study))
  expect_equal(get_cluster(mod_A), study_fac)
  expect_equal(get_cluster(mod_B), study_fac)
  expect_equal(get_cluster(mod_C), study_fac)
  expect_equal(get_cluster(mod_D), study_fac)
  expect_equal(get_cluster(mod_E), study_fac)
  expect_equal(get_cluster(mod_F), study_fac)
})

test_that("get_boot_F_f() works for robu objects.", {

  y_boot <- oswald2013$yi

  Cmat_int <- constrain_zero(1, coefs = coef(mod_A))
  f_A <- get_boot_F_f(mod_A, C_mat = Cmat_int)

  expect_equal(
    f_A(y_boot, get_cluster(mod_A)),
    Wald_test(mod_A, constraints = Cmat_int, vcov = "CR0", test = "Naive-F")$Fstat
  )

  f_B <- get_boot_F_f(mod_B, C_mat = Cmat_A, type = "CR1", test = "HTZ")
  expect_equal(
    f_B(y_boot, get_cluster(mod_B)),
    Wald_test(mod_B, constraints = Cmat_A, vcov = "CR1", test = "HTZ")$Fstat
  )

  f_C1 <- get_boot_F_f(mod_C, C_mat = Cmat_B, type = "CR2", test = "EDT")
  expect_equal(
    f_C1(y_boot, get_cluster(mod_C)),
    Wald_test(mod_C, constraints = Cmat_B, vcov = "CR2", test = "EDT")$Fstat
  )

  f_C2 <- get_boot_F_f(mod_C, C_mat = Cmat_D, type = "CR3", test = "chi-sq")
  expect_equal(
    f_C2(y_boot, get_cluster(mod_C)),
    Wald_test(mod_C, constraints = Cmat_D, vcov = "CR3", test = "chi-sq")$Fstat
  )

  f_C3 <- get_boot_F_f(mod_C, C_mat = Cmat_E, type = "CR2", test = "Naive-Fp")
  expect_equal(
    f_C3(y_boot, get_cluster(mod_C)),
    Wald_test(mod_C, constraints = Cmat_E, vcov = "CR2", test = "Naive-Fp")$Fstat
  )

  f_C4 <- get_boot_F_f(mod_C, C_mat = Cmat_F, type = "CR0", test = "HTA")
  expect_equal(
    f_C4(y_boot, get_cluster(mod_C)),
    Wald_test(mod_C, constraints = Cmat_F, vcov = "CR0", test = "HTA")$Fstat
  )

})

test_that("run_cwb options work for robu objects.", {

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

  boot_yi_C <- run_cwb(mod_C, cluster = get_cluster(mod_C), R = 9,
                       auxiliary_dist = "Rademacher")
  fit_C <- get_fitted(mod_C)
  abs_res_C <- abs(get_res(mod_C))
  C_check <- sapply(boot_yi_C, function(x) all.equal(abs(x - fit_C), abs_res_C))
  expect_true(all(C_check))

})

test_that("Wald_test_cwb() results do not depend on sort order.", {

  ord <- sample(1:nrow(oswald2013))
  oswald_scram <- oswald2013[ord,]

  scram_A <- robu(yi ~ 1,
                  studynum = Study, var.eff.size = vi,
                  modelweights = "CORR",
                  data = oswald_scram)
  expect_equal(coef(mod_A), coef(scram_A))
  expect_equal(get_res(mod_A)[ord], get_res(scram_A))
  expect_equal(get_fitted(mod_A)[ord], get_fitted(scram_A))
  expect_equal(get_cluster(mod_A)[ord], get_cluster(scram_A))

  scram_B <- robu(yi ~ 0 + Crit.Cat,
                  studynum = Study, var.eff.size = vi,
                  modelweights = "CORR",
                  data = oswald_scram)
  expect_equal(coef(mod_B), coef(scram_B))
  expect_equal(get_res(mod_B)[ord], get_res(scram_B))
  expect_equal(get_fitted(mod_B)[ord], get_fitted(scram_B))
  expect_equal(get_cluster(mod_B)[ord], get_cluster(scram_B))

  scram_C <- robu(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring,
                  studynum = Study, var.eff.size = vi,
                  modelweights = "CORR",
                  data = oswald_scram)
  expect_equal(coef(mod_C), coef(scram_C))
  expect_equal(get_res(mod_C)[ord], get_res(scram_C))
  expect_equal(get_fitted(mod_C)[ord], get_fitted(scram_C))
  expect_equal(get_cluster(mod_C)[ord], get_cluster(scram_C))

  scram_D <- robu(yi ~ Crit.Domain + IAT.Focus + Scoring,
                  studynum = Study, var.eff.size = vi,
                  modelweights = "CORR",
                  data = oswald_scram)
  expect_equal(coef(mod_D), coef(scram_D))
  expect_equal(get_res(mod_D)[ord], get_res(scram_D))
  expect_equal(get_fitted(mod_D)[ord], get_fitted(scram_D))
  expect_equal(get_cluster(mod_D)[ord], get_cluster(scram_D))

  scram_E <- robu(yi ~ Crit.Cat + IAT.Focus + Scoring,
                  studynum = Study, var.eff.size = vi,
                  modelweights = "CORR",
                  data = oswald_scram)
  expect_equal(coef(mod_E), coef(scram_E))
  expect_equal(get_res(mod_E)[ord], get_res(scram_E))
  expect_equal(get_fitted(mod_E)[ord], get_fitted(scram_E))
  expect_equal(get_cluster(mod_E)[ord], get_cluster(scram_E))

  scram_F <- robu(yi ~ Crit.Cat + Crit.Domain + IAT.Focus,
                  studynum = Study, var.eff.size = vi,
                  modelweights = "CORR",
                  data = oswald_scram)
  expect_equal(coef(mod_F), coef(scram_F))
  expect_equal(get_res(mod_F)[ord], get_res(scram_F))
  expect_equal(get_fitted(mod_F)[ord], get_fitted(scram_F))
  expect_equal(get_cluster(mod_F)[ord], get_cluster(scram_F))

  orig_A <- Wald_test_cwb(mod_B, constraints = Cmat_A,
                               R = 19,
                               auxiliary_dist = "Rademacher",
                               adjust = "CR0",
                               type = "CR0",
                               test = "Naive-F",
                               seed = 1)

  scram_A <- Wald_test_cwb(scram_B, constraints = Cmat_A,
                           R = 19,
                           auxiliary_dist = "Rademacher",
                           adjust = "CR0",
                           type = "CR0",
                           test = "Naive-F",
                           seed = 1)


  expect_equal(attr(orig_A, "original"), attr(scram_A, "original"))
  expect_equal(attr(orig_A, "bootstraps"), attr(scram_A, "bootstraps"))


  orig_D <- Wald_test_cwb(mod_C, constraints = Cmat_D,
                          R = 12,
                          auxiliary_dist = "Mammen",
                          adjust = "CR1",
                          type = "CR1",
                          test = "Naive-Fp",
                          seed = 2)

  scram_D <- Wald_test_cwb(scram_C, constraints = Cmat_D,
                           R = 12,
                           auxiliary_dist = "Mammen",
                           adjust = "CR1",
                           type = "CR1",
                           test = "Naive-Fp",
                           seed = 2)

  expect_equal(attr(orig_D, "original"), attr(scram_D, "original"))
  expect_equal(attr(orig_D, "bootstraps"), attr(scram_D, "bootstraps"))

  orig_E <- Wald_test_cwb(mod_C, constraints = Cmat_E,
                          R = 5,
                          auxiliary_dist = "Webb six",
                          adjust = "CR2",
                          type = "CR2",
                          test = "HTZ",
                          seed = 3)

  scram_E <- Wald_test_cwb(scram_C, constraints = Cmat_E,
                           R = 5,
                           auxiliary_dist = "Webb six",
                           adjust = "CR2",
                           type = "CR2",
                           test = "HTZ",
                           seed = 3)

  expect_equal(attr(orig_E, "original"), attr(scram_E, "original"))
  expect_equal(attr(orig_E, "bootstraps"), attr(scram_E, "bootstraps"))

  orig_F <- Wald_test_cwb(mod_C, constraints = Cmat_F,
                          R = 8,
                          auxiliary_dist = "uniform",
                          adjust = "CR3",
                          type = "CR3",
                          test = "EDT",
                          seed = 4)

  scram_F <- Wald_test_cwb(scram_C, constraints = Cmat_F,
                           R = 8,
                           auxiliary_dist = "uniform",
                           adjust = "CR3",
                           type = "CR3",
                           test = "EDT",
                           seed = 4)

  expect_equal(attr(orig_F, "original"), attr(scram_F, "original"))
  expect_equal(attr(orig_F, "bootstraps"), attr(scram_F, "bootstraps"))

})

compare_mod_results <- function(mod1, mod2, tol = 1e-6) {
  expect_equal(coef(mod1), coef(mod2), tolerance = tol)
  expect_equal(get_res(mod1), get_res(mod2), tolerance = tol)
  expect_equal(get_fitted(mod1), get_fitted(mod2), tolerance = tol)
  expect_equal(get_cluster(mod1), get_cluster(mod2), tolerance = tol)
}


test_that("Wald_test_cwb() works for models with missing observations.", {


  tol <- 1e-5

  set.seed(20230212)
  oswald_missX <- oswald2013
  oswald_missX$IAT.Focus[sample(1:nrow(oswald2013), size = 40)] <- NA
  oswald_completeX <- subset(oswald_missX, !is.na(IAT.Focus))
  oswald_missY <- oswald_missX
  oswald_missY$yi[sample(1:nrow(oswald2013), size = 40)] <- NA
  oswald_completeY <- subset(oswald_missY, !is.na(yi) & !is.na(IAT.Focus))

  mod_completeY <- robu(yi ~ 0 + IAT.Focus + Crit.Domain, V = V_complete,
                       studynum = Study, var.eff.size = vi,
                       modelweights = "CORR",
                       data = oswald_completeY)

  mod_missY <- robu(yi ~ 0 + IAT.Focus + Crit.Domain, V = V_complete,
                    studynum = Study, var.eff.size = vi,
                    modelweights = "CORR",
                    data = oswald_missY)

  compare_mod_results(mod_completeY, mod_missY)

  mod_completeX <- robu(yi ~ 0 + IAT.Focus + Crit.Domain, V = V_complete,
                        studynum = Study, var.eff.size = vi,
                        modelweights = "CORR",
                        data = oswald_completeX)

  mod_missX <- robu(yi ~ 0 + IAT.Focus + Crit.Domain, V = V_complete,
                    studynum = Study, var.eff.size = vi,
                    modelweights = "CORR",
                    data = oswald_missX)

  compare_mod_results(mod_completeX, mod_missX)

  test_completeY <- Wald_test_cwb(mod_completeY,
                                 constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                                 R = 5,
                                 auxiliary_dist = "Rademacher",
                                 adjust = "CR0",
                                 type = "CR0",
                                 test = "Naive-F",
                                 seed = 19)

  test_missY <- Wald_test_cwb(mod_missY,
                            constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                            R = 5,
                            auxiliary_dist = "Rademacher",
                            adjust = "CR0",
                            type = "CR0",
                            test = "Naive-F",
                            seed = 19)

  expect_equal(attr(test_completeY, "original"), attr(test_missY, "original"), tolerance = tol)
  expect_equal(attr(test_completeY, "bootstraps"), attr(test_missY, "bootstraps"), tolerance = tol)

  test_completeX <- Wald_test_cwb(mod_completeX,
                              constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                              R = 5,
                              auxiliary_dist = "Rademacher",
                              adjust = "CR0",
                              type = "CR0",
                              test = "Naive-F",
                              seed = 19)

  test_missX <- Wald_test_cwb(mod_missX,
                              constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                              R = 5,
                              auxiliary_dist = "Rademacher",
                              adjust = "CR0",
                              type = "CR0",
                              test = "Naive-F",
                              seed = 19)

  expect_equal(attr(test_completeX, "original"), attr(test_missX, "original"), tolerance = tol)
  expect_equal(attr(test_completeX, "bootstraps"), attr(test_missX, "bootstraps"), tolerance = tol)

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

  test_missY2 <- Wald_test_cwb(mod_missY,
                             constraints = constrain_equal("IAT.Focus", reg_ex = TRUE),
                             cluster = oswald_missY$Study,
                             R = 5,
                             auxiliary_dist = "Rademacher",
                             adjust = "CR0",
                             type = "CR0",
                             test = "Naive-F",
                             seed = 19)

  expect_equal(attr(test_missY2, "original"), attr(test_missY, "original"), tolerance = tol)
  expect_equal(attr(test_missY2, "bootstraps"), attr(test_missY, "bootstraps"), tolerance = tol)


})


test_that("Wald_test_cwb() works with user-weighted robu models.", {

  oswald2013$wt <- 1 + rpois(nrow(oswald2013), lambda = 1)
  table(oswald2013$wt)

  mod_wt <- robu(yi ~ 0 + Crit.Cat + Crit.Domain + IAT.Focus + Scoring,
                 studynum = Study, var.eff.size = vi,
                 userweights = wt, modelweights = "CORR",
                 data = oswald2013)

  test_wt <- Wald_test_cwb(mod_wt,
                           constraints = constrain_equal("Crit.Cat", reg_ex = TRUE),
                           R = 19,
                           seed = 19)

  expect_s3_class(test_wt, "Wald_test_wildmeta")

})
