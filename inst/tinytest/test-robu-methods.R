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


# estimate_null() works for robu objects.

mod_A_con <- estimate_null(mod_B, C_mat = Cmat_A)
tinytest::expect_equal(get_fitted(mod_A), get_fitted(mod_A_con))
tinytest::expect_equal(get_res(mod_A), get_res(mod_A_con))

mod_D_con <- estimate_null(mod_C, Cmat_D)
tinytest::expect_equal(get_fitted(mod_D), get_fitted(mod_D_con))
tinytest::expect_equal(get_res(mod_D), get_res(mod_D_con))

mod_E_con <- estimate_null(mod_C, Cmat_E)
tinytest::expect_equal(get_fitted(mod_E), get_fitted(mod_E_con))
tinytest::expect_equal(get_res(mod_E), get_res(mod_E_con))

mod_F_con <- estimate_null(mod_C, Cmat_F)
tinytest::expect_equal(get_fitted(mod_F), get_fitted(mod_F_con))
tinytest::expect_equal(get_res(mod_F), get_res(mod_F_con))

mod_B_con <- estimate_null(mod_C, Cmat_B)
tinytest::expect_equal(get_fitted(mod_B), get_fitted(mod_B_con))
tinytest::expect_equal(get_res(mod_B), get_res(mod_B_con))


# get_cluster() works for robu objects.

study_fac <- as.numeric(as.factor(oswald2013$Study))
tinytest::expect_equal(get_cluster(mod_A), study_fac)
tinytest::expect_equal(get_cluster(mod_B), study_fac)
tinytest::expect_equal(get_cluster(mod_C), study_fac)
tinytest::expect_equal(get_cluster(mod_D), study_fac)
tinytest::expect_equal(get_cluster(mod_E), study_fac)
tinytest::expect_equal(get_cluster(mod_F), study_fac)


# get_boot_F() works for robu objects.

y_boot <- oswald2013$yi

Cmat_int <- constrain_zero(1, coefs = coef(mod_A))
tinytest::expect_equal(
  get_boot_F(mod_A, y_boot = y_boot, C_mat = Cmat_int, cluster = get_cluster(mod_A)),
  Wald_test(mod_A, constraints = Cmat_int, vcov = "CR0", test = "Naive-F")$Fstat
)

tinytest::expect_equal(
  get_boot_F(mod_B, y_boot = y_boot, C_mat = Cmat_A, cluster = get_cluster(mod_B),
             type = "CR1", test = "HTZ"),
  Wald_test(mod_B, constraints = Cmat_A, vcov = "CR1", test = "HTZ")$Fstat
)

tinytest::expect_equal(
  get_boot_F(mod_C, y_boot = y_boot, C_mat = Cmat_B, cluster = get_cluster(mod_C),
             type = "CR2", test = "EDT"),
  Wald_test(mod_C, constraints = Cmat_B, vcov = "CR2", test = "EDT")$Fstat
)

tinytest::expect_equal(
  get_boot_F(mod_C, y_boot = y_boot, C_mat = Cmat_D, cluster = get_cluster(mod_C),
             type = "CR3", test = "chi-sq"),
  Wald_test(mod_C, constraints = Cmat_D, vcov = "CR3", test = "chi-sq")$Fstat
)

tinytest::expect_equal(
  get_boot_F(mod_C, y_boot = y_boot, C_mat = Cmat_E, cluster = get_cluster(mod_C),
             type = "CR2", test = "Naive-Fp"),
  Wald_test(mod_C, constraints = Cmat_E, vcov = "CR2", test = "Naive-Fp")$Fstat
)

tinytest::expect_equal(
  get_boot_F(mod_C, y_boot = y_boot, C_mat = Cmat_F, cluster = get_cluster(mod_C),
             type = "CR0", test = "HTA"),
  Wald_test(mod_C, constraints = Cmat_F, vcov = "CR0", test = "HTA")$Fstat
)


# run_cwb options work for robu objects.

boot_yi_A <- run_cwb(mod_A, cluster = get_cluster(mod_A), R = 9,
                     auxiliary_dist = "Rademacher")
fit_A <- get_fitted(mod_A)
abs_res_A <- abs(get_res(mod_A))
A_check <- sapply(boot_yi_A, function(x) all.equal(abs(x - fit_A), abs_res_A))
tinytest::expect_true(all(A_check))

boot_yi_B <- run_cwb(mod_B, cluster = get_cluster(mod_B), R = 9,
                     auxiliary_dist = "Rademacher")
fit_B <- get_fitted(mod_B)
abs_res_B <- abs(get_res(mod_B))
B_check <- sapply(boot_yi_B, function(x) all.equal(abs(x - fit_B), abs_res_B))
tinytest::expect_true(all(B_check))

boot_yi_C <- run_cwb(mod_C, cluster = get_cluster(mod_C), R = 9,
                     auxiliary_dist = "Rademacher")
fit_C <- get_fitted(mod_C)
abs_res_C <- abs(get_res(mod_C))
C_check <- sapply(boot_yi_C, function(x) all.equal(abs(x - fit_C), abs_res_C))
tinytest::expect_true(all(C_check))


# Wald_test_cwb() results do not depend on sort order.

ord <- sample(1:nrow(oswald2013))
oswald_scram <- oswald2013[ord,]

scram_A <- robu(yi ~ 1,
                studynum = Study, var.eff.size = vi,
                modelweights = "CORR",
                data = oswald_scram)
tinytest::expect_equal(coef(mod_A), coef(scram_A))
tinytest::expect_equal(get_res(mod_A)[ord], get_res(scram_A))
tinytest::expect_equal(get_fitted(mod_A)[ord], get_fitted(scram_A))
tinytest::expect_equal(get_cluster(mod_A)[ord], get_cluster(scram_A))

scram_B <- robu(yi ~ 0 + Crit.Cat,
                studynum = Study, var.eff.size = vi,
                modelweights = "CORR",
                data = oswald_scram)
tinytest::expect_equal(coef(mod_B), coef(scram_B))
tinytest::expect_equal(get_res(mod_B)[ord], get_res(scram_B))
tinytest::expect_equal(get_fitted(mod_B)[ord], get_fitted(scram_B))
tinytest::expect_equal(get_cluster(mod_B)[ord], get_cluster(scram_B))

scram_C <- robu(yi ~ Crit.Cat + Crit.Domain + IAT.Focus + Scoring,
                studynum = Study, var.eff.size = vi,
                modelweights = "CORR",
                data = oswald_scram)
tinytest::expect_equal(coef(mod_C), coef(scram_C))
tinytest::expect_equal(get_res(mod_C)[ord], get_res(scram_C))
tinytest::expect_equal(get_fitted(mod_C)[ord], get_fitted(scram_C))
tinytest::expect_equal(get_cluster(mod_C)[ord], get_cluster(scram_C))

scram_D <- robu(yi ~ Crit.Domain + IAT.Focus + Scoring,
                studynum = Study, var.eff.size = vi,
                modelweights = "CORR",
                data = oswald_scram)
tinytest::expect_equal(coef(mod_D), coef(scram_D))
tinytest::expect_equal(get_res(mod_D)[ord], get_res(scram_D))
tinytest::expect_equal(get_fitted(mod_D)[ord], get_fitted(scram_D))
tinytest::expect_equal(get_cluster(mod_D)[ord], get_cluster(scram_D))

scram_E <- robu(yi ~ Crit.Cat + IAT.Focus + Scoring,
                studynum = Study, var.eff.size = vi,
                modelweights = "CORR",
                data = oswald_scram)
tinytest::expect_equal(coef(mod_E), coef(scram_E))
tinytest::expect_equal(get_res(mod_E)[ord], get_res(scram_E))
tinytest::expect_equal(get_fitted(mod_E)[ord], get_fitted(scram_E))
tinytest::expect_equal(get_cluster(mod_E)[ord], get_cluster(scram_E))

scram_F <- robu(yi ~ Crit.Cat + Crit.Domain + IAT.Focus,
                studynum = Study, var.eff.size = vi,
                modelweights = "CORR",
                data = oswald_scram)
tinytest::expect_equal(coef(mod_F), coef(scram_F))
tinytest::expect_equal(get_res(mod_F)[ord], get_res(scram_F))
tinytest::expect_equal(get_fitted(mod_F)[ord], get_fitted(scram_F))
tinytest::expect_equal(get_cluster(mod_F)[ord], get_cluster(scram_F))

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


tinytest::expect_equal(attr(orig_A, "original"), attr(scram_A, "original"))
tinytest::expect_equal(attr(orig_A, "bootstraps"), attr(scram_A, "bootstraps"))


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

tinytest::expect_equal(attr(orig_D, "original"), attr(scram_D, "original"))
tinytest::expect_equal(attr(orig_D, "bootstraps"), attr(scram_D, "bootstraps"))

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

tinytest::expect_equal(attr(orig_E, "original"), attr(scram_E, "original"))
tinytest::expect_equal(attr(orig_E, "bootstraps"), attr(scram_E, "bootstraps"))

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

tinytest::expect_equal(attr(orig_F, "original"), attr(scram_F, "original"))
tinytest::expect_equal(attr(orig_F, "bootstraps"), attr(scram_F, "bootstraps"))


# Wald_test_cwb() works with user-weighted robu models.

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

tinytest::expect_inherits(test_wt, "Wald_test_wildmeta")

