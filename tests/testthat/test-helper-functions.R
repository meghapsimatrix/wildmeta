suppressPackageStartupMessages(library(clubSandwich))

data("AchievementAwardsRCT", package = "clubSandwich")

test_that("constrain_predictors works as expected.", {

  lm_full <- lm(awarded ~ 0 + school_type + treated + year + sex + siblings + immigrant,
                data = AchievementAwardsRCT)
  Xmat <- model.matrix(lm_full)

  Cmat <- constrain_equal("school_type", reg_ex = TRUE, coef(lm_full))
  X_reduced <- constrain_predictors(Xmat, Cmat)
  a <- residuals(lm.fit(X_reduced, AchievementAwardsRCT$awarded))
  b <- residuals(update(lm_full, . ~ . + 1 - school_type))
  expect_equal(a, b, ignore_attr = TRUE)

  C_wrong <- cbind(diag(nrow = 3), diag(rep(0, 3)))
  expect_error(constrain_predictors(Xmat, C_wrong))

  Cmat <- constrain_zero("year", reg_ex = TRUE, coef(lm_full))
  X_reduced <- constrain_predictors(Xmat, Cmat)
  a <- residuals(lm.fit(X_reduced, AchievementAwardsRCT$awarded))
  b <- residuals(update(lm_full, . ~ . - year))
  expect_equal(a, b, ignore_attr = TRUE)

  Cmat <- constrain_zero(c(4,8,9,10), coef(lm_full))
  X_reduced <- constrain_predictors(Xmat, Cmat)
  a <- residuals(lm.fit(X_reduced, AchievementAwardsRCT$awarded))
  b <- residuals(update(lm_full, . ~ school_type + year))
  expect_equal(a, b, ignore_attr = TRUE)

  lm_skinny <- lm(awarded ~ 0 + school_type, data = AchievementAwardsRCT)
  Xmat <- model.matrix(lm_skinny)

  Cmat <- constrain_equal("school_type", reg_ex = TRUE, coef(lm_skinny))
  X_reduced <- constrain_predictors(Xmat, Cmat)
  a <- residuals(lm.fit(X_reduced, AchievementAwardsRCT$awarded))
  b <- residuals(update(lm_skinny, . ~ 1))
  expect_equal(a, b, ignore_attr = TRUE)

  Cmat <- constrain_equal(1:2, coef(lm_skinny))
  X_reduced <- constrain_predictors(Xmat, Cmat)
  a <- residuals(lm.fit(X_reduced, AchievementAwardsRCT$awarded))
  b <- residuals(update(lm_skinny, . ~ I(school_type=="Secular")))
  expect_equal(a, b, ignore_attr = TRUE)

})

moments <- function(x) {
  c(mean(x), mean(x^2), mean(x^3), mean(x^4))
}

test_that("return_wts returns weights with specified distributions.", {

  set.seed(20220112)

  n_clusters <- 100000
  dists <- c("Rademacher","Mammen","Webb six","uniform","standard normal")
  wt_samples <- lapply(dists, wild_wts, n_clusters = n_clusters)

  expect_true(all(lengths(wt_samples) == n_clusters))

  # Check Rademacher
  Rad_vals <- sort(unique(wt_samples[[1]]))
  expect_identical(Rad_vals, c(-1, 1))
  expect_true(all(moments(wt_samples[[1]]) - c(0,1,0,1) < .01))

  # Check Mammen
  Mam_vals <- sort(unique(wt_samples[[2]]))
  expect_identical(Mam_vals, c(-(sqrt(5) - 1)/2, (sqrt(5) + 1)/2))
  expect_true(all(moments(wt_samples[[2]]) - c(0,1,1,2) < .01))

  # Check Webb
  Webb_vals <- sort(unique(wt_samples[[3]]))
  expect_equal(Webb_vals, c(-sqrt((3:1) / 2),sqrt((1:3)/2)))
  expect_true(all(moments(wt_samples[[3]]) - c(0,1,0,7/6) < .01))

  # Check uniform
  expect_true(all(wt_samples[[4]] > -sqrt(3) & wt_samples[[4]] < sqrt(3)))
  expect_true(all(moments(wt_samples[[4]]) - c(0,1,0,9/5) < .01))

  # Check normal
  expect_true(all(moments(wt_samples[[5]]) - c(0,1,0,3) < .01))

})
