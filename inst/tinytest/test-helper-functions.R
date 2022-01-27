suppressPackageStartupMessages(library(clubSandwich))

data("AchievementAwardsRCT", package = "clubSandwich")

# constrain_predictors works as expected.

lm_full <- lm(awarded ~ 0 + school_type + treated + year + sex + siblings + immigrant,
              data = AchievementAwardsRCT)
Xmat <- model.matrix(lm_full)

Cmat <- constrain_equal("school_type", reg_ex = TRUE, coef(lm_full))
X_reduced <- wildmeta:::constrain_predictors(Xmat, Cmat)
a <- residuals(lm.fit(X_reduced, AchievementAwardsRCT$awarded))
b <- residuals(update(lm_full, . ~ . + 1 - school_type))
tinytest::expect_equivalent(a, b)

Cmat <- constrain_zero("year", reg_ex = TRUE, coef(lm_full))
X_reduced <- wildmeta:::constrain_predictors(Xmat, Cmat)
a <- residuals(lm.fit(X_reduced, AchievementAwardsRCT$awarded))
b <- residuals(update(lm_full, . ~ . - year))
tinytest::expect_equivalent(a, b)

Cmat <- constrain_zero(c(4,8,9,10), coef(lm_full))
X_reduced <- wildmeta:::constrain_predictors(Xmat, Cmat)
a <- residuals(lm.fit(X_reduced, AchievementAwardsRCT$awarded))
b <- residuals(update(lm_full, . ~ school_type + year))
tinytest::expect_equivalent(a, b)

lm_skinny <- lm(awarded ~ 0 + school_type, data = AchievementAwardsRCT)
Xmat <- model.matrix(lm_skinny)

Cmat <- constrain_equal("school_type", reg_ex = TRUE, coef(lm_skinny))
X_reduced <- wildmeta:::constrain_predictors(Xmat, Cmat)
a <- residuals(lm.fit(X_reduced, AchievementAwardsRCT$awarded))
b <- residuals(update(lm_skinny, . ~ 1))
tinytest::expect_equivalent(a, b)

Cmat <- constrain_equal(1:2, coef(lm_skinny))
X_reduced <- wildmeta:::constrain_predictors(Xmat, Cmat)
a <- residuals(lm.fit(X_reduced, AchievementAwardsRCT$awarded))
b <- residuals(update(lm_skinny, . ~ I(school_type=="Secular")))
tinytest::expect_equivalent(a, b)


# return_wts returns weights with specified distributions.

moments <- function(x) {
  c(mean(x), mean(x^2), mean(x^3), mean(x^4))
}

set.seed(20220112)

n_clusters <- 100000
dists <- c("Rademacher","Mammen","Webb six","uniform","standard normal")
wt_samples <- lapply(dists, wildmeta:::wild_wts, n_clusters = n_clusters)

tinytest::expect_true(all(lengths(wt_samples) == n_clusters))

# Check Rademacher
Rad_vals <- sort(unique(wt_samples[[1]]))
tinytest::expect_identical(Rad_vals, c(-1, 1))
tinytest::expect_true(all(moments(wt_samples[[1]]) - c(0,1,0,1) < .01))

# Check Mammen
Mam_vals <- sort(unique(wt_samples[[2]]))
tinytest::expect_identical(Mam_vals, c(-(sqrt(5) - 1)/2, (sqrt(5) + 1)/2))
tinytest::expect_true(all(moments(wt_samples[[2]]) - c(0,1,1,2) < .01))

# Check Webb
Webb_vals <- sort(unique(wt_samples[[3]]))
tinytest::expect_equal(Webb_vals, c(-sqrt((3:1) / 2),sqrt((1:3)/2)))
tinytest::expect_true(all(moments(wt_samples[[3]]) - c(0,1,0,7/6) < .01))

# Check uniform
tinytest::expect_true(all(wt_samples[[4]] > -sqrt(3) & wt_samples[[4]] < sqrt(3)))
tinytest::expect_true(all(moments(wt_samples[[4]]) - c(0,1,0,9/5) < .01))

# Check normal
tinytest::expect_true(all(moments(wt_samples[[5]]) - c(0,1,0,3) < .01))
