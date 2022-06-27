suppressPackageStartupMessages(library(clubSandwich))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(robumeta))

SATcoaching_full <- subset(SATcoaching, !is.na(hrs), !is.na(test))

robu_mod <- robu(d ~ hrs + test,
                  studynum = study,
                  var.eff.size = V,
                  small = FALSE,
                  data = SATcoaching_full)

rma_mod <- rma.mv(d ~ hrs + test,
                  V = V,
                  random = ~ 1 | study,
                  data = SATcoaching_full)

rma_sub <- rma.mv(d ~ hrs + test,
                  V = V,
                  random = ~ 1 | study,
                  data = SATcoaching,
                  subset = !is.na(hrs) & !is.na(test))

suppressWarnings(
  rma_miss <- rma.mv(d ~ hrs + test,
                     V = V,
                     random = ~ 1 | study,
                     data = SATcoaching)
)

test_that("run_cwb works without setting a future plan.",{

  robu_res <- run_cwb(robu_mod,
                      cluster = get_cluster(robu_mod),
                      R = 10, seed = 10, simplify = TRUE)

  rma_res <- run_cwb(rma_mod,
                     cluster = get_cluster(rma_mod),
                     R = 8, seed = 12, simplify = TRUE)

  expect_identical(ncol(robu_res), 10L)
  expect_identical(ncol(rma_res), 8L)

  robu_test <- Wald_test_cwb(robu_mod,
                             constraints = constrain_zero(2:3),
                             R = 10)

  expect_s3_class(robu_test, "Wald_test_wildmeta")
  expect_true(!is.na(robu_test$p_val))

  rma_test <- Wald_test_cwb(rma_mod,
                            constraints = constrain_zero(2:3),
                            R = 21, seed = 5)

  expect_s3_class(rma_test, "Wald_test_wildmeta")
  expect_true(!is.na(rma_test$p_val))

  sub_test <- Wald_test_cwb(rma_sub,
                            constraints = constrain_zero(2:3),
                            R = 21, seed = 5)

  expect_s3_class(sub_test, "Wald_test_wildmeta")
  expect_true(!is.na(sub_test$p_val))

  suppressWarnings(
    miss_test <- Wald_test_cwb(rma_miss,
                               constraints = constrain_zero(2:3),
                               R = 21, seed = 5)
  )

  expect_s3_class(miss_test, "Wald_test_wildmeta")
  expect_true(!is.na(miss_test$p_val))

  expect_equal(rma_test, sub_test)
  expect_equal(rma_test, miss_test)

})

test_that("run_cwb returns the same results with plan(sequential) and plan(multisession).", {

  skip_on_cran()

  skip_if_not_installed("future")
  skip_if_not_installed("parallelly")
  skip_if_not_installed("future.apply")

  library(future)

  f <- function(x, cluster, time = 0.01) {
    Sys.sleep(time)
    x
  }

  plan(sequential)

  time_robu_seq <- system.time(
    robu_seq <- run_cwb(robu_mod,
                        cluster = get_cluster(robu_mod),
                        R = 60,
                        seed = 10, simplify = TRUE),
  )

  time_rma_seq <- system.time(
    rma_seq <- run_cwb(rma_mod,
                       cluster = get_cluster(rma_mod),
                       R = 60,
                       seed = 12, simplify = TRUE)
  )

  if (parallelly::supportsMulticore()) {
    plan(multicore)
  } else {
    plan(multisession)
  }

  time_robu_multi <- system.time(
    robu_multi <- run_cwb(robu_mod,
                        cluster = get_cluster(robu_mod),
                        R = 60,
                        f = f,
                        time = 0.5,
                        seed = 10, simplify = TRUE)
  )

  time_rma_multi <- system.time(
    rma_multi <- run_cwb(rma_mod,
                       cluster = get_cluster(rma_mod),
                       R = 60,
                       f = f,
                       time = 0.5,
                       seed = 12, simplify = TRUE)
  )

  expect_equal(robu_seq, robu_multi)
  expect_equal(rma_seq, rma_multi)

  expect_lt(time_robu_multi[3], 30)
  expect_lt(time_rma_multi[3], 30)

})

test_that("Wald_test_cwb() returns the same results with plan(sequential) and plan(multisession).", {

  skip_on_cran()

  skip_if_not_installed("future")
  skip_if_not_installed("parallelly")
  skip_if_not_installed("future.apply")

  library(future)

  plan(sequential)

  robu_seq <- Wald_test_cwb(robu_mod,
                            constraints = constrain_zero(2:3),
                            R = 12, seed = 100)

  rma_seq <- Wald_test_cwb(rma_mod,
                           constraints = constrain_zero(2:3),
                           R = 18, seed = 101)

  sub_seq <- Wald_test_cwb(rma_sub,
                           constraints = constrain_zero(2:3),
                           R = 18, seed = 101)

  suppressWarnings(
    mis_seq <- Wald_test_cwb(rma_miss,
                             constraints = constrain_zero(2:3),
                             R = 18, seed = 101)
  )

  expect_equal(rma_seq, sub_seq)
  expect_equal(rma_seq, mis_seq)

  if (parallelly::supportsMulticore()) {
    plan(multicore)
  } else {
    plan(multisession)
  }

  robu_multi <- Wald_test_cwb(robu_mod,
                              constraints = constrain_zero(2:3),
                              R = 12, seed = 100)

  rma_multi <- Wald_test_cwb(rma_mod,
                             constraints = constrain_zero(2:3),
                             R = 18, seed = 101)

  sub_multi <- Wald_test_cwb(rma_sub,
                             constraints = constrain_zero(2:3),
                             R = 18, seed = 101)

  suppressWarnings(
    mis_multi <- Wald_test_cwb(rma_miss,
                               constraints = constrain_zero(2:3),
                               R = 18, seed = 101)
  )

  expect_equal(robu_seq, robu_multi)
  expect_equal(rma_seq, rma_multi)
  expect_equal(sub_seq, sub_multi)
  expect_equal(mis_seq, mis_multi)

})

test_that("run_cwb uses future_args if specified.", {

  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

})
