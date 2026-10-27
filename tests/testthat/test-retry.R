skip_if_not_installed("future")
skip_if_not_installed("parallelly")
skip_if_not_installed("future.apply")
library(future)

suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(clubSandwich))


# failed_bootstrap <- readRDS("tests/testthat/testdata/cwb_failure_debug.rds")
failed_bootstrap <- readRDS("testdata/cwb_failure_debug.rds")
dat <- failed_bootstrap$dat

V_list <- vcalc(vi = var_g_j, cluster = study_id, rho = 0.7, obs = esid, sparse = TRUE, data = dat)

rma_fit <- rma.mv(
  yi =  g, V = V_list,
  mods = ~ moderator,
  random = ~ 1 | study_id / esid,
  data = dat,
  sparse = TRUE,
  control = list(optimizer="nlminb")
)

C <- length(unique(dat$moderator))



test_that("Wald_test_cwb() works with `retry > 1L` when future.apply package is available.", {

  if (parallelly::supportsMulticore()) {
    plan(multicore)
  } else {
    plan(multisession)
  }

  wald_CWB_results <- Wald_test_cwb(
    full_model = rma_fit,
    constraints = constrain_zero(2:C),
    R = 69,
    seed = 459191506
  )
  full_booties <- attr(wald_CWB_results, "bootstraps")
  na_vals <- is.na(full_booties)
  expect_gt(sum(na_vals), 0L)

  Wald_CWB_retries <- Wald_test_cwb(
    full_model = rma_fit,
    constraints = constrain_zero(2:C),
    R = 69,
    retry = 4L,
    seed = 459191506
  )

  retry_booties <- attr(Wald_CWB_retries, "bootstraps")
  na_retries <- sum(is.na(retry_booties))

  expect_identical(na_retries, 0L)

  expect_identical(full_booties[!na_vals], retry_booties[!na_vals])

})

test_that("Wald_test_cwb() works with `retry > 1L` when future.apply package is absent.", {

  wald_CWB_results <- testthat::with_mocked_bindings(
    Wald_test_cwb(
      full_model = rma_fit,
      constraints = constrain_zero(2:C),
      R = 69,
      seed = 459191508
    ),
    future_apply_is_available = \(x) FALSE
  )

  full_booties <- attr(wald_CWB_results, "bootstraps")
  na_vals <- is.na(full_booties)
  expect_gt(sum(na_vals), 0L)

  Wald_CWB_retries <- testthat::with_mocked_bindings(
    Wald_test_cwb(
      full_model = rma_fit,
      constraints = constrain_zero(2:C),
      R = 69,
      retry = 4L,
      seed = 459191506
    ),
    future_apply_is_available = \(x) FALSE
  )

  retry_booties <- attr(Wald_CWB_retries, "bootstraps")
  na_retries <- sum(is.na(retry_booties))

  expect_identical(na_retries, 0L)

  expect_identical(full_booties[!na_vals], retry_booties[1:sum(!na_vals)])
})
