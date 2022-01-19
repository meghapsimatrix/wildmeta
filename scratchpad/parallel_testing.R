library(future.apply)
library(clubSandwich)
library(robumeta)

source("R/helpers.R")
source("R/robu.R")
source("R/S3_methods.R")
source("R/robu_update.R")
source("R/plot_wildmeta.R")

set.seed(12102020)

full_model <- robu(d ~ 0 + study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   small = FALSE,
                   data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = full_model$b.r)


run_cwb <- function(model,
                    cluster,
                    R,
                    f = NULL,
                    ...,
                    auxiliary_dist = "Rademacher",
                    adjust = "CR0",
                    simplify = FALSE) {


  # coerce cluster variable to factor
  if (!is.factor(cluster)) cluster <- as.factor(cluster)

  # # residuals and predicted values ------------------------------------------

  res <- get_res(model)
  pred <- get_fitted(model)

  # Adjust ------------------------------------------------------------------

  if (adjust %in% c("CR1", "CR2", "CR3", "CR4")) {
    split_res <- split(res, cluster)
    B_j <- attr(clubSandwich::vcovCR(model,
                                     cluster = cluster,
                                     type = adjust), "adjustments")
    res_list <- Map(function(x, y) as.vector(x %*% y), x = B_j, y = split_res)
    res <- unsplit(res_list, cluster)
  }

  # bootstrap ---------------------------------------------------------------
  num_cluster <- unique(cluster)

  plan(multisession, workers = 4)

  bootstraps <- future_replicate(n = 12, {

    wts <- return_wts(auxiliary_dist = auxiliary_dist)
    eta <- wts[cluster]
    y_boot <- pred + res * eta


  }, simplify = simplify & is.null(f))

  if (is.null(f)) {
    return(bootstraps)
  }

  boot_stats <- sapply(bootstraps, f, cluster = cluster, ..., simplify = simplify)

  return(boot_stats)
}




Wald_test_cwb <- function(full_model,
                          constraints,
                          R,
                          cluster = NULL,
                          auxiliary_dist = "Rademacher",
                          adjust = "CR0",
                          type = "CR0",
                          test = "Naive-F") {

  if (inherits(constraints, "function")) {
    constraints <- constraints(stats::coef(full_model))
  }

  # added the null model
  null_model <- estimate_null(full_model,
                              C_mat = constraints)

  # for run_cwb_new need to pull out the clusters
  # will there be an issue with missing data in clusters for rma.mv?
  if (is.null(cluster)) cluster <- get_cluster(null_model)


  boots <- run_cwb(null_model,
                   cluster = cluster,
                   R = R,
                   f = get_boot_F,  # this goes to sapply
                   full_model = full_model,
                   C_mat = constraints,
                   type = type,
                   test = test,
                   auxiliary_dist = auxiliary_dist,
                   adjust = adjust,
                   simplify = TRUE)

  full_vcov <- clubSandwich::vcovCR(full_model, type = type, cluster = cluster)
  org_F <- clubSandwich::Wald_test(full_model,
                                   constraints = constraints,
                                   vcov = full_vcov,
                                   test = test)

  org_F <- org_F$Fstat

  p_val <- mean(boots > org_F, na.rm = TRUE)
  boot_test <- if (adjust != "CR0") "CWB Adjusted" else "CWB"

  p_boot <- data.frame(
    Test = boot_test,
    Adjustment = adjust,
    `CR type` = type,
    `Statistic` = test,
    R = R,
    p_val = p_val
  )

  class(p_boot) <- c("Wald_test_wildmeta", class(p_boot))
  attr(p_boot, "bootstraps") <- boots
  attr(p_boot, "original") <- org_F

  return(p_boot)

}


Wald_test_cwb(full_model = full_model,
              constraints = C_mat,
              R = 10000)
