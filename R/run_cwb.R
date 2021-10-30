run_cwb <- function(model,
                        cluster,
                        f = NULL,
                        ...,
                        R = 10,
                        auxiliary_dist = "Rademacher",
                        adjust = "CR0",
                        simplify = FALSE) {

  # coerce cluster variable to factor
  if (!is.factor(cluster)) cluster <- as.factor(cluster)

  # residuals and predicted values ------------------------------------------

  res <- stats::residuals(model)
  pred <- stats::fitted.values(model)

  # Adjust ------------------------------------------------------------------

  if (adjust %in% c("CR1","CR2","CR3","CR4")) {
    split_res <- split(res, cluster)
    B_j <- attr(clubSandwich::vcovCR(model,
                                     cluster = cluster,
                                     type = adjust), "adjustments")
    res_list <- purrr::map2(B_j, split_res, ~ as.vector(.x %*% .y))
    res <- unsplit(res_list, cluster)
  }

  # bootstrap ---------------------------------------------------------------
  num_cluster <- unique(cluster)
  k_j <- as.numeric(table(cluster))

  bootstraps <- replicate(n = R, {

    wts <- return_wts(auxiliary_dist = auxiliary_dist, cluster_var = num_cluster)
    eta <- wts[cluster]
    y_boot <- pred + res * eta

    y_new <- rep(NA, length = nrow(model$X.f))
    y_new[model$not.na] <- y_boot

  }, simplify = simplify & is.null(f))

  if (is.null(f)) {
    return(bootstraps)
  }

  boot_stats <- lapply(bootstraps, f, ..., simplify = simplify)

  return(boot_stats)
}
