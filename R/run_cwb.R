run_cwb <- function(model,
                    cluster,
                    f = NULL,
                    ...,
                    R,
                    auxiliary_dist = "Rademacher",
                    adjust = "CR0",
                    simplify = FALSE) {


  # coerce cluster variable to factor
  if (!is.factor(cluster)) cluster <- as.factor(cluster)

  # # residuals and predicted values ------------------------------------------

  # need to figure out the robu situation
  # model$fitted.values <- fitted.robu(model)
  # model$residuals <- residuals.robu(model)

  res <- get_res(model)
  pred <- get_fitted(model)

  # Adjust ------------------------------------------------------------------

  if (adjust %in% c("CR1", "CR2", "CR3", "CR4")) {
    split_res <- split(res, cluster)
    B_j <- attr(clubSandwich::vcovCR(model,
                                     cluster = cluster,
                                     type = adjust), "adjustments")
    res_list <- purrr::map2(B_j, split_res, ~ as.vector(.x %*% .y))
    res <- unsplit(res_list, cluster)
  }

  # bootstrap ---------------------------------------------------------------
  num_cluster <- unique(cluster)

  bootstraps <- replicate(n = R, {

    wts <- return_wts(auxiliary_dist = auxiliary_dist, cluster_var = num_cluster)
    eta <- wts[cluster]
    y_boot <- pred + res * eta

    #y_new <- rep(NA, length = nrow(model$X.f))  # how is this working with robu?
    #y_new[model$not.na] <- y_boot

  }, simplify = simplify & is.null(f))

  if (is.null(f)) {
    return(bootstraps)
  }

  boot_stats <- sapply(bootstraps, f, ..., simplify = simplify) # lapply doesn;t have simplify as argument

  return(boot_stats)
}
