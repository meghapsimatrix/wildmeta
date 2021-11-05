Wald_test_cwb <- function(full_model,
                          constraint_matrix,
                          R,
                          auxiliary_dist = "Rademacher",
                          adjust = "CR0"){

  # added the null model
  null_model <- estimate_null(full_model,
                              C_mat = constraint_matrix,
                              R = R)

  # for run_cwb_new need to pull out the clusters
  cluster <- get_cluster(null_model)

  boots <- run_cwb(null_model,
                   cluster = cluster,
                   R = R)

  # then need to calculate the bootstrapped F here


  org_F <- clubSandwich::Wald_test(full_model,
                                   constraints = C_mat,
                                   vcov = clubSandwich::vcovCR(full_model, type = "CR1"),
                                   test = "Naive-F")

  org_F <- org_F$Fstat

  p_val <- mean(Fstat > org_F)
  test <- "CWB"


  if (adjust == TRUE) {
    test <- "CWB Adjusted"
  }

  p_boot <- data.frame(test = test, p_val = p_val)

  attr(p_boot, "bootstraps") <- boots
  attr(p_boot, "test") <- test

  class(p_boot) <- "Wald_test_wildmeta"

  return(p_boot)

}
