Wald_test_cwb <- function(full_model,
                      constraints,
                      R,
                      auxiliary_dist = "Rademacher",
                      adjust = FALSE){

  boots <- run_cwb(full_model,
                   C_mat,
                   R = R,
                   auxiliary_dist = auxiliary_dist,
                   adjust = adjust)


  org_F <- clubSandwich::Wald_test(full_model,
                                   constraints = C_mat,
                                   vcov = clubSandwich::vcovCR(full_model, type = "CR1"),
                                   test = "Naive-F")

  org_F <- org_F$Fstat


  # JEP: You don't need the dplyr for any of this!
  p_val <- mean(Fstat > org_F)
  test <- "CWB"


  if (adjust == TRUE) {
    test <- "CWB Adjusted"
  }

  p_boot <- list(p_val, test)


  # output boots too somehow - list
  # and figure out print?
  attr(p_boot, "bootstraps") <- boots
  attr(p_boot, "test") <- test

  class(p_boot) <- "Wald_test_wildmeta"

  return(p_boot)

}
