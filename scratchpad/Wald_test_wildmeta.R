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
                                   test = "Naive-F") %>%
    dplyr::pull(Fstat)


  p_boot <- data.frame(Fstat = as.numeric(boots)) %>%
    dplyr::summarize(p_val = mean(Fstat > org_F)) %>%
    dplyr::mutate(test = "CWB") %>%
    dplyr::select(test, p_val)

  if(adjust == TRUE){
    p_boot$test <- "CWB Adjusted"
  }

  p_boot <- p_boot %>%
    dplyr::select(test,  p_val)

  # output boots too somehow - list
  # and figure out print?

  #class(p_boot) <- "Wald_test_wildmeta"

  return(p_boot)

}
