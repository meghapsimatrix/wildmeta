do_CWB <- function(null_model, clusters, boots, auxiliary_dist, adjusted) {

  # Set up the bootstrap process

  # Generate a list of bootstrapped outcome vectors

}


# this is going to be what the users use right?
Wald_test_CWB <- function(full_mod, constraints, …) {

  # Compute the null mod based on full_mod, constraints

  # Compute the test statistic based on full_mod

  yi_boots <- do_CWB(null_mod, …)

  # Compute the test statistic based on bootstrapped outcomes

  # Compute p-values
  # what do you want the output to look like

}
