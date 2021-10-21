run_CWB <- function(full_model,
                    C_mat,
                    R,
                    auxiliary_dist = "Rademacher",
                    adjust = FALSE) {

  UseMethod(do_CWB)

}


Wald_test_CWB <- function(full_mod, constraints, …) {

  UseMethod(Wald_test_CWB)
}
