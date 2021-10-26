run_cwb <- function(full_model,
                    C_mat,
                    R,
                    auxiliary_dist = "Rademacher",
                    adjust = FALSE) {

  UseMethod("run_CWB")

}


Wald_test_cwb <- function(full_mod, constraints, ...) {

  UseMethod(Wald_test_CWB)
}
