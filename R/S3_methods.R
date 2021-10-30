# run_cwb <- function(full_model,
#                     C_mat,
#                     R,
#                     auxiliary_dist = "Rademacher",
#                     adjust = FALSE) {
#
#   UseMethod("run_cwb")
#
# }
#
#
# Wald_test_cwb <- function(full_mod, constraints, ...) {
#
#   UseMethod("Wald_test_cwb")
# }


estimate_null <- function(full_model,
                          C_mat,
                          R) {

  UseMethod("estimate_null")
}
