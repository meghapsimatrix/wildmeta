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
                          C_mat) {

  UseMethod("estimate_null")
}


get_cluster <- function(full_model){

  UseMethod("get_cluster")

}


get_boot_F <- function(full_model,
                       y_boot,
                       C_mat) {

  UseMethod("get_boot_F")
}


get_fitted <- function(model){

  UseMethod("get_fitted")

}

get_res <- function(model){

  UseMethod("get_res")

}

