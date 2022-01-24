
estimate_null <- function(full_model,
                          C_mat) {

  UseMethod("estimate_null")
}


get_cluster <- function(full_model) {

  UseMethod("get_cluster")

}


#' @importFrom clubSandwich vcovCR
#' @importFrom clubSandwich Wald_test

get_boot_F <- function(full_model,
                       y_boot,
                       C_mat,
                       cluster,
                       type,
                       test) {

  UseMethod("get_boot_F")
}


get_fitted <- function(model) {

  UseMethod("get_fitted")

}

get_res <- function(model) {

  UseMethod("get_res")

}
