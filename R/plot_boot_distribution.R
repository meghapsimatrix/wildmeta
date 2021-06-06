#' @title Plot distribution of bootstrap replicates
#'
#' @description Plots distributions of bootstrap replicates, i.e, meta-regression coefficients obtained from different bootstrap samples.
#'
#' @param boot_dat results from the cwb() function.
#'
#'
#' @return A density plot showing the distribution of bootstrap replicates.
#'
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' library(clubSandwich)
#' library(robumeta)
#'
#' full <- robu(d ~ study_type,
#'              studynum = study,
#'              var.eff.size = V,
#'              small = FALSE,
#'              data = SATcoaching)
#'
#' cwb_dat <- cwb(full_model = full,
#'                indices = 2:3,
#'                R = 99)
#'
#' plot_boot_distribution(boot_dat = cwb_dat)
#'

plot_boot_distribution <- function(boot_dat, ...){

  bootstraps <- tibble::tibble(boot_F = unlist(boot_dat$boot_F))

  ggplot2::ggplot(bootstraps, ggplot2::aes(x = boot_F)) +
    ggplot2::geom_density(...) +
    ggplot2::labs(x = "F_statistic", y = "Density") +
    ggplot2::theme_bw()



}
