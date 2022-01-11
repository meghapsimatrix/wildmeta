#' @title Plot distribution of bootstrap test statistics
#'
#' @description Creates a density plot showing the distribution of bootstrap test statistics.
#'
#' @param results Results from Wald_test_cwb function
#' @param ... Any other arguments to be passed to \code{ggplot2::geom_density())}
#'
#'
#' @return A ggplot2 density plot.
#'
#' @export
#'
#' @examples
#' library(clubSandwich)
#' library(robumeta)
#'
#' model <- robu(d ~ 0 + study_type + hrs + test,
#'              studynum = study,
#'               var.eff.size = V,
#'               small = FALSE,
#'               data = SATcoaching)
#'
#' res <- Wald_test_cwb(full_model = full_model,
#'                      constraint_matrix = C_mat,
#'                      R = 12)
#'
#' plot(res, fill = "darkred", alpha = 0.5)
#'


#' @export


plot.Wald_test_wildmeta <- function(results, ...) {

  boots <- attributes(results)$bootstraps
  org_F <- attributes(results)$original

  bootstraps <- data.frame(boot_F = boots)

  ggplot2::ggplot(bootstraps, ggplot2::aes(x = boot_F)) +
    ggplot2::geom_density(...) +
    ggplot2::geom_vline(xintercept = org_F, linetype = "dashed") +
    ggplot2::labs(x = "F_statistic", y = "Density") +
    ggplot2::scale_x_continuous() +
    ggplot2::scale_y_continuous() +
    ggplot2::theme_bw()

}
