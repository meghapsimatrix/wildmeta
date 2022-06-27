#' @title Plot distribution of bootstrap test statistics
#'
#' @description Creates a density plot showing the distribution of bootstrap test statistics.
#'
#' @param x Results from Wald_test_cwb function
#' @param ... Any other arguments to be passed to \code{ggplot2::geom_density()}
#'
#'
#' @return A ggplot2 density plot.
#'
#' @export
#'
#' @examples
#' data("SATcoaching", package = "clubSandwich")
#' library(clubSandwich)
#' library(robumeta)
#'
#' full_model <- robu(d ~ 0 + study_type + hrs + test,
#'                    studynum = study,
#'                    var.eff.size = V,
#'                    small = FALSE,
#'                    data = SATcoaching)
#'
#'
#' res <- Wald_test_cwb(full_model = full_model,
#'                      constraints = constrain_equal(1:3),
#'                      R = 99)
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#' plot(res, fill = "darkred", alpha = 0.5)
#' }
#'


plot.Wald_test_wildmeta <- function(x, ...) {

  if (ggplot2_is_missing()) stop("The plot() function requires the ggplot2 package. Please install it.", call. = FALSE)

  boots <- attributes(x)$bootstraps
  org_F <- attributes(x)$original

  bootstraps <- data.frame(boot_F = boots)
  stat <- paste(x$Statistic, "statistic")

  ggplot2::ggplot(bootstraps, ggplot2::aes_string(x = "boot_F")) +
    ggplot2::geom_density(...) +
    ggplot2::geom_vline(xintercept = org_F, linetype = "dashed") +
    ggplot2::labs(x = stat, y = "Density") +
    ggplot2::scale_x_continuous() +
    ggplot2::scale_y_continuous() +
    ggplot2::theme_minimal()

}

ggplot2_is_missing <- function() {
  !requireNamespace("ggplot2", quietly = TRUE)
}
