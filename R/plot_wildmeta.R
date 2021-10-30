plot.wildmeta <- function(boots, ...){

  bootstraps <- tibble::tibble(boot_F = boots)

  ggplot2::ggplot(bootstraps, ggplot2::aes(x = boot_F)) +
    ggplot2::geom_density(...) +
    ggplot2::labs(x = "F_statistic", y = "Density") +
    ggplot2::scale_x_continuous() +
    ggplot2::scale_y_continuous() +
    ggplot2::theme_bw()

}
