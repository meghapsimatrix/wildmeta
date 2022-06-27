data("SATcoaching", package = "clubSandwich")
suppressPackageStartupMessages(library(robumeta))



full_model <- robu(d ~ 0 + study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   small = FALSE,
                   data = SATcoaching)


res <- Wald_test_cwb(full_model = full_model,
                     constraints = constrain_equal(1:3),
                     R = 99)

test_that("plot() returns a ggplot2 object when run on a Wald_test_wildemeta",{

  skip_if_not_installed("ggplot2")

  x <- plot(res)
  y <-
    plot(res, fill = "purple", alpha = 0.5) +
    ggplot2::theme_light()

  expect_s3_class(x, "ggplot")
  expect_s3_class(y, "ggplot")

})


test_that("plot() throws an error if ggplot2 is not installed.", {

  skip_if_not_installed("mockery")

  mockery::stub(plot.Wald_test_wildmeta, 'ggplot2_is_missing', TRUE)
  expect_error(plot(res))
})
