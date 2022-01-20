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
  x <- plot(res)
  y <-
    plot(res, fill = "purple", alpha = 0.5) +
    ggplot2::theme_light()

  expect_s3_class(x, "ggplot")
  expect_s3_class(y, "ggplot")

})
