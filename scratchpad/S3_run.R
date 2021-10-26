library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

source("R/S3_methods.R")
source("R/helpers.R")
source("scratchpad/rma_mv.R")  # check this function :D
source("scratchpad/plot_wildmeta.R")
source("scratchpad/robu.R") # check this function too :D
source("scratchpad/Wald_test_wildmeta.R")


# metafor -----------------------------------------------------------------

# no intercept

full_model<- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = coef(full_model))

<<<<<<< HEAD
# JAMES sometimes this throws convergence issues -
# so something like safely or something? - and output the # of bootstraps successfully run in waldtest
=======
# sometimes i get convergence issues
>>>>>>> 745f0f3ecebcf9cbb8db20a5d32d3cb0cb0e8640
boots <- run_cwb(full_model,
                 C_mat,
                 R = 99)

plot(boots, fill = "darkred", alpha = 0.6)


# need to figure out how to make Wald_test_cwb talk to run_cwb :D
Wald_test_cwb(full_model,
              C_mat,
              R = 99)

# compare to club CR2
Wald_test(full_model,
          C_mat,
          vcov = "CR2")


#intercept

full_model <- rma.mv(yi = d ~ study_type + hrs + test,
                    V = V,
                    random = ~ study_type| study,
                    data = SATcoaching)

C_mat <- constrain_zero(2:3, coefs = coef(full_model))

# sometimes this throws convergence issues -
# so something like safely or something?
boots <- run_cwb(full_model,
                 C_mat,
                 R = 99)

plot(boots, fill = "darkred", alpha = 0.6)

Wald_test_cwb(full_model,
              C_mat,
              R = 99)


# robumeta ----------------------------------------------------------------

# no intercept

full_model <- robu(d ~ 0 + study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   small = FALSE,
                   data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = full_model$b.r)


boots <- run_cwb(full_model,
                 C_mat,
                 R = 99)

plot(boots, fill = "darkred", alpha = 0.6)

Wald_test_cwb(full_model,
              C_mat,
              R = 99)

Wald_test(full_model,
          C_mat,
          vcov = "CR2")

# intercept
full_model <- robu(d ~ study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   small = FALSE,
                   data = SATcoaching)

C_mat <- constrain_zero(2:3, coefs = full_model$b.r)


boots <- run_cwb(full_model,
                 C_mat,
                 R = 99)

plot(boots, fill = "darkred", alpha = 0.6)


Wald_test_cwb(full_model,
              C_mat,
              R = 99)
