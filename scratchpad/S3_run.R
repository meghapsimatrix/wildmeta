library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

# JEP: CTRL + Shift + L to load all package functions
# (equivalent to sourcing entire R/ directory)

source("R/S3_methods.R")
source("R/helpers.R")
source("R/cwb_rma_mv.R")  # check this function :D
source("R/plot_wildmeta.R")
source("R/cwb_robu.R") # check this function too :D
source("R/Wald_test_wildmeta.R")


# metafor -----------------------------------------------------------------

# no intercept

full_model<- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = coef(full_model))

# JAMES sometimes this throws convergence issues -
# so something like safely or something? - and output the # of bootstraps successfully run in waldtest
boots <- run_cwb(full_model,
                 C_mat,
                 R = 99)

# JEP: Show the observed test statistic also
plot(boots, fill = "darkred", alpha = 0.6)


# need to figure out how to make Wald_test_cwb talk to run_cwb :D
# JEP: One way to do this would be to make the results of Wald_test_cwb()
# include the bootstrap distribution, the test statistic, the p-value,
# with a special class of CWB_Wald or something.
# And then write a print.CWB_Wald() method to show just the results of the test,
# a plot.CWB_Wald() method to make the graph, etc.

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
                 R = 12)



plot(boots, fill = "darkred", alpha = 0.6)

check <- Wald_test_cwb(full_model,
                       C_mat,
                       R = 12)

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
