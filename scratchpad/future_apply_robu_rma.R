library(future.apply)
library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

devtools::load_all()

# a robu model ------------------------------------------------------------

robu_model <- robu(d ~ 0 + study_type + hrs + test,
              studynum = study,
              var.eff.size = V,
              small = FALSE,
              data = SATcoaching)

system.time(robu_res <- Wald_test_cwb(full_model = robu_model,
                          constraints = constrain_equal(1:3),
                          R = 99))

robu_res


# metafor -----------------------------------------------------------------

rma_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                    V = V,
                    random = ~ study_type | study,
                    data = SATcoaching)

system.time(rma_res <- Wald_test_cwb(full_model = rma_model,
                                     constraints = constrain_equal(1:3),
                                     R = 99))

rma_res

# breaks when doing parallel? ---------------------------------------------
plan(multisession)

system.time(rma_res <- Wald_test_cwb(full_model = rma_model,
                         constraints = constrain_equal(1:3),
                         R = 99))


rma_res
