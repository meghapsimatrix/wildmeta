library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

devtools::load_all()


source("R/S3_methods.R")
source("R/helpers.R")
source("R/run_cwb.R")  # check this function :D
source("R/plot_wildmeta.R")
source("R/Wald_test_wildmeta.R") # check this function too :D



full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = coef(full_model))
constraint_matrix <- C_mat

# See R/helpers.R for constrain_predictors
null_model <- estimate_null(full_model, C_mat)

cluster_id <- clubSandwich:::findCluster.rma.mv(null_model)

bootstraps <- run_cwb(
  null_model,
  cluster =  cluster_id,
  R = 12,
  adjust = "CR2",
  simplify = FALSE
)

boot_stats <- lapply(bootstraps,
                     FUN = get_boot_F.rma.mv,
                     full_model = full_model,
                     C_mat = constraint_matrix)



# robumeta ----------------------------------------------------------------

model <- robu(d ~ 0 + study_type + hrs + test,
              studynum = study,
              var.eff.size = V,
              small = FALSE,
              data = SATcoaching)


bootstraps <- run_cwb(
  model = full_model,
  cluster =  full_model$data.full$study,
  R = 2,
  adjust = "CR2",
  simplify = FALSE
)

bootstraps

boot_stats <- lapply(bootstraps,
                     FUN = get_boot_F.robu,
                     full_model = full_model,
                     C_mat = constraint_matrix)
