library(clubSandwich)
library(metafor)
library(robumeta)

devtools::load_all()

source("R/helpers.R")

get_boot_stats <- function(y_boot,
                           model,
                           C_mat){


  y_new <- rep(NA, length = nrow(model$X.f))
  y_new[model$not.na] <- y_boot


  boot_mod <- update(model, formula = y_new ~ .)


  cov_mat <- clubSandwich::vcovCR(boot_mod, type = "CR1")

  wald_res <- clubSandwich::Wald_test(boot_mod,
                                      constraints = C_mat,
                                      vcov = cov_mat,
                                      test = "Naive-F")

  #resolution <- wald_res$Fstat



  return(wald_res)



}



full <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ 1| study,
                     data = SATcoaching)

constraints <- constrain_equal(1:3, coefs = coef(full))

# null_model <- estimate_null(full_model,
#                             C_mat)
#
# cluster_id <- clubSandwich:::findCluster.rma.mv(full_model)

# boots <- run_cwb(null_model,
#                  cluster = cluster_id,
#                  R = 12,
#                  simplify = FALSE)
#
# save(boots, file = "scratchpad/boots.Rdata")

load("scratchpad/boots.Rdata")

rm(res)
#y_boot <- boots[[12]]

#why is it spitting out the same number

get_boot_stats(y_boot = boots[[2]], model = full, C_mat = constraints)
get_boot_stats(y_boot = boots[[12]], model = full, C_mat = constraints)



