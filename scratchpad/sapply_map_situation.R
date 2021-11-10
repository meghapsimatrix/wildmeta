library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

devtools::load_all()

source("R/helpers.R")
source("R/plot_wildmeta.R")


load("scratchpad/boots.Rdata")

boots

full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = coef(full_model))

sapply(boots,
       FUN = get_boot_F,
       full_model = full_model,
       C_mat = C_mat)


sapply(boots,
       FUN = get_boot_F.rma.mv,
       full_model = full_model,
       C_mat = C_mat)

y_boot <- boots[[3]]

y_new <- rep(NA, length = nrow(full_model$X.f))
y_new[full_model$not.na] <- y_boot


boot_mod <- tryCatch(update(full_model, formula = y_new ~ .),
                     error = function(e) NA)

if(inherits(boot_mod, "rma.mv")){

  cov_mat <- clubSandwich::vcovCR(boot_mod, type = "CR1")

  res <- clubSandwich::Wald_test(boot_mod,
                                 constraints = C_mat,
                                 vcov = cov_mat,
                                 test = "Naive-F")

  res <- res$Fstat


} else{

  res <- NA
}

res





