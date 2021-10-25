library(clubSandwich)
library(robumeta)
library(tidyverse)

source("R/helpers.R")

set.seed(12102020)


full_model <- robu(d ~ 0 + study_type + hrs + test,
             studynum = study,
             var.eff.size = V,
             small = FALSE,
             data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = full_model$b.r)

dep <- full_model$modelweights

X_mat <- full_model$X.full %>%
  select(-1) %>%
  as.matrix()

effect_size <- full_model$data.full$effect.size
v <- full_model$data.full$var.eff.size
study <- full_model$data.full$study

Xnull <- constrain_predictors(Xmat = X_mat, Cmat = C_mat)

# doesn't like update
update(full_model, formula = effect_size ~ Xnull)

# some issue here
robumeta::robu(effect_size ~ Xnull,
               studynum = study,
               var.eff.size = v,
               small = FALSE,
               modelweights = dep)

# this works
robumeta::robu(effect_size ~ hrs + test,
               studynum = study,
               var.eff.size = v,
               small = FALSE,
               modelweights = dep,
               data = dat)


