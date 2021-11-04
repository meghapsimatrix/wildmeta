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
#Error in solve.default(sumXWX) :
#  system is computationally singular: reciprocal condition number = 3.96487e-17

dat <- tibble(effect_size = effect_size,
              study = study,
              v = v)

dat <- bind_cols(dat, as.data.frame(Xnull))

null_formula <- paste("effect_size ~ 0 + ", paste(colnames(as.data.frame(Xnull)), collapse = " + "))


null_mod <- robumeta::robu(as.formula(null_formula),
               studynum = study,
               var.eff.size = v,
               small = FALSE,
               modelweights = dep,
               data = dat)

# this works
hrs <- full_model$X.full$hrs
test <- full_model$X.full$test

dat <- dat %>%
  mutate(hrs = hrs,
         test = test)

null_mod_hand <- robumeta::robu(effect_size ~ 0 + hrs + test,
               studynum = study,
               var.eff.size = v,
               small = FALSE,
               modelweights = dep,
               data = dat)

all.equal(residuals(null_model), residuals(null_mod_hand))




