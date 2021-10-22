library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

# Drop rows with missing predictors
dat <- SATcoaching[complete.cases(SATcoaching),]

full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = dat)
summary(full_model)

X <- full_model$X
Cmat <- constrain_equal(1:3, coefs = coef(full_model))

# See R/helpers.R for constrain_predictors
Xnull <- constrain_predictors(X, Cmat)

null_auto <- update(full_model, yi = full_model$yi, mods = ~ 0 + Xnull)
null_auto

# Check that null_auto gives the same fitted values as directly specifying the null model
null_hand <- update(full_model, yi = d ~ hrs + test)
all.equal(predict(null_auto), predict(null_hand))


# Things get a little bit trickier if there are rows with missing predictors
full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = SATcoaching)
summary(full_model)


X <- full_model$X
Cmat <- constrain_equal(1:3, coefs = coef(full_model))
Xnull <- constrain_predictors(X, Cmat)

update(full_model, yi = full_model$yi, mods = ~ 0 + Xnull)
# Doesn't work because the random effects in full_model involve all
# nrow(SATcoaching) observations, but nrow(X) only includes the
# observations with non-missing values

# This will work for rma.mv() models, I think.
Xnull_f <- matrix(NA, nrow = nrow(full_model$X.f), ncol = ncol(Xnull))
Xnull_f[full_model$not.na,] <- Xnull
null_with_NAs <- update(full_model, yi = full_model$yi.f, mods = ~ 0 + Xnull_f)
# Check that predictions agree with above
all.equal(predict(null_auto)$pred, predict(null_with_NAs)$pred)


indices <- 2

y_dat <- dplyr::bind_cols(full_model$yi,
                          full_model$vi,
                          .name_repair = ~ vctrs::vec_as_names(..., repair = "unique",
                                                               quiet = TRUE)) %>%
  dplyr::rename(effect_size = 1,
                v = 2)

x_dat <- tibble::as_tibble(full_model$X) %>%
  dplyr::select(-1)

study <- tibble::tibble(study = clubSandwich:::findCluster.rma.mv(full_model))

dat <- dplyr::bind_cols(y_dat, x_dat, study)
dat$study <- as.character(dat$study)



null_formula <- paste(rownames(full_model$beta)[ - indices], collapse = " + ")
null_formula <- stringr::str_replace(null_formula, "intrcpt", "1")



full_formula <- paste(rownames(full_model$beta), collapse = " + ")
full_formula <- stringr::str_replace(full_formula, "intrcpt", "1")


null_model <- metafor::rma.mv(yi = stats::as.formula(paste("effect_size ~ ", null_formula)),
                              V = v,
                              random = ~ 1 | study,
                              data = dat)


# this doesn't work - and I think we need the null model to have the same random effects structure and what not
null_model <- update(full_model,
                     yi = stats::as.formula(paste("effect_size ~ ", null_formula)))  # but what if there are other things in the formula


null_model <- update(full_model,
                     yi = stats::as.formula(paste("effect_size ~ ", null_formula)))  # but what if there are other things in the formula
null_model


# it works if i do this but...
effect_size <- dat$effect_size
null_model <- update(full_model,
                     yi = stats::as.formula(paste("effect_size ~ ", null_formula)))  # but what if there are other things in the formula

null_model



# this works fine
yi_boot <- rnorm(nrow(dat))
boot_model <- update(full_model,
                     yi = yi_boot)
