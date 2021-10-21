library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)


full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = SATcoaching)
summary(full_model)

X <- full_model$X
Cmat <- constrain_equal(1:3, coefs = coef(full_model))

q <- nrow(Cmat)
p <- ncol(Cmat)
XtX_inv <- chol2inv(chol(crossprod(X)))
Cnull <- diag(nrow = p) - XtX_inv %*% t(Cmat) %*% chol2inv(chol(Cmat %*% XtX_inv %*% t(Cmat))) %*% Cmat
qr(X %*% Cnull)
Xnull <- qr.X(qr(X %*% Cnull), ncol = p - q)

null_model <- update(full_model, yi = full_model$yi, mods = ~ 0 + Xnull)
null_model


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
