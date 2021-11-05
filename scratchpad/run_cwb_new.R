library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

<<<<<<< HEAD
devtools::load_all()

#--------------------------------------------------------------
# Using robu model

full_model <- robu(d ~ 0 + study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   small = FALSE,
                   data = SATcoaching)

# Not working yet because of some sort of namespace thing
fitted.values(full_model)
residuals(full_model)

full_model$fitted.values <- fitted.robu(full_model)
full_model$residuals <- residuals.robu(full_model)
fitted.values(full_model)
residuals(full_model)

# No adjustment
run_cwb(
  full_model,
  cluster = full_model$data.full$study,
  R = 12
)

# With CR2 adjustment
run_cwb(
  full_model,
  cluster = full_model$data.full$study,
  R = 12,
  adjust = "CR2"
)


C_mat <- constrain_equal(1:3, coefs = full_model$b.r)
constraint_matrix <- C_mat
R <- 12
adjust <- "CR0"
auxiliary_dist <- "Rademacher"
simplify <- FALSE
f <- get_boot_F

# added the null model
null_model <- estimate_null.robu(full_model,
                            C_mat = constraint_matrix,
                            R = R)
=======
# Drop rows with missing predictors
dat <- SATcoaching[complete.cases(SATcoaching),]
>>>>>>> c4b3dc51601af04f1b04e057c2a84640ffb37123

full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = dat)
summary(full_model)

X <- full_model$X
Cmat <- constrain_equal(1:3, coefs = coef(full_model))
C_mat <- constrain_equal(1:3, coefs = coef(full_model))

# See R/helpers.R for constrain_predictors
Xnull <- constrain_predictors(X, Cmat)

null_auto <- update(full_model, yi = full_model$yi, mods = ~ 0 + Xnull)
null_auto

bootstraps <- run_cwb(
  null_model,
  cluster = full_model$data.full$study,
  R = 12,
  adjust = "CR2"
)

boot_stats <- lapply(bootstraps,
                     FUN = get_boot_F.robu,
                     full_model = full_model,
                     C_mat = constraint_matrix)
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

<<<<<<< HEAD
full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                      V = V,
                      random = ~ study_type| study,
                      data = SATcoaching)

fitted.values(full_model)
residuals(full_model)

cluster_id <- clubSandwich:::findCluster.rma.mv(full_model)

run_cwb(
  full_model,
  cluster = cluster_id,
  R = 12
)

# With CR2 adjustment
run_cwb(
  full_model,
  cluster = cluster_id,
  adjust = "CR2",
  R = 12
)

# Verify wild bootstrap process
bootstraps <- run_cwb(full_model, cluster = cluster_id, R = 12)

bootstraps %>%
  # back out the auxiliary random variables
  map_dfc(~ (.x - fitted.values(full_model)) / residuals(full_model)) %>%
  # check that auxiliaries are constant within cluster
  mutate(cluster = cluster_id) %>%
  group_by(cluster) %>%
  summarize(
    across(everything(), ~ diff(range(.x))),
    .groups = "drop"
  ) %>%
  select(-cluster) %>%
  unlist() %>%
  sd()

constraint_matrix <- constrain_equal(1:3, coefs = as.numeric(full_model$b))

get_boot_F.rma.mv(y_boot = bootstraps[[1]], full_model = full_model, C_mat = constraint_matrix)

boot_stats <- lapply(bootstraps,
                     FUN = get_boot_F.rma.mv,
                     full_model = full_model,
                     C_mat = constraint_matrix)
=======
#all.equal(predict(null_hand)$pred, predict(null_model)$pred)


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
>>>>>>> c4b3dc51601af04f1b04e057c2a84640ffb37123
