library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

devtools::load_all()

source("R/helpers.R")
source("R/plot_wildmeta.R") #why won't this get loaded?

# robumeta ----------------------------------------------------------------

full_model <- robu(d ~ 0 + study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   small = FALSE,
                   data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = full_model$b.r)

null_model <- estimate_null(full_model,
                            C_mat)

cluster_id <- null_model$data.full$study


boots <- run_cwb(null_model,
                 cluster = cluster_id,
                 R = 12,
                 simplify = FALSE)

sapply(boots,
       FUN = get_boot_F,
       full_model = full_model,
       C_mat = C_mat)

res <- Wald_test_cwb(full_model = full_model,
              constraint_matrix = C_mat,
              R = 12)

# need to do the whole print thing

res$p_val

str(res)

attributes(res)$bootstraps

plot_cwb(res, fill = "darkred", alpha = 0.5)


# metafor -----------------------------------------------------------------


full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = coef(full_model))

null_model <- estimate_null(full_model,
                            C_mat)

cluster_id <- clubSandwich:::findCluster.rma.mv(full_model)

#set.seed(11092021)  # does set seed even work with the bootstraps? -
                    # do I need to add an argument where i generate the weights?
boots <- run_cwb(null_model,
                 cluster = cluster_id,
                 R = 12,
                 simplify = FALSE)


# why is it all the same?
sapply(boots,
       FUN = get_boot_F,
       full_model = full_model,
       C_mat = C_mat)

# is it only returning the last value?


get_boot_F.rma.mv(full_model,
                  y_boot = boots[[2]],
                  C_mat = C_mat)




res <- Wald_test_cwb(full_model = full_model,
                     constraint_matrix = C_mat,
                     R = 12)

res

plot_cwb(res)

# i think it's just returning the last F?


y_boot <- boots[[1]]

y_new <- rep(NA, length = nrow(full_model$X.f))
y_new[full_model$not.na] <- y_boot


boot_mod <- tryCatch(update(full_model, formula = y_new ~ .),
                     error = function(e) NA)


cov_mat <- clubSandwich::vcovCR(boot_mod, type = "CR1")
res <- clubSandwich::Wald_test(boot_mod,
                               constraints = C_mat,
                               vcov = cov_mat,
                               test = "Naive-F")

res <- res$Fstat
res




y_boot <- boots[[12]]

y_new <- rep(NA, length = nrow(full_model$X.f))
y_new[full_model$not.na] <- y_boot


boot_mod <- tryCatch(update(full_model, formula = y_new ~ .),
                     error = function(e) NA)


cov_mat <- clubSandwich::vcovCR(boot_mod, type = "CR1")
res <- clubSandwich::Wald_test(boot_mod,
                               constraints = C_mat,
                               vcov = cov_mat,
                               test = "Naive-F")

res <- res$Fstat
res




