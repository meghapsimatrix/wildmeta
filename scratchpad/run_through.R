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

res

res$p_val
res$test

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

null_hand <- update(full_model, yi = d ~ hrs + test)
all.equal(predict(null_model), predict(null_hand))

cluster_id <- clubSandwich:::findCluster.rma.mv(full_model)
cluster_id

#set.seed(11092021)  # does set seed even work with the bootstraps? -
                    # do I need to add an argument where i generate the weights?
boots <- run_cwb(null_model,
                 cluster = cluster_id,
                 R = 12,
                 simplify = FALSE)

boots

#save(boots, file = "scratchpad/boots.Rdata")

# works now with <<-
sapply(boots,
       FUN = get_boot_F,
       full_model = full_model,
       C_mat = C_mat)


map_dbl(boots, .f = get_boot_F, full_model = full_model, C_mat = C_mat)

res <- Wald_test_cwb(full_model = full_model,
                     constraint_matrix = C_mat,
                     R = 12)

res

plot_cwb(res)

