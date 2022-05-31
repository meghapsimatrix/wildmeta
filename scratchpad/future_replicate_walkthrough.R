library(future.apply)
library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)


devtools::load_all()
source("R/helpers.R")

# need to figure out how to set plan(multisession) within run_cwb()
# any special input needed for future apply
# need the check if below code is running in parallel

# a robu model ------------------------------------------------------------

model <- robu(d ~ 0 + study_type + hrs + test,
              studynum = study,
              var.eff.size = V,
              small = FALSE,
              data = SATcoaching)


# stuff before running the bootstraps -------------------------------------

cluster <- get_cluster(model)
if (!is.factor(cluster)) cluster <- as.factor(cluster)

res <- get_res(model)
pred <- get_fitted(model)

num_cluster <- unique(cluster)

auxiliary_dist <- "Rademacher"
f <- get_boot_F
simplify <- TRUE



# future replicate --------------------------------------------------------

plan(multisession)


bootstraps <- future.apply::future_replicate(n = 10, {

  wts <- return_wts(auxiliary_dist = auxiliary_dist, cluster_var = num_cluster)
  eta <- wts[cluster]
  y_boot <- pred + res * eta


}, simplify = simplify & is.null(f))

if (is.null(f)) {
  return(bootstraps)
}



# future sapply to run get boot F -----------------------------------------

C_mat <- constrain_equal(1:3, coefs = model$b.r)


boot_stats <- future_sapply(bootstraps,
                            FUN = get_boot_F,
                            cluster = cluster,
                            full_model = model,
                            C_mat = C_mat)

boot_stats
