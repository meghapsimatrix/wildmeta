library(future.apply)
library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)
library(wildmeta)



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

n_clusters <- length(unique(cluster))


# future replicate --------------------------------------------------------

plan(multisession)


system.time(bootstraps <- future.apply::future_replicate(n = 10, {

  wts <- wild_wts(auxiliary_dist = auxiliary_dist, n_clusters = n_clusters)
  eta <- wts[cluster]
  y_boot <- pred + res * eta


}, simplify = simplify & is.null(f)))



system.time(bootstraps <- replicate(n = 10, {

  wts <- wild_wts(auxiliary_dist = auxiliary_dist, n_clusters = n_clusters)
  eta <- wts[cluster]
  y_boot <- pred + res * eta


}, simplify = simplify & is.null(f)))

if (is.null(f)) {
  return(bootstraps)
}



# future sapply to run get boot F -----------------------------------------

C_mat <- constrain_equal(1:3, coefs = model$b.r)


system.time(boot_stats <- future_sapply(bootstraps,
                            FUN = get_boot_F,
                            cluster = cluster,
                            full_model = model,
                            C_mat = C_mat))

system.time(boot_stats <- sapply(bootstraps,
                                 FUN = get_boot_F,
                                 cluster = cluster,
                                 full_model = model,
                                 C_mat = C_mat))

boot_stats




# check function -----------------------------------------------------------

robu_res <- Wald_test_cwb(full_model = model,
                          constraints = C_mat,
                          R = 10)

robu_res


# metafor -----------------------------------------------------------------



rma_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                    V = V,
                    random = ~ study_type | study,
                    data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = coef(rma_model))
model <- rma_model
full_model <- rma_model

rma_res <- Wald_test_cwb(full_model = rma_model,
                         constraints = constrain_equal(1:3),
                        R = 100)

rma_res

constraints <- C_mat

if (inherits(constraints, "function")) {
  constraints <- constraints(stats::coef(full_model))
}

# compute the null model
null_model <- estimate_null(full_model,
                            C_mat = constraints)

if (is.null(cluster)) cluster <- get_cluster(null_model)

seed <- 20220606
adjust <- "CR1"

metafor_boots <-
  run_cwb(null_model,
           cluster = cluster,
           R = R,
           f = get_boot_F,  # this goes to sapply
           full_model = full_model,
           C_mat = constraints,
           type = type,
           test = test,
           auxiliary_dist = auxiliary_dist,
           adjust = adjust,
           simplify = TRUE,
           seed = seed)



#if (!is.null(seed)) set.seed(seed)

# coerce cluster variable to factor
if (!is.factor(cluster)) cluster <- as.factor(cluster)

# # residuals and predicted values ------------------------------------------

res <- get_res(model)
pred <- get_fitted(model)

# Adjust ------------------------------------------------------------------



# bootstrap ---------------------------------------------------------------
n_clusters <- length(unique(cluster))

R <- 10
bootstraps <- future.apply::future_replicate(n = R, {

  wts <- wild_wts(auxiliary_dist = auxiliary_dist, n_clusters = n_clusters)
  eta <- wts[cluster]
  y_boot <- pred + res * eta


}, simplify = simplify & is.null(f))

if (is.null(f)) {
  return(bootstraps)
}

bootstraps

# use future sapply
system.time(boot_stats <- sapply(bootstraps,
                                 FUN = get_boot_F,
                                 cluster = cluster,
                                 full_model = model,
                                 C_mat = C_mat))


type <- "CR0"
test <- "Naive-F"

full_vcov <- clubSandwich::vcovCR(full_model, type = type, cluster = cluster)
org_F <- clubSandwich::Wald_test(full_model,
                                 constraints = constraints,
                                 vcov = full_vcov,
                                 test = test)

org_F <- org_F$Fstat

p_val <- mean(boots > org_F, na.rm = TRUE)
boot_test <- if (adjust != "CR0") "CWB Adjusted" else "CWB"

p_boot <- data.frame(
  Test = boot_test,
  Adjustment = adjust,
  CR_type = type,
  Statistic = test,
  R = R,
  p_val = p_val
)

class(p_boot) <- c("Wald_test_wildmeta", class(p_boot))
attr(p_boot, "bootstraps") <- boots
attr(p_boot, "original") <- org_F

