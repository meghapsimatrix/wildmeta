library(future.apply)
library(wildmeta)
library(clubSandwich)
library(robumeta)

source("R/helpers.R")
source("R/robu.R")
source("R/S3_methods.R")
source("R/robu_update.R")
source("R/plot_wildmeta.R")




y <- future_replicate(100, mean(rexp(10)))



seed <- 12102020



full_model <- robu(d ~ 0 + study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   small = FALSE,
                   data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = full_model$b.r)

cluster <- get_cluster(full_model)
model <- full_model
adjust <- "CR0"
auxiliary_dist <- "Rademacher"
simplify <- FALSE
f <- NULL


if (!is.null(seed)) set.seed(seed)

# coerce cluster variable to factor
if (!is.factor(cluster)) cluster <- as.factor(cluster)

# # residuals and predicted values ------------------------------------------

res <- get_res(model)
pred <- get_fitted(model)

# Adjust ------------------------------------------------------------------

if (adjust %in% c("CR1", "CR2", "CR3", "CR4")) {
  split_res <- split(res, cluster)
  B_j <- attr(clubSandwich::vcovCR(model,
                                   cluster = cluster,
                                   type = adjust), "adjustments")
  res_list <- Map(function(x, y) as.vector(x %*% y), x = B_j, y = split_res)
  res <- unsplit(res_list, cluster)
}

# bootstrap ---------------------------------------------------------------
n_clusters <- length(unique(cluster))


# some kind of like if parallel then do following if not then use replicate or something?

plan(multisession, workers = 4)

bootstraps <- future_replicate(n = 12, {

  wts <- wild_wts(auxiliary_dist = auxiliary_dist, n_clusters = n_clusters)
  eta <- wts[cluster]
  y_boot <- pred + res * eta


}, simplify = simplify & is.null(f))

if (is.null(f)) {
  return(bootstraps)
}

boot_stats <- sapply(bootstraps, f, cluster = cluster, ..., simplify = simplify)
