library(metafor)
library(robumeta)
library(clubSandwich)

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


# added the null model
null_model <- estimate_null(full_model,
                            C_mat = constraints,
                            R = R)

null_model

null_model$fitted.values <- fitted.robu(null_model)
null_model$residuals <- residuals.robu(null_model)

cluster <- get_cluster(null_model)

run_cwb(
  null_model,
  cluster = cluster,
  R = 12
)

# Verify wild bootstrap process
bs <- run_cwb(full_model, cluster = full_model$data.full$study, R = 20)

library(tidyverse)
bs %>%
  # back out the auxiliary random variables
  map_dfc(~ (.x - fitted.values(full_model)) / residuals(full_model)) %>%
  # check that auxiliaries are constant within cluster
  mutate(cluster = full_model$data.full$study) %>%
  group_by(cluster) %>%
  summarize(
    across(everything(), ~ diff(range(.x))),
    .groups = "drop"
  ) %>%
  select(-cluster) %>%
  unlist() %>%
  sd()

#--------------------------------------------------------------
# Using rma.mv model

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
bs <- run_cwb(full_model, cluster = cluster_id, R = 12)

bs %>%
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
