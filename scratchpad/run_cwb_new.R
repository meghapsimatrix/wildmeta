library(metafor)
library(robumeta)

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
run_cwb_new(
  full_model,
  cluster = full_model$data.full$study,
  R = 12
)

# With CR2 adjustment
run_cwb_new(
  full_model,
  cluster = full_model$data.full$study,
  R = 12,
  adjust = "CR2"
)

# Verify wild bootstrap process
bs <- run_cwb_new(full_model, cluster = full_model$data.full$study, R = 20)

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

run_cwb_new(
  full_model,
  cluster = cluster_id,
)

# With CR2 adjustment
run_cwb_new(
  full_model,
  cluster = cluster_id,
  adjust = "CR2"
)

# Verify wild bootstrap process
bs <- run_cwb_new(full_model, cluster = cluster_id)

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
