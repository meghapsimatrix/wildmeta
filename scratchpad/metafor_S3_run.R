library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

source("R/S3_methods.R")
source("R/helpers.R")
source("scratchpad/rma_mv.R")  # check this function :D
source("scratchpad/plot_wildmeta.R")

full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = coef(full_model))


boots <- run_CWB(full_model,
                 C_mat,
                 R = 99)

plot(boots, fill = "darkred", alpha = 0.6)
