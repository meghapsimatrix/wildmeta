library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = SATcoaching)


posslm <-  possibly(.f = lm, otherwise = "Error")

posslm(lm(d ~ 1, data = SATcoaching))


possrmamv <-  possibly(.f = rma.mv, otherwise = "Error")

full_model <- possrmamv(rma.mv(yi = d ~ 0 + study_type + hrs + test,
                               V = V,
                               random = ~ study_type| study,
                               data = SATcoaching))




safermamv <-  safely(.f = rma.mv)

full_model <- safermamv(rma.mv(yi = d ~ 0 + study_type + hrs + test,
                               V = V,
                               random = ~ study_type| study,
                               data = SATcoaching))

full_model$error
