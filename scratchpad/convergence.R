library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = ~ study_type| study,
                     data = SATcoaching)

full_model

tryCatch(rma.mv(yi = d ~ 0 + study_type + hrs + test,
                V = V,
                random = ~ study_type| study,
                data = SATcoaching), error = function(e) NA)

tryCatch(rma.mv(yi = d ~ what,
                V = V,
                random = ~ study_type | study,
                data = SATcoaching), error = function(e) NA)


possibly_rma <- possibly(rma.mv, otherwise = NA, quiet = FALSE)

possibly_rma(yi = d ~ 0 + study_type + hrs + test,
             V = V,
             random = ~ study_type | study,
             data = SATcoaching)


possibly_lm <- possibly(lm, otherwise = NA)

possibly_lm(d ~ 0 +  study_type + hrs + test, data = SATcoaching)
possibly_lm(d ~ 0 + what, data = SATcoaching)


possibly_robu <- possibly(robu, otherwise = NA)

full_model <- robu(d ~ 0 + study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   small = FALSE,
                   data = SATcoaching)

possibly_robu(robu(d ~ 0 + study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   small = FALSE,
                   data = SATcoaching))



safermamv <-  safely(.f = rma.mv)

full_model <- safermamv(rma.mv(yi = d ~ 0 + study_type + hrs + test,
                               V = V,
                               random = ~ study_type| study,
                               data = SATcoaching))

full_model$result


SATcoaching <- SATcoaching %>%
  group_by(study) %>%
  mutate(esid = row_number()) %>%
  ungroup()


full_model <- rma.mv(yi = d ~ 0 + study_type + hrs + test,
                     V = V,
                     random = list(~ study_type| study, ~ study_type|esid),
                     struct = c("DIAG", "DIAG"),
                     data = SATcoaching)
