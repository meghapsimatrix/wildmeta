library(metafor)

dat <- dat.konstantopoulos2011

mod <- rma.mv(yi = yi ~ year,
              V = vi,
              random = ~ 1 | district / study,
              data = dat)
mod
y_new <- rnorm(nrow(dat))

update(mod, formula = y_new ~ .)


mod <- rma.mv(yi = yi, mods =  ~ year,
              V = vi,
              random = ~ 1 | district / study,
              data = dat)

update(mod, formula = y_new ~ .)
