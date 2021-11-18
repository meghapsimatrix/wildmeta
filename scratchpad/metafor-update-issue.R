library(metafor)

update_with_noise <- function(mod, noise_sd) {
  y_noise <- rnorm(n = nobs(mod), sd = noise_sd)
  mod_e <- new.env(parent = parent.frame(mod$call$yi))
  mod_e$y_noise <- y_noise
  environment(mod$call$yi) <- mod_e
  update(mod, formula = y_noise ~ .)
}

dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
res_uni <- rma(yi ~ year, vi, data=dat, method="REML")

update_with_noise(res_uni, noise_sd = 0.5)


dat.long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
levels(dat.long$group) <- c("exp", "con")
dat.long$group <- relevel(dat.long$group, ref="con")
dat.long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat.long)

res_mv <- rma.mv(yi, vi,
                  mods = ~ year + group,
                  random = ~ group | study, struct="UN",
                  data=dat.long)

update_with_noise(res_mv, noise_sd = 0.5)
