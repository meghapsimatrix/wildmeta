library(metafor)

update_with_noise <- function(mod, noise_sd, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  new_dat <- eval(mod$call$data)

  if (is.null(mod$formula.yi)) {
    yi_name <- as.character(mod$call$yi)
  } else {
    yi_name <- as.character(mod$formula.yi[[2]])
  }

  new_dat[[yi_name]] <- rnorm(n = length(mod$yi), sd = noise_sd)
  update(mod, data = new_dat)
}

# univariate examples

dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
mod <- rma(yi ~ year, vi, data=dat, method="REML")
A <- update_with_noise(mod, noise_sd = 0.5, seed = 4)

mod <- rma(yi = yi, mods = ~ year, vi, data=dat, method="REML")
B <- update_with_noise(res_uni, noise_sd = 0.5, seed = 4)

floating_vector <- dat$yi
mod <- rma(yi = floating_vector, mods = ~ year, vi, data=dat, method="REML")
C <- update_with_noise(res_uni, noise_sd = 0.5, seed = 4)

mod <- rma(yi = floating_vector ~ year, vi, data=dat, method="REML")
D <- update_with_noise(res_uni, noise_sd = 0.5, seed = 4)

identical(A$beta, B$beta)
identical(A$beta, C$beta)
identical(A$beta, D$beta)

# multivariate examples
dat.long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
levels(dat.long$group) <- c("exp", "con")
dat.long$group <- relevel(dat.long$group, ref="con")
dat.long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat.long)

mod <- rma.mv(yi ~ year + group, vi,
                      random = ~ group | study, struct="UN",
                      data=dat.long)

A <- update_with_noise(res_mv_frml, noise_sd = 0.5, seed = 5)

mod <- rma.mv(yi, vi,
                 mods = ~ year + group,
                 random = ~ group | study, struct="UN",
                 data=dat.long)

B <- update_with_noise(res_mv_mods, noise_sd = 0.5, seed = 5)

identical(A$b, B$b)
