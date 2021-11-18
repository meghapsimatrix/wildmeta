fitted.robu <- function(object, ...) {
  as.numeric(object$data.full$pred.r)
}

residuals.robu <- function(object, ...) {
  as.numeric(object$data.full$e.r)
}
