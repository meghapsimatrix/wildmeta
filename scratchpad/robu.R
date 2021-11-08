fitted.robu <- function(object, ...) {
  as.numeric(object$data.full$pred)
}

residuals.robu <- function(object, ...) {
  as.numeric(object$data.full$e)
}
