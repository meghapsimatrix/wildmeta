#------------------------------------------------------
# Matrix functions
#------------------------------------------------------

# split a matrix
mat_split <- function(mat, f) {
  lapply(levels(f), function(l) mat[f==l,,drop=FALSE])
}

# trace product
trace_product <- function(A, B) {

  a_vec <- as.vector(t(A))
  b_vec <- as.vector(B)
  sum(a_vec * b_vec)

}

#------------------------------------------------------
# robu model-fitting functions
#------------------------------------------------------
#' @importFrom stats lm.wfit
#' @importFrom stats residuals

CE_handmade <- function(X, y, v, cluster, rho, vcov = NULL) {

  k_j <- as.numeric(table(cluster))
  m <- length(k_j)

  sigma_sq_j <- tapply(v, cluster, mean)

  w_tilde <- 1 / (k_j * sigma_sq_j)
  w_tilde_j <- as.numeric(w_tilde[cluster])

  mod_prelim <- stats::lm.wfit(x = X, y = y, w = w_tilde_j)

  resid <- stats::residuals(mod_prelim)

  # calculate weighted residual sum of squares
  QE <- sum(w_tilde_j * resid^2)

  # Create M tilde ----------------------------------------------------------
  X_j <- mat_split(X, cluster)
  M_tilde <- chol2inv(chol(crossprod(X, w_tilde_j * X)))
  p_j <- lapply(X_j, colSums)

  # trace products ----------------------------------------------------------

  # the first B for the trace product numerator
  w_over_k <- w_tilde / k_j
  B_num_all_1 <- Map(function(w, X) w * crossprod(X), w = w_over_k, X = X_j)
  B_num_1 <- Reduce(`+`, B_num_all_1)

  # the second B for the trace product numerator
  XJX_j <- lapply(p_j, tcrossprod)
  B_num_all_2 <- Map(function(w, xjx) w * xjx, w = w_over_k, xjx = XJX_j)
  B_num_2 <- Reduce(`+`, B_num_all_2)

  # the B for the trade product denominator
  B_den_all <- Map(function(w, xjx) w^2 * xjx, w = w_tilde, xjx = XJX_j)
  B_den <- Reduce(`+`, B_den_all)

  num_minus <- m - (1 - rho) * trace_product(M_tilde, B_num_1) - rho * trace_product(M_tilde, B_num_2)
  den <- sum(k_j * w_tilde) - trace_product(M_tilde, B_den)

  tau_sq <- (QE - num_minus) / den
  tau_sq <- ifelse(tau_sq < 0, 0, tau_sq)

  # CE weights
  w_j <- 1 / (k_j * (sigma_sq_j + tau_sq))
  w_ij <- as.numeric(w_j[cluster])

  # fit WLS regression
  mod_CE <- stats::lm.wfit(x = X, y = y, w = w_ij)

  res <- mod_CE[c("coefficients","residuals","fitted.values","weights")]
  res$tau_sq <- tau_sq

  res$X <- X
  res$w_ij <- w_ij
  res$cluster <- cluster
  res$k_j <- k_j
  res$sigma_sq_j <- sigma_sq_j
  res$w_tilde_j <- w_tilde_j
  res$num_minus <- num_minus
  res$den <- den
  res$nobs <- sum(k_j)
  res$vcov_type <- vcov

  class(res) <- c("handmade.robu.CE","handmade.robu")

  if (!is.null(vcov)) {
    res$vcov <- clubSandwich::vcovCR(res, cluster = cluster, type = vcov)
  }

  return(res)

}

HE_handmade <- function(X, y, v, cluster, vcov = NULL) {

  w_ij <- 1 / v

  mod_prelim <- stats::lm.wfit(x = X, y = y, w = w_ij)

  resid <- stats::residuals(mod_prelim)

  # calculate sums of squares
  QE <- sum(w_ij * resid^2)
  Q1 <- sum(tapply(resid, cluster, sum)^2)

  # other constants used in variance component calculations
  X_j <- mat_split(X, cluster)
  k_j <- as.numeric(table(cluster))
  w_j <- split(w_ij, cluster)
  WX <- w_ij * X
  WX_j <- mat_split(WX, cluster)
  XWX <- crossprod(X, WX)
  M_tilde <- chol2inv(chol(XWX))
  JX_j <- do.call(rbind, sapply(X_j, colSums, simplify = FALSE))
  JWX_j <- do.call(rbind, sapply(WX_j, colSums, simplify = FALSE))
  XJX <- crossprod(JX_j)
  XWJWX <- crossprod(JWX_j)
  XWWX <- crossprod(WX)
  trW <- sum(w_ij)
  N <- nrow(X)
  p <- ncol(X)

  XJWX <- crossprod(JX_j, JWX_j)
  kjXJWX <- crossprod(k_j * JX_j, JWX_j)
  VXJX <- M_tilde %*% XJX

  A1 <- sum(k_j^2) -
    trace_product(M_tilde, kjXJWX + t(kjXJWX)) +
    trace_product(VXJX, M_tilde %*% XWJWX)
  B1 <- N -
    trace_product(M_tilde, XJWX + t(XJWX)) +
    trace_product(VXJX, M_tilde %*% XWWX)
  C1 <- sum(v) - trace_product(M_tilde, XJX)
  A2 <- trW - trace_product(M_tilde, XWJWX)
  B2 <- trW - trace_product(M_tilde, XWWX)
  C2 <- N - p

  # variance component estimates
  omega_sq <- (A2 * (Q1 - C1) - A1 * (QE - C2)) / (B1 * A2 - B2 * A1)
  omega_sq <- max(0, omega_sq)
  tau_sq <- (QE - C2) / A2 - omega_sq * B2 / A2
  tau_sq <- max(0, tau_sq)

  # HE weights
  w_ij <- 1 / (v + omega_sq + tau_sq)

  # fit WLS regression
  mod_HE <- stats::lm.wfit(x = X, y = y, w = w_ij)

  res <- mod_HE[c("coefficients","residuals","fitted.values","weights")]
  res$tau_sq <- tau_sq
  res$omega_sq <- omega_sq

  res$X <- X
  res$v <- v
  res$w_ij <- w_ij
  res$cluster <- cluster
  res$constants <- list(A1 = A1, B1 = B1, C1 = C1, A2 = A2, B2 = B2, C2 = C2, p = p)
  res$nobs <- N
  res$vcov_type <- vcov

  class(res) <- c("handmade.robu.HE","handmade.robu")

  if (!is.null(vcov)) {
    res$vcov <- clubSandwich::vcovCR(res, cluster = cluster, type = vcov)
  }

  return(res)

}

user_handmade <- function(X, y, v, weights, cluster, vcov = NULL) {

  # fit WLS regression
  mod_user <- stats::lm.wfit(x = X, y = y, w = weights)

  res <- mod_user[c("coefficients","residuals","fitted.values","weights")]

  res$X <- X
  res$v <- v
  res$w_ij <- weights
  res$cluster <- cluster
  res$nobs <- nrow(X)
  res$vcov_type <- vcov

  class(res) <- c("handmade.robu.user","handmade.robu")

  if (!is.null(vcov)) {
    res$vcov <- clubSandwich::vcovCR(res, cluster = cluster, type = vcov)
  }

  return(res)

}


#-------------------------------------------------------------------------------
# Methods for handmade.robu and handmade_robu objects

#' @importFrom stats model.matrix
#' @export

model.matrix.handmade.robu <- function(object, ...) {
  object$X
}

#' @export

model.matrix.handmade_robu <- function(object, ...) {
  object$X
}

#' @importFrom sandwich bread
#' @export
#'

bread.handmade.robu <- function(x, ...) {
  x$nobs * chol2inv(chol(crossprod(x$X, x$w_ij * x$X)))
}

#' @export

bread.handmade_robu <- function(x, ...) {
  x$nobs * chol2inv(chol(crossprod(x$X, x$w_ij * x$X)))
}

#' @importFrom clubSandwich vcovCR
#' @export

vcovCR.handmade.robu <- function(obj, cluster = obj$cluster,
                                 type = obj$vcov_type,
                                 target = NULL,
                                 inverse_var = NULL,
                                 form = "sandwich", ...) {

  if (is.null(inverse_var)) {
    inverse_var <- is.null(target) & (!inherits(obj, "handmade.robu.user"))
  }

  if (is.null(target) & inherits(obj, "handmade.robu.user")) {
      V <- as.numeric(tapply(obj$v, obj$cluster, mean)[obj$cluster])
      target <- tapply(V, cluster, function(x) diag(x, nrow = length(x)))
  }

  class(obj) <- "handmade_robu"
  clubSandwich::vcovCR(
    obj, cluster = cluster, type = type,
    target = target, inverse_var = inverse_var,
    form = form
  )
}

#-------------------------------------------------------------------------------
# Updating methods for robu objects
# Re-fitted model is a handmade.robu object

update_robu <- function(mod, y, vcov = NULL) {
  UseMethod("update_robu")
}

#' @export

update_robu.default <- function(mod, y, vcov = NULL) {

  lab <- mod$mod_label[[1]]
  if (!is.null(lab)) {
    modelweights <- switch(lab,
                           "RVE: Correlated Effects Model" = "CORR",
                           "RVE: Hierarchical Effects Model" = "HIER",
                           "RVE: User Specified Weights" = "user",
                           "missing")
  }

  if (is.null(lab) || modelweights == "missing") stop("mod must be a robu object.")

  cluster <- as.factor(mod$study_orig_id)
  resort <- order(order(cluster))
  X <- mod$Xreg[resort,,drop=FALSE]
  v <- mod$data.full$var.eff.size[resort]

  # need to fix that y might be longer than X due to missing values

  if (modelweights == "CORR") {
    res <- CE_handmade(X = X, y = y, v = v,
                       cluster = cluster,
                       rho = mod$mod_info$rho, vcov = vcov)
  } else if (modelweights == "HIER") {
    res <- HE_handmade(X = X, y = y, v = v,
                       cluster = cluster, vcov = vcov)
  } else if (modelweights == "user") {
    weights <- mod$data.full$userweights[resort]
    res <- user_handmade(X = X, y = y, v = v,
                         weights = weights,
                         cluster = cluster, vcov = vcov)
  }

  res
}

#-------------------------------------------------------------------------------
# Updating methods for handmade robu objects

#' @export

update_robu.handmade.robu.CE <- function(mod, y, vcov = mod$vcov_type) {

  mod_prelim <- stats::lm.wfit(x = mod$X, y = y, w = mod$w_tilde_j)
  resid <- stats::residuals(mod_prelim)
  QE <- sum(mod$w_tilde_j * resid^2)
  tau_sq <- (QE - mod$num_minus) / mod$den
  tau_sq <- ifelse(tau_sq < 0, 0, tau_sq)  # added this

  # CE weights
  w_j <- 1 / (mod$k_j * (mod$sigma_sq_j + tau_sq))
  w_ij <- as.numeric(w_j[mod$cluster])

  # fit WLS regression
  mod_CE <- stats::lm.wfit(x = mod$X, y = y, w = w_ij)

  res <- mod
  res[c("coefficients","residuals","fitted.values","weights")] <- mod_CE[c("coefficients","residuals","fitted.values","weights")]
  res$tau_sq <- tau_sq
  res$w_ij <- w_ij

  if (!is.null(vcov)) {
    res$vcov <- clubSandwich::vcovCR(res, cluster = mod$cluster, type = vcov)
  }

  return(res)
}

#' @export

update_robu.handmade.robu.HE <- function(mod, y, vcov = mod$vcov_type) {

  w_ij <- 1 / mod$v
  mod_prelim <- stats::lm.wfit(x = mod$X, y = y, w = w_ij)

  resid <- stats::residuals(mod_prelim)

  # calculate sums of squares
  QE <- sum(w_ij * resid^2)
  Q1 <- sum(tapply(resid, mod$cluster, sum)^2)

  # other constants used in variance component calculations
  A1 <- mod$constants$A1
  B1 <- mod$constants$B1
  C1 <- mod$constants$C1
  A2 <- mod$constants$A2
  B2 <- mod$constants$B2
  C2 <- mod$constants$C2

  # variance component estimates
  omega_sq <- (A2 * (Q1 - C1) - A1 * (QE - C2)) / (B1 * A2 - B2 * A1)
  omega_sq <- max(0, omega_sq)
  tau_sq <- (QE - C2) / A2 - omega_sq * B2 / A2
  tau_sq <- max(0, tau_sq)

  # HE weights
  w_ij <- 1 / (mod$v + omega_sq + tau_sq)

  # fit WLS regression
  mod_HE <- stats::lm.wfit(x = mod$X, y = y, w = w_ij)

  res <- mod
  res[c("coefficients","residuals","fitted.values","weights")] <- mod_HE[c("coefficients","residuals","fitted.values","weights")]
  res$tau_sq <- tau_sq
  res$omega_sq <- omega_sq
  res$w_ij <- w_ij

  if (!is.null(vcov)) {
    res$vcov <- clubSandwich::vcovCR(res, cluster = mod$cluster, type = vcov)
  }

  return(res)

}

#' @export

update_robu.handmade.robu.user <- function(mod, y, vcov = mod$vcov_type) {
  user_handmade(X = mod$X, y = y, v = mod$v,
                       weights = mod$weights,
                       cluster = mod$cluster, vcov = vcov)
}
