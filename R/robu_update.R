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

CE_handmade <- function(X, y, v, cluster, rho) {

  k_j <- as.numeric(table(cluster))
  m <- length(k_j)

  sigma_sq_j <- tapply(v, cluster, mean)

  w_tilde <- 1 / (k_j * sigma_sq_j)
  w_tilde_j <- as.numeric(w_tilde[cluster])

  mod_prelim <- lm.wfit(x = X, y = y, w = w_tilde_j)

  res <- residuals(mod_prelim)

  # calculate weighted residual sum of squares
  QE <- sum(w_tilde_j * res^2)

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
  mod_CE <- lm.wfit(x = X, y = y, w = w_ij)

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

  class(res) <- c("handmade.robu.CE","handmade.robu")

  return(res)

}

HE_handmade <- function(X, y, v, cluster) {

  w_ij <- 1 / v

  mod_prelim <- lm.wfit(x = X, y = y, w = w_ij)

  res <- residuals(mod_prelim)

  # calculate sums of squares
  QE <- sum(w_ij * res^2)
  Q1 <- sum(tapply(res, cluster, sum)^2)

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
  mod_HE <- lm.wfit(x = X, y = y, w = w_ij)

  res <- mod_HE[c("coefficients","residuals","fitted.values","weights")]
  res$tau_sq <- tau_sq
  res$omega_sq <- omega_sq

  res$X <- X
  res$v <- v
  res$w_ij <- w_ij
  res$cluster <- cluster
  res$QE <- QE
  res$Q1 <- Q1
  res$constants <- list(A1 = A1, B1 = B1, C1 = C1, A2 = A2, B2 = B2, C2 = C2, p = p)
  res$more <- list(XJX = XJX, XWJWX = XWJWX, XJWX = XJWX,
                   kjXJWX = kjXJWX, VXJX = VXJX,
                   XWWX = XWWX, trW = trW)
  res$nobs <- N

  class(res) <- c("handmade.robu.HE","handmade.robu")

  return(res)

}

#-------------------------------------------------------------------------------
# Methods for handmade.robu objects

model.matrix.handmade.robu <- function(object, ...) object$X

bread.handmade.robu <- function(x, ...) {
  x$nobs * chol2inv(chol(crossprod(x$X, x$w_ij * x$X)))
}

#-------------------------------------------------------------------------------
# Updating methods for robu objects
# Re-fitted model is a handmade.robu object

update_robu <- function(mod, y) {
  UseMethod("update_robu")
}

#' @export

update_robu.default <- function(mod, y) {

  modelweights <- switch(mod$mod_label[[1]],
                         "RVE: Correlated Effects Model" = "CORR",
                         "RVE: Hierarchical Effects Model" = "HIER",
                         "RVE: User Specified Weights" = "user",)

  cluster <- as.factor(mod$study_orig_id)
  resort <- order(order(cluster))
  X <- mod$Xreg[resort,,drop=FALSE]
  v <- mod$data.full$var.eff.size[resort]

  # need to fix that y might be longer than X due to missing values

  if (modelweights == "CORR") {
    res <- CE_handmade(X = X, y = y, v = v,
                       cluster = cluster,
                       rho = mod$mod_info$rho)
  } else if (modelweights == "HIER") {
    res <- HE_handmade(X = X, y = y, v = v,
                       cluster = cluster)
  } else if (modelweights == "user") {

  }

  res
}

#-------------------------------------------------------------------------------
# Updating methods for handmade robu objects

#' @export

update_robu.handmade.robu.CE <- function(mod, y) {

  mod_prelim <- lm.wfit(x = mod$X, y = y, w = mod$w_tilde_j)
  res <- residuals(mod_prelim)
  QE <- sum(mod$w_tilde_j * res^2)
  tau_sq <- (QE - mod$num_minus) / mod$den
  tau_sq <- ifelse(tau_sq < 0, 0, tau_sq)  # added this

  # CE weights
  w_j <- 1 / (mod$k_j * (mod$sigma_sq_j + tau_sq))
  w_ij <- rep(w_j, mod$k_j)

  # fit WLS regression
  mod_CE <- lm.wfit(x = mod$X, y = y, w = w_ij)

  res <- mod_CE[c("coefficients","residuals","fitted.values","weights")]
  res$tau_sq <- tau_sq

  if (calc_CR0) {
    res$vcov <- calc_CR0(X = mod$X, w_ij = w_ij,
                         resid = res$residuals,
                         cluster = mod$cluster)
  }

  return(res)
}

#' @export

update_robu.handmade.robu.HE <- function(mod, y) {

}

#' @export

update_robu.handmade.robu.user <- function(mod, y) {

}
