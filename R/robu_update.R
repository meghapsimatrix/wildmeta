#------------------------------------------------------
# Model-fitting/estimation/testing functions
#------------------------------------------------------

# split a matrix
mat_split <- function(mat, f) {
  lapply(unique(f), function(l) mat[f==l,,drop=FALSE])
}

# trace product 
trace_product <- function(A, B) {
  
  a_vec <- as.vector(t(A))
  b_vec <- as.vector(B)
  sum(a_vec * b_vec)
  
}


robu_handmade <- function(X, y, v, cluster, rho = .8, calc_vcov = NULL) {
  
  k_j <- as.numeric(table(cluster))
  sigma_sq_j <- tapply(v, cluster, mean)
  
  m <- length(k_j)
  
  w_tilde <- 1 / (k_j * sigma_sq_j)
  w_tilde_j <- rep(w_tilde, k_j)
  
  mod_prelim <- lm.wfit(x = X, y = y, w = w_tilde_j)
  
  res <- residuals(mod_prelim) 
  
  # calculate weighted residual sum of squares
  QE <- sum(w_tilde_j * res^2)
  
  # split the design matrix by study
  X_j <- mat_split(X, cluster)
  
  # Create M tilde ----------------------------------------------------------
  M_tilde <- chol2inv(chol(crossprod(X, w_tilde_j * X)))
  p_j <- lapply(X_j, colSums)
  
  # trace products ----------------------------------------------------------
  
  # the first B for the trace product numerator
  w_over_k <- w_tilde / k_j
  B_num_all_1 <- Map(function(w, X) w * crossprod(X), w = w_over_k, X = X_j)
  B_num_1 <- reduce(B_num_all_1, `+`)
  
  # the second B for the trace product numerator
  XJX_j <- lapply(p_j, tcrossprod)
  B_num_all_2 <- Map(function(w, xjx) w * xjx, w = w_over_k, xjx = XJX_j)
  B_num_2 <- reduce(B_num_all_2, `+`)
  
  # the B for the trade product denominator
  B_den_all <- Map(function(w, xjx) w^2 * xjx, w = w_tilde, xjx = XJX_j)
  B_den <- reduce(B_den_all, `+`)
  
  num_minus <- m - (1 - rho) * trace_product(M_tilde, B_num_1) - rho * trace_product(M_tilde, B_num_2)
  den <- sum(k_j * w_tilde) - trace_product(M_tilde, B_den)
  
  tau_sq <- (QE - num_minus) / den
  tau_sq <- ifelse(tau_sq < 0, 0, tau_sq)  # added this like robu
  
  # CE weights
  w_j <- 1 / (k_j * (sigma_sq_j + tau_sq))
  w_ij <- rep(w_j, k_j)
  
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
  
  class(res) <- "handmade"
  
  if (!is.null(calc_vcov)) {
    res$vcov <- vcovCR(res, cluster = cluster, type = calc_vcov, inverse_var = TRUE)
  }
  
  return(res)
  
}

model.matrix.handmade <- function(object, ...) object$X

bread.handmade <- function(x, ...) {
  x$nobs * chol2inv(chol(crossprod(x$X, x$w_ij * x$X)))
}

update_robu <- function(mod, y, fix_tau_sq = FALSE, calc_CR0 = TRUE) {
  
  if (fix_tau_sq) {
    tau_sq <- mod$tau_sq
  } else {
    mod_prelim <- lm.wfit(x = mod$X, y = y, w = mod$w_tilde_j)
    res <- residuals(mod_prelim) 
    QE <- sum(mod$w_tilde_j * res^2)
    tau_sq <- (QE - mod$num_minus) / mod$den
    tau_sq <- ifelse(tau_sq < 0, 0, tau_sq)  # added this
  }
  
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
