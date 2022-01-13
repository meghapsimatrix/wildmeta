# Constrain Predictors ----------------------------------------------------

constrain_predictors <- function(Xmat, Cmat) {

  q <- nrow(Cmat)
  p <- ncol(Cmat)
  if (ncol(Xmat) != ncol(Cmat)) stop("Constraint matrix must have same number of columns as predictor matrix.")

  XtX_inv <- chol2inv(chol(crossprod(Xmat)))
  Cnull <- diag(nrow = p) - XtX_inv %*% t(Cmat) %*% chol2inv(chol(Cmat %*% XtX_inv %*% t(Cmat))) %*% Cmat
  Xnull <- qr.X(qr(Xmat %*% Cnull), ncol = p - q)

  return(Xnull)

}



# Auxiliary distribution --------------------------------------------------

#' @importFrom stats runif
#' @importFrom stats rnorm

wild_wts <- function(auxiliary_dist, n_clusters) {

  auxiliary_dist <- match.arg(auxiliary_dist,
                              c("Rademacher","Mammen","Webb six",
                                "uniform","standard normal"),
                              several.ok = FALSE)

  switch(
    auxiliary_dist,
    Rademacher = sample(c(-1, 1),
                        size = n_clusters,
                        replace = TRUE),
    Mammen = sample(c(-(sqrt(5) - 1)/2, (sqrt(5) + 1)/2),
                    size = n_clusters,
                    replace = TRUE,
                    prob = c((sqrt(5) + 1) /(2 * sqrt(5)), (sqrt(5) - 1) /(2 * sqrt(5)))),
    `Webb six` = sample(c(-sqrt(3:1),sqrt(1:3)) / sqrt(2),
                        size = n_clusters,
                        replace = TRUE),
    uniform = stats::runif(n = n_clusters,
                           min = -sqrt(3),
                           max = sqrt(3)),
    `standard normal` = stats::rnorm(n = n_clusters)
  )

}
