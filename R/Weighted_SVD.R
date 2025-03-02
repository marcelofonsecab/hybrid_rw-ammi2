#' Perform Weighted (or Robust) Singular Value Decomposition
#'
#' This function performs a weighted singular value decomposition (SVD) on a data matrix.
#' Optionally, a robust SVD method (\code{rSVDdpd}) can be used.
#'
#' @param data A numeric matrix representing the data.
#' @param weight A numeric vector representing weights for the SVD. Default is \code{NULL}.
#' @param Ncomp The number of components to retain. Default is 2.
#' @param robust Logical; if \code{TRUE}, performs robust SVD using \code{rSVDdpd}. Default is \code{FALSE}.
#'
#' @return A matrix representing the reconstructed data from the weighted (or robust) SVD.
#'
#' @examples
#' # Create a random matrix and perform weighted SVD with 2 components
#' mat <- matrix(rnorm(100), nrow = 10)
#' W <- rep(1, length(mat))
#' svd_res <- WeightedSvdRes(mat, weight = W, Ncomp = 2)
#' dim(svd_res)
#'
#' @importFrom rsvddpd rSVDdpd
#' @export
WeightedSvdRes <- function(data, weight = NULL, Ncomp = 2, robust = FALSE) {
  Ngen <- nrow(data)
  Nenv <- ncol(data)
  W <- c(weight)
  D <- data
  X <- matrix(0, Ngen, Nenv)
  aux <- matrix(1, Ngen, Nenv)
  Xold <- Inf * aux
  Err <- Inf
  eps <- 1e-5

  if (!robust) {
    wsvd <- svd(W * D + (1 - W) * X)
  } else {
    wsvd <- rsvddpd::rSVDdpd(W * D + (1 - W) * X, alpha = 0.3, maxiter = 100,
                             initu = svd(W * D + (1 - W) * X)$u,
                             initv = svd(W * D + (1 - W) * X)$v,
                             eps = 1e-4)
  }
  U <- wsvd$u
  S <- diag(wsvd$d)
  V <- wsvd$v

  if (Ncomp + 1 < length(wsvd$d)) {
    S[(Ncomp + 1):length(wsvd$d), (Ncomp + 1):length(wsvd$d)] <- 0
  }

  X <- U %*% S %*% t(V)

  m <- 0
  while (m < 100) {
    m <- m + 1
    Xold <- X
    if (!robust) {
      wsvd <- svd(W * D + (1 - W) * X)
    } else {
      wsvd <- rsvddpd::rSVDdpd(W * D + (1 - W) * X, alpha = 0.3, maxiter = 100,
                               initu = svd(W * D + (1 - W) * X)$u,
                               initv = svd(W * D + (1 - W) * X)$v,
                               eps = 1e-4)
    }
    U <- wsvd$u
    S <- diag(wsvd$d)
    V <- wsvd$v
    if (Ncomp + 1 < length(wsvd$d)) {
      S[(Ncomp + 1):length(wsvd$d), (Ncomp + 1):length(wsvd$d)] <- 0
    }
    X <- U %*% S %*% t(V)
    Err <- sum((X - Xold)^2)
    if (Err < eps) break
  }

  return(X)
}
