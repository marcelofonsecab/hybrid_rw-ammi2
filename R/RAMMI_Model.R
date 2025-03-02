#' Fit a Robust AMMI (R-AMMI) Model
#'
#' This function fits a robust AMMI model using robust linear regression.
#' It applies the regression to median-transformed data and then performs a robust SVD on the residuals.
#'
#' @param dataframe A data frame with genotype-environment data.
#' @param weight A numeric vector of weights. If \code{NULL}, the regression is run without weights.
#' @param Ncomp Number of SVD components to retain. Default is 2.
#'
#' @return A list containing:
#'   \describe{
#'     \item{residual}{The residual matrix from the robust regression.}
#'     \item{matrix.fitted}{The matrix of fitted values from the robust regression.}
#'     \item{reduced.SVD}{The robust SVD reconstruction using \code{Ncomp} components.}
#'     \item{SVD}{The robust SVD object.}
#'   }
#'
#' @examples
#' sim_data <- sim_amb(seed = 123)
#' rammi_res <- rammi_model(sim_data, weight = rep(1, nrow(sim_data)), Ncomp = 2)
#' names(rammi_res)
#'
#' @importFrom MASS rlm
#' @importFrom rsvddpd rSVDdpd
#' @export
rammi_model <- function(dataframe, weight = NULL, Ncomp = 2) {
  data <- transform_usable_data(dataframe, median, "dataframe")
  Ngen <- nlevels(data$gen)
  Nenv <- nlevels(data$env)
  weight <- c(weight)

  if (is.null(weight)) {
    model <- MASS::rlm(yield ~ gen + env, data = data)
  } else {
    model <- MASS::rlm(yield ~ gen + env, data = data, w = weight)
  }

  residual.matrix <- matrix(model$residuals, ncol = Nenv, nrow = Ngen)
  fitted.values <- matrix(model$fitted.values, ncol = Nenv, nrow = Ngen)
  SVD <- rsvddpd::rSVDdpd(residual.matrix, alpha = 0.3, maxiter = 100,
                          initu = svd(residual.matrix)$u,
                          initv = svd(residual.matrix)$v,
                          eps = 1e-4)
  SVD.red <- SVD$u[, 1:Ncomp] %*% diag(SVD$d[1:Ncomp]) %*% t(SVD$v[, 1:Ncomp])

  return(list(residual = residual.matrix,
              matrix.fitted = fitted.values,
              reduced.SVD = SVD.red,
              SVD = SVD))
}
