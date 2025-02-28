#' Fit a Robust Weighted AMMI (RW-AMMI) Model
#'
#' This function fits a robust weighted AMMI model by applying robust linear regression
#' with weights and then performing a weighted (or robust) SVD on the residuals.
#'
#' @param dataframe A data frame with genotype-environment data.
#' @param weight A numeric vector of weights. Must be provided.
#' @param Ncomp Number of SVD components to retain. Default is 2.
#' @param rSVD Logical; if \code{TRUE}, performs robust SVD. Default is \code{FALSE}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{residual}{The residual matrix from the robust weighted regression.}
#'     \item{matrix.fitted}{The matrix of fitted values from the regression.}
#'     \item{reduced.SVD}{The weighted SVD reconstruction.}
#'     \item{SVD}{The SVD object of the reduced matrix.}
#'   }
#'
#' @examples
#' sim_data <- sim.amb(seed = 123)
#' rwammi_res <- rwammi.model(sim_data, weight = rep(1, nrow(sim_data)), Ncomp = 2, rSVD = FALSE)
#' names(rwammi_res)
#'
#' @export
rwammi.model <- function(dataframe, weight = NULL, Ncomp = 2, rSVD = FALSE) {
  data <- transform_usable_data(dataframe, median, "dataframe")
  Ngen <- nlevels(data$gen)
  Nenv <- nlevels(data$env)
  weight <- c(weight)

  model <- MASS::rlm(yield ~ gen + env, data = data, w = weight)
  residual.matrix <- matrix(model$residuals, ncol = Nenv, nrow = Ngen)
  fitted.values <- matrix(model$fitted.values, ncol = Nenv, nrow = Ngen)

  SVD.aux <- WeightedSvdRes(residual.matrix, weight = weight, Ncomp = Ncomp, robust = rSVD)
  colnames(SVD.aux) <- colnames(residual.matrix)
  rownames(SVD.aux) <- rownames(residual.matrix)
  SVD.red <- SVD.aux
  SVD.fin <- svd(SVD.red)

  return(list(residual = residual.matrix,
              matrix.fitted = fitted.values,
              reduced.SVD = SVD.red,
              SVD = SVD.fin))
}
