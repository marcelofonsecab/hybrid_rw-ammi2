#' Fit a Weighted AMMI (W-AMMI) Model
#'
#' This function fits a weighted AMMI model by performing weighted linear regression
#' and then applying a weighted SVD on the residuals.
#'
#' @param data A data frame with genotype-environment data.
#' @param weight A numeric vector of weights. Must be provided.
#' @param Ncomp Number of SVD components to retain. Default is 2.
#'
#' @return A list containing:
#'   \describe{
#'     \item{dataframe}{The transformed data with an added column of weighted residuals.}
#'     \item{residual}{The residual matrix from the weighted regression.}
#'     \item{matrix.fitted}{The matrix of fitted values from the weighted regression.}
#'     \item{reduced.SVD}{The weighted SVD reconstruction.}
#'     \item{SVD}{The SVD object of the weighted residual matrix.}
#'   }
#'
#' @examples
#' sim_data <- sim.amb(seed = 123)
#' # Assume weight vector is defined appropriately
#' weight_vec <- rep(1, nrow(sim_data))
#' wammi_res <- wammi.model(sim_data, weight = weight_vec, Ncomp = 2)
#' names(wammi_res)
#'
#' @importFrom reshape2 melt
#' @export
wammi_model <- function(data, weight = NULL, Ncomp = 2) {
  dataframe <- transform_usable_data(data, mean, type = "dataframe")
  Nenv <- nlevels(dataframe$env)

  if (is.null(weight)) {
    stop("This function requires a weight vector")
  }

  weight <- c(weight)
  weighted.lm <- lm(yield ~ gen + env, weights = weight, data = dataframe)
  fitted.values <- matrix(weighted.lm$fitted.values, ncol = Nenv)
  residual.matrix <- matrix(weighted.lm$residuals, ncol = Nenv)
  colnames(residual.matrix) <- levels(dataframe$env)
  rownames(residual.matrix) <- levels(dataframe$gen)

  SVD.aux <- WeightedSvdRes(residual.matrix, weight = weight, Ncomp = Ncomp)
  colnames(SVD.aux) <- colnames(residual.matrix)
  rownames(SVD.aux) <- rownames(residual.matrix)

  SVD.red <- SVD.aux
  SVD.fin <- svd(SVD.red)

  dataframe <- cbind(dataframe, "W.residuals" = reshape2::melt(residual.matrix)$value)

  return(list(dataframe = dataframe,
              residual = residual.matrix,
              matrix.fitted = fitted.values,
              reduced.SVD = SVD.red,
              SVD = SVD.fin))
}
