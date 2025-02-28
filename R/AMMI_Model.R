#' Fit an AMMI Model
#'
#' This function fits the Additive Main effects and Multiplicative Interaction (AMMI) model.
#' It calculates fitted values, residuals, and performs singular value decomposition on the residuals.
#'
#' @param dataframe A data frame with genotype-environment data.
#' @param Ncomp Number of components to retain in the SVD. Default is 2.
#'
#' @return A list containing:
#'   \describe{
#'     \item{estimated.data}{The fitted data including main and multiplicative effects.}
#'     \item{anovatable}{The ANOVA table of the fitted linear model.}
#'     \item{SVD}{The SVD object from the residual matrix.}
#'     \item{residual}{The residual matrix from the fitted model.}
#'     \item{matrix.fitted}{The matrix of fitted values from the linear model.}
#'     \item{reduced.SVD}{The SVD reconstruction using only the first \code{Ncomp} components.}
#'   }
#'
#' @examples
#' sim_data <- sim.amb(seed = 123)
#' ammi_res <- ammi.model(sim_data, Ncomp = 2)
#' names(ammi_res)
#'
#' @export
ammi.model <- function(dataframe, Ncomp = 2) {
  df.mean <- transform_usable_data(dataframe, mean, type = "dataframe")
  Ngen <- nlevels(df.mean$gen)
  Nenv <- nlevels(df.mean$env)
  env.names <- levels(df.mean$env)
  gen.names <- levels(df.mean$gen)

  model.GEI <- lm(yield ~ gen + env, data = df.mean)
  matrix.GEI <- matrix(residuals(model.GEI), ncol = Nenv, nrow = Ngen)
  svd.GEI <- svd(matrix.GEI)
  svd.u <- svd.GEI$u
  svd.d <- diag(svd.GEI$d)
  svd.v <- svd.GEI$v

  aux <- matrix(model.GEI$fitted.values, ncol = Nenv)
  if (Ncomp >= 1) {
    for (i in 1:Ncomp) {
      aux <- aux + (svd.u[, i] * svd.d[i, i]) %*% t(svd.v[, i])
    }
  }

  colnames(aux) <- env.names
  rownames(aux) <- gen.names
  matrix.fitted <- matrix(model.GEI$fitted.values, ncol = Nenv)

  SVD.red <- if (Ncomp == 1) 0 else svd.u[, 1:Ncomp] %*% (svd.d[1:Ncomp, 1:Ncomp]) %*% t(svd.v[, 1:Ncomp])
  colnames(SVD.red) <- env.names
  rownames(SVD.red) <- gen.names

  return(list(estimated.data = aux,
              anovatable = anova(model.GEI),
              SVD = svd.GEI,
              residual = matrix.GEI,
              matrix.fitted = matrix.fitted,
              reduced.SVD = SVD.red))
}
