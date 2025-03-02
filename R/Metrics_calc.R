#' Calculate GEI Metrics for Different Models
#'
#' This function calculates various metrics for evaluating GEI models, such as
#' mean percentage of explained variance (MPEV), mean squared error on singular values (MSE.Singvals),
#' and other related metrics based on the SVD of residuals.
#'
#' @param allSVDs A list of SVD results from models (e.g., output from \code{All_SVDS()}).
#' @param Ncomp Number of SVD components used in the models. Default is 2.
#'
#' @return A list of metrics for each model.
#'
#' @examples
#' # Assuming 'models' is the output from All_SVDS(...)
#' metrics <- GEI_Metrics(models, Ncomp = 2)
#' metrics$AMMI
#'
#' @export
GEI_Metrics <- function(allSVDs, Ncomp = 2) {
  Metrics <- list()
  realresidual <- allSVDs$RealSVD$residual
  realdf <- allSVDs$RealSVD$matrix.fitted + realresidual

  ammiresidual <- allSVDs$AMMI$reduced.SVD
  ammidf <- allSVDs$AMMI$matrix.fitted + ammiresidual

  Ngen <- nrow(realresidual)
  Nenv <- ncol(realresidual)

  real.singvals <- svd(realresidual)$d
  ammi.singvals <- svd(ammiresidual)$d
  aux.mse <- (ammi.singvals[1:Ncomp] - real.singvals[1:Ncomp])^2

  MPEV <- (sum(ammi.singvals[1:Ncomp]^2) / sum(real.singvals^2)) * 100
  MSE.Singvals <- round(aux.mse, digits = 5)
  MtSPE <- round(sum(sort((ammidf - realdf)^2)[1:(Ngen * Nenv * 0.9)]) / (Ngen * Nenv * 0.9), 5)
  I <- svd(realresidual)$v[, 1:Ncomp]
  P <- svd(ammiresidual)$v[, 1:Ncomp]
  small.eigen <- min(eigen(t(I) %*% P %*% t(P) %*% I)$values)
  Maxsub <- round(acos(sqrt(small.eigen)) / (pi / 2), 5)

  Metrics[["AMMI"]] <- list("MPEV" = MPEV,
                            "MSE.Singvals1" = MSE.Singvals[1],
                            "MSE.Singvals2" = MSE.Singvals[2],
                            "MtSPE" = MtSPE,
                            "Maxsub" = Maxsub)

  if (length(allSVDs) > 2) {
    for (i in 3:length(allSVDs)) {
      for (k in 1:length(allSVDs[[i]])) {
        aux.matrix <- allSVDs[[i]][[k]]$reduced.SVD
        aux.svd <- svd(aux.matrix)
        aux.singvals <- aux.svd$d
        aux.fitted <- allSVDs[[i]][[k]]$matrix.fitted
        aux.est <- aux.fitted + aux.matrix

        MPEV <- (sum(aux.singvals[1:Ncomp]^2) / sum(real.singvals^2)) * 100
        MSE.Singvals <- round((aux.singvals[1:Ncomp] - real.singvals[1:Ncomp])^2, digits = 5)
        MtSPE <- round(sum(sort((aux.est - realdf)^2)[1:(Ngen * Nenv * 0.9)]) / (Ngen * Nenv * 0.9), 5)
        I <- svd(realresidual)$v[, 1:Ncomp]
        P <- aux.svd$v[, 1:Ncomp]
        small.eigen <- min(eigen(t(I) %*% P %*% t(P) %*% I)$values)
        Maxsub <- round(acos(sqrt(small.eigen)) / (pi / 2), 5)

        Metrics[[paste0(names(allSVDs[[i]])[k])]] <- list("MPEV" = MPEV,
                                                          "MSE.Singvals1" = MSE.Singvals[1],
                                                          "MSE.Singvals2" = MSE.Singvals[2],
                                                          "MtSPE" = MtSPE,
                                                          "Maxsub" = Maxsub)
      }
    }
  }

  return(Metrics)
}
