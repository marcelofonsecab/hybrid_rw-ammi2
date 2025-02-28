#' Create a Biplot from SVD Results
#'
#' This function creates a biplot using scores and loadings extracted from an SVD object.
#' The scores and loadings are scaled by the square roots of the singular values.
#'
#' @param SVD An SVD object (list with elements \code{u}, \code{d}, and \code{v}).
#' @param env.names Optional character vector for environment names. Defaults to generated names.
#' @param gen.names Optional character vector for genotype names. Defaults to numeric indices.
#' @param nms A character vector of length 2 for labeling the axes. Default is \code{c("IPC1", "IPC2")}.
#' @param lims Logical; if \code{FALSE}, axis limits are automatically calculated. Default is \code{FALSE}.
#' @param lims.U Optional numeric vector for x-axis limits.
#' @param lims.V Optional numeric vector for y-axis limits.
#'
#' @return Produces a biplot on the current graphic device.
#'
#' @examples
#' sim_data <- sim.amb(seed = 123)
#' ammi_res <- ammi.model(sim_data, Ncomp = 2)
#' BiplotCreation(ammi_res$SVD, env.names = levels(sim_data$env), gen.names = levels(sim_data$gen))
#'
#' @export
BiplotCreation <- function(SVD, env.names = NULL, gen.names = NULL,
                           nms = c("IPC1", "IPC2"), lims = FALSE,
                           lims.U, lims.V) {
  Ncomp <- 2
  diag_vals <- if (is.matrix(SVD$d)) SVD$d else diag(SVD$d)
  scores <- SVD$u[, 1:Ncomp] %*% (diag_vals[1:Ncomp, 1:Ncomp] ^ 0.5)
  loadings <- SVD$v[, 1:Ncomp] %*% (diag_vals[1:Ncomp, 1:Ncomp] ^ 0.5)

  if (is.null(env.names)) {
    env.names <- paste0("Env ", 1:nrow(loadings))
  }
  if (is.null(gen.names)) {
    gen.names <- 1:nrow(scores)
  }

  expand <- 1.01
  if (!lims) {
    lims.U <- c(min(loadings[, 1] * expand, scores[, 1] * expand),
                max(loadings[, 1] * expand, scores[, 1] * expand))
    lims.V <- c(min(loadings[, 2] * expand, scores[, 2] * expand),
                max(loadings[, 2] * expand, scores[, 2] * expand))
  }

  par(pty = "s")
  plot(scores, ylim = lims.V, xlim = lims.U, type = "n",
       xlab = nms[1], ylab = nms[2])
  abline(h = 0)
  abline(v = 0)
  points(scores, cex = 1.2, pch = 21, bg = "blue", col = "black")

  par(new = TRUE)
  plot(0, 0, ylim = lims.V, xlim = lims.U, type = "n",
       axes = FALSE, xlab = "", ylab = "")
  arrows(0, 0, loadings[, 1], loadings[, 2], length = 0.1, col = "forestgreen", lwd = 2)
  text(loadings[, 1] * 1.01, loadings[, 2] * 1.01,
       labels = env.names, col = "darkgreen", cex = 1.4)
}
