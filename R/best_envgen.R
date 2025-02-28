#' Select the Best Genotype for Each Environment
#'
#' This function identifies the best (or second best, etc.) genotype for every environment
#' based on projections derived from SVD results.
#'
#' @param x An object containing an SVD result (with elements \code{SVD}).
#' @param quant An integer indicating which best genotype to select (1 for best, 2 for second best, etc.). Default is 1.
#' @param gen_names0 Optional original genotype names. Default is \code{NULL} (will generate numeric names).
#' @param env_names Optional environment names. Default is \code{NULL} (will generate default names).
#' @param colors Logical; if \code{TRUE}, returns colored output. Default is \code{FALSE}.
#'
#' @return A matrix (or character matrix if \code{colors} is TRUE) indicating the best genotype for each environment.
#'
#' @examples
#' sim_data <- sim.amb(seed = 123)
#' ammi_res <- ammi.model(sim_data, Ncomp = 2)
#' best_geno <- best_envgen(ammi_res, quant = 1)
#' best_geno
#'
#' @export
best_envgen <- function(x, quant = 1, gen_names0 = NULL, env_names = NULL, colors = FALSE) {
  Ncomp <- 2
  diag_vals <- if (is.matrix(x$SVD$d)) x$SVD$d else diag(x$SVD$d)
  scores <- x$SVD$u[, 1:Ncomp] %*% (diag_vals[1:Ncomp, 1:Ncomp] ^ 0.5)
  loadings <- x$SVD$v[, 1:Ncomp] %*% (diag_vals[1:Ncomp, 1:Ncomp] ^ 0.5)

  gen_names <- paste0(1:nrow(scores))
  if (is.null(env_names)) {
    env_names <- paste0("ENV", 1:nrow(loadings))
  }

  rownames(scores) <- gen_names
  rownames(loadings) <- env_names

  df_aux <- matrix("0", nrow = length(env_names), ncol = 1)
  rownames(df_aux) <- env_names
  colnames(df_aux) <- paste0(quant)

  data_aux <- df_aux

  for (env in env_names) {
    line_projection <- c(-loadings[env, 2],
                         loadings[env, 1],
                         loadings[env, 1],
                         loadings[env, 2])
    line_projection <- matrix(line_projection, nrow = 2)
    gen_pos <- NULL
    for (i in 1:length(gen_names)) {
      vec_coord <- matrix(c(0, loadings[env, 1] * scores[i, 1] + loadings[env, 2] * scores[i, 2]), ncol = 1)
      solve_matrix <- solve(line_projection, vec_coord)
      gen_pos <- rbind(gen_pos, t(solve_matrix))
    }

    values <- NULL
    for (i in 1:length(gen_names)) {
      values <- c(values, sum(loadings[env, ] * gen_pos[i, ]) / sqrt(sum(loadings[env, ]^2)))
    }

    aux <- data.frame(values, gen_names)
    best_quant <- aux[order(aux[,"values"], decreasing = TRUE), 2][quant]
    data_aux[env, ] <- best_quant
  }

  if (is.null(gen_names0)) {
    gen_names0 <- paste0(1:nrow(scores))
  }
  gen_names0 <- gen_names0[as.numeric(data_aux)]

  if (colors) {
    # Placeholder for color conversion (implementation-dependent)
    data_aux <- paste0(gen_names0, " [colored]")
  } else {
    data_aux[, 1] <- gen_names0
  }

  return(data_aux)
}
