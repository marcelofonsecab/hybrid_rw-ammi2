#' Simulate Genotype-Environment Data
#'
#' This function simulates a new dataset for genotype-environment interaction (GEI) analysis.
#' It generates data based on a global mean, genotype effects, environment effects, and random
#' multiplicative effects via singular value decomposition.
#'
#' @param seed Optional seed for reproducibility. If \code{NULL}, a random seed is generated.
#' @param Ngen Number of genotypes. Default is 100.
#' @param Nenv Number of environments. Default is 8.
#' @param Ncomp Number of principal components to use for the multiplicative effects. Default is 2.
#' @param effectGlobal A numeric vector of length 2 (mean and standard deviation) for the global effect.
#'   Default is \code{c(mean = 15, sd = sqrt(3))}.
#' @param effectGen A numeric vector of length 2 (mean and standard deviation) for genotype effects.
#'   Default is \code{c(mean = 5, sd = 1)}.
#' @param effectEnv A numeric vector of length 2 (mean and standard deviation) for environment effects.
#'   Default is \code{c(mean = 8, sd = sqrt(2))}.
#' @param k A numeric vector of weights for each principal component. Default is \code{rep(1, Nenv)}.
#'
#' @return A data frame containing the simulated GEI data with the following columns:
#' \describe{
#'   \item{gen}{Genotype identifier (factor).}
#'   \item{env}{Environment identifier (factor).}
#'   \item{rep}{Replication number (factor, default is 1).}
#'   \item{yield}{Simulated response variable.}
#' }
#'
#' @examples
#' # Simulate data with default parameters
#' sim_data <- sim.amb(seed = 123)
#' head(sim_data)
#'
#' @export
sim.amb <- function(seed = NULL, Ngen = 100, Nenv = 8, Ncomp = 2,
                    effectGlobal = c(mean = 15, sd = sqrt(3)),
                    effectGen = c(mean = 5, sd = 1),
                    effectEnv = c(mean = 8, sd = sqrt(2)),
                    k = rep(1, Nenv)) {

  # Set a seed for reproducibility
  seed.aux <- ifelse(is.null(seed), sample(1:2^16, 1), seed)
  set.seed(seed.aux)

  # Generate overall global, genotype, and environment effects
  globalMean <- rnorm(1, effectGlobal[1], effectGlobal[2])
  alpha <- rnorm(Ngen, mean = effectGen[1], sd = effectGen[2])
  beta <- rnorm(Nenv, mean = effectEnv[1], sd = effectEnv[2])

  # Generate a random matrix and perform singular value decomposition (SVD)
  rand.mat <- matrix(runif(Ngen * Nenv, min = -0.5, max = 0.5), ncol = Nenv)
  rand.svd <- svd(rand.mat)

  # Calculate the basic response matrix using main effects
  simulated.amb <- matrix(1, nrow = Ngen, ncol = Nenv) * globalMean +
    alpha %*% t(rep(1, Nenv)) +
    rep(1, Ngen) %*% t(beta)

  # Add multiplicative effects using Ncomp principal components
  for (j in 1:Ncomp) {
    simulated.amb <- simulated.amb +
      ((rand.svd$u[, j] * rand.svd$d[j]) %*% t(rand.svd$v[, j])) * k[j]
  }

  # Organize the simulated data into a data frame
  aux <- data.frame(
    gen   = as.factor(rep(paste0("G", sprintf('%03d', 1:Ngen)), each = Nenv)),
    env   = as.factor(rep(paste0("E", 1:Nenv), times = Ngen)),
    rep   = factor(1),
    yield = as.vector(simulated.amb)
  )

  return(aux)
}
