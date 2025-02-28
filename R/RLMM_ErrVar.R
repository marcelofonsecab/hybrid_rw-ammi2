#' Calculate Error Variance from a Robust Linear Mixed Model (RLMM)
#'
#' This function computes the error variance using a robust linear mixed model (RLMM)
#' for a given dataset. The model adjusts based on whether there is a single environment or genotype.
#'
#' @param data A data frame containing genotype-environment data.
#'
#' @return A numeric value representing the error variance.
#'
#' @examples
#' sim_data <- sim.amb(seed = 123)
#' rep_data <- Create.Replications(sim_data, reps = 2, sig = 1)
#' err_var <- RLMM.Error_variance(rep_data)
#' err_var
#'
#' @export
RLMM.Error_variance <- function(data) {
  quant.gen <- length(unique(data$gen))
  quant.env <- length(unique(data$env))

  if (quant.env == 1) {
    aux <- robustlmm::rlmer(yield ~ rep + (1 | gen), data = data)
    tmp <- attr(VarCorr(aux), "sc")^2
  } else if (quant.gen == 1) {
    aux <- robustlmm::rlmer(yield ~ rep + (1 | env), data = data)
    tmp <- attr(VarCorr(aux), "sc")^2
  }

  return(tmp)
}
