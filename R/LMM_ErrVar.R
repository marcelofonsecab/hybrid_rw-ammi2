#' Calculate Error Variance from a Linear Mixed Model (LMM)
#'
#' This function computes the error variance using a linear mixed model (LMM) for a given dataset.
#' Depending on whether there is a single environment or a single genotype, the model is adjusted accordingly.
#'
#' @param data A data frame containing genotype-environment data.
#'
#' @return A numeric value representing the error variance.
#'
#' @examples
#' sim_data <- sim.amb(seed = 123)
#' rep_data <- Create_Replications(sim_data, reps = 2, sig = 1)
#' err_var <- LMM_Error_variance(rep_data)
#' err_var
#'
#' @importFrom lme4 lmer VarCorr
#' @export
LMM_Error_variance <- function(data) {
  quant.gen <- length(unique(data$gen))
  quant.env <- length(unique(data$env))

  if (quant.env == 1) {
    aux <- lme4::lmer(yield ~ rep + (1 | gen), data = data)
    tmp <- attr(VarCorr(aux), "sc")^2
  } else if (quant.gen == 1) {
    aux <- lme4::lmer(yield ~ rep + (1 | env), data = data)
    tmp <- attr(VarCorr(aux), "sc")^2
  }

  return(tmp)
}
