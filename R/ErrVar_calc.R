#' Calculate Error Variances by Group
#'
#' This function calculates the error variances for each genotype and environment using
#' both a linear mixed model (LMM) and a robust linear mixed model (RLMM). Optionally,
#' parallel processing can be used via the \code{cluster} argument.
#'
#' @param data A data frame containing genotype-environment data.
#' @param cluster Optional cluster object for parallel processing. Default is \code{NULL}.
#'
#' @return A list with four elements containing error variances:
#'   \describe{
#'     \item{GenotypeLMM}{A tibble with error variances for each genotype (LMM).}
#'     \item{EnvironmentLMM}{A tibble with error variances for each environment (LMM).}
#'     \item{GenotypeRLMM}{A tibble with error variances for each genotype (RLMM).}
#'     \item{EnvironmentRLMM}{A tibble with error variances for each environment (RLMM).}
#'   }
#'
#' @examples
#' sim_data <- sim.amb(seed = 123)
#' rep_data <- Create.Replications(sim_data, reps = 2, sig = 1)
#' err_vars <- Error_Var(rep_data)
#' err_vars
#'
#' @export
Error_Var <- function(data, cluster = NULL) {
  tmp1 <- data %>%
    dplyr::group_by(gen) %>%
    dplyr::group_split() %>%
    pbapply::pblapply(LMM.Error_variance, cl = cluster) %>%
    unlist() %>%
    tibble::tibble(Group = unique(data$gen), Error_Variance = .)

  tmp2 <- data %>%
    dplyr::group_by(env) %>%
    dplyr::group_split() %>%
    pbapply::pblapply(LMM.Error_variance, cl = cluster) %>%
    unlist() %>%
    tibble::tibble(Group = unique(data$env), Error_Variance = .)

  tmp3 <- data %>%
    dplyr::group_by(gen) %>%
    dplyr::group_split() %>%
    pbapply::pblapply(RLMM.Error_variance, cl = cluster) %>%
    unlist() %>%
    tibble::tibble(Group = unique(data$gen), Error_Variance = .)

  tmp4 <- data %>%
    dplyr::group_by(env) %>%
    dplyr::group_split() %>%
    pbapply::pblapply(RLMM.Error_variance, cl = cluster) %>%
    unlist() %>%
    tibble::tibble(Group = unique(data$env), Error_Variance = .)

  tmp <- list(GenotypeLMM = tmp1, EnvironmentLMM = tmp2,
              GenotypeRLMM = tmp3, EnvironmentRLMM = tmp4)
  return(tmp)
}
