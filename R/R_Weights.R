#' Calculate Model Weights for LMM and RLMM
#'
#' This function calculates various weighting schemes for the data based on error variances
#' obtained from LMM and RLMM analyses, as well as weights derived from robust linear models.
#'
#' @param data A data frame containing genotype-environment data.
#' @param errorvariances A list containing error variances for genotypes and environments,
#'   as computed by \code{Error_Var()}.
#'
#' @return A list containing sub-lists for LMM and RLMM weights with several weighting schemes.
#'
#' @examples
#' sim_data <- sim_amb(seed = 123)
#' rep_data <- Create_Replications(sim_data, reps = 2, sig = 1)
#' err_vars <- Error_Var(rep_data)
#' weights <- R_Weights(rep_data, errorvariances = err_vars)
#' names(weights)
#'
#' @export
R_Weights <- function(data, errorvariances) {
  Ngen <- nlevels(data$gen)
  Nenv <- nlevels(data$env)
  Nrep <- max(data$rep)

  rep.weights <- tapply(data[,"rep"], data[,c("gen", "env")], max) / Nrep

  lmm.env.weights.tmp <- 1 / errorvariances$EnvironmentLMM$Error_Variance
  lmm.Weights.env.tmp <- matrix(rep(lmm.env.weights.tmp / max(lmm.env.weights.tmp), Ngen),
                                ncol = Nenv, byrow = TRUE)
  lmm.Weights.env <- lmm.Weights.env.tmp * rep.weights

  lmm.gen.weight.tmp <- 1 / errorvariances$GenotypeLMM$Error_Variance
  lmm.Weights.gen.tmp <- matrix(rep(lmm.gen.weight.tmp / max(lmm.gen.weight.tmp), Nenv),
                                nrow = Ngen, byrow = FALSE)
  lmm.Weights.gen <- lmm.Weights.gen.tmp * rep.weights

  data.rlm <- transform_usable_data(data, median, "dataframe")
  model.rlm.tmp <- MASS::rlm(yield ~ gen + env, data = data.rlm)
  rlm.weights.tmp <- model.rlm.tmp$w
  Weights.rlm <- matrix(rlm.weights.tmp, ncol = Nenv, byrow = FALSE)

  lmm.Weights.gendotenv <- lmm.Weights.gen * lmm.Weights.env
  lmm.Weights.genplusenv <- ((lmm.Weights.gen + lmm.Weights.env) / 2)
  lmm.Weights.rlmdotenv <- Weights.rlm * lmm.Weights.env
  lmm.Weights.rlmdotgen <- Weights.rlm * lmm.Weights.gen
  lmm.Weights.rlmdotgendotenv <- Weights.rlm * lmm.Weights.gen * lmm.Weights.env
  lmm.Weights.rlmdotgenplusenv <- Weights.rlm * ((lmm.Weights.gen + lmm.Weights.env) / 2)

  rlmm.env.weights.tmp <- 1 / errorvariances$EnvironmentRLMM$Error_Variance
  rlmm.Weights.env.tmp <- matrix(rep(rlmm.env.weights.tmp / max(rlmm.env.weights.tmp), Ngen),
                                 ncol = Nenv, byrow = TRUE)
  rlmm.Weights.env <- rlmm.Weights.env.tmp * rep.weights

  rlmm.gen.weight.tmp <- 1 / errorvariances$GenotypeRLMM$Error_Variance
  rlmm.Weights.gen.tmp <- matrix(rep(rlmm.gen.weight.tmp / max(rlmm.gen.weight.tmp), Nenv),
                                 nrow = Ngen, byrow = FALSE)
  rlmm.Weights.gen <- rlmm.Weights.gen.tmp * rep.weights

  rlmm.Weights.gendotenv <- rlmm.Weights.gen * rlmm.Weights.env
  rlmm.Weights.genplusenv <- ((rlmm.Weights.gen + rlmm.Weights.env) / 2)
  rlmm.Weights.rlmdotenv <- Weights.rlm * rlmm.Weights.env
  rlmm.Weights.rlmdotgen <- Weights.rlm * rlmm.Weights.gen
  rlmm.Weights.rlmdotgendotenv <- Weights.rlm * rlmm.Weights.gen * rlmm.Weights.env
  rlmm.Weights.rlmdotgenplusenv <- Weights.rlm * ((rlmm.Weights.gen + rlmm.Weights.env) / 2)

  weights_list <- list(
    LMM = list(
      lmm.Weights.Env = lmm.Weights.env,
      lmm.Weights.Gen = lmm.Weights.gen,
      lmm.Weights.GendotEnv = lmm.Weights.gendotenv,
      lmm.Weights.GenplusEnv = lmm.Weights.genplusenv,
      Weights.Rlm = Weights.rlm,
      lmm.Weights.RlmdotEnv = lmm.Weights.rlmdotenv,
      lmm.Weights.RlmdotGen = lmm.Weights.rlmdotgen,
      lmm.Weights.RlmdotGendotEnv = lmm.Weights.rlmdotgendotenv,
      lmm.Weights.RlmdotGenplusEnv = lmm.Weights.rlmdotgenplusenv
    ),
    RLMM = list(
      rlmm.Weights.Env = rlmm.Weights.env,
      rlmm.Weights.Gen = rlmm.Weights.gen,
      rlmm.Weights.GendotEnv = rlmm.Weights.gendotenv,
      rlmm.Weights.GenplusEnv = rlmm.Weights.genplusenv,
      Weights.Rlm = Weights.rlm,
      rlmm.Weights.RlmdotEnv = rlmm.Weights.rlmdotenv,
      rlmm.Weights.RlmdotGen = rlmm.Weights.rlmdotgen,
      rlmm.Weights.RlmdotGendotEnv = rlmm.Weights.rlmdotgendotenv,
      rlmm.Weights.RlmdotGenplusEnv = rlmm.Weights.rlmdotgenplusenv
    )
  )

  return(weights_list)
}
