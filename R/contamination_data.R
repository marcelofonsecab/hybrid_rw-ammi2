#' Contaminate Genotype-Environment Data
#'
#' This function introduces contamination into the first replication of a genotype-environment dataset.
#' Different types of contamination can be applied: shifting values, creating point mass contamination,
#' or inflating variances.
#'
#' @param data A data frame containing the original genotype-environment data.
#' @param percentage Percentage of data to contaminate. If greater than 1, it is treated as a percentage.
#'   Default is 5.
#' @param seed Seed value for reproducibility. Default is 1.
#' @param type Type of contamination. Must be one of \code{"shift"}, \code{"pointmass"}, or \code{"varinflated"}.
#'   Default is \code{"shift"}.
#' @param k Shift factor for contamination if \code{type} is "shift". Default is 7.
#' @param c Scaling factor for contamination if \code{type} is "pointmass" or "varinflated". Default is 10.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{data.not_contaminated}{The original data.}
#'     \item{data.contaminated}{The contaminated dataset.}
#'   }
#'
#' @examples
#' sim_data <- sim_amb(seed = 123)
#' contaminated <- Data_Contamination(sim_data, percentage = 5, seed = 1, type = "shift")
#' str(contaminated)
#'
#' @export
Data_Contamination <- function(data, percentage = 5, seed = 1,
                               type = "shift", k = 7, c = 10) {
  set.seed(seed)
  types_aux <- c("shift", "pointmass", "varinflated")
  match.arg(type, types_aux)
  percentage_aux <- ifelse(percentage > 1, percentage / 100, percentage)
  df_aux <- data[data$rep == 1, ]
  matrix_aux <- tapply(data$yield, data[, c("gen", "env")], mean)
  Means_Deviations <- cbind(mean = apply(matrix_aux, 2, mean),
                            deviation = apply(matrix_aux, 2, sd))

  Nenv <- length(levels(df_aux$env))
  Ngen <- length(levels(df_aux$gen))
  Nrow <- nrow(df_aux)

  sample_aux <- sample(1:Nrow, Nrow * percentage_aux, replace = FALSE)
  df_tmp <- df_aux[sample_aux, ]

  for (i in 1:Nenv) {
    sample_tmp <- which(df_tmp$env == rownames(Means_Deviations)[i])
    df_tmp[sample_tmp, "yield"] <- rnorm(length(sample_tmp),
                                         mean = Means_Deviations[i, 1] +
                                           ifelse(type == "shift", k, 0) * Means_Deviations[i, 2],
                                         sd = Means_Deviations[i, 2] / ifelse(type %in% c("pointmass", "varinflated"), sqrt(c), 1))
  }

  df_aux[sample_aux, ] <- df_tmp
  df_aux <- rbind(df_aux, data[data$rep != 1, ])
  data.full <- list(data.not_contaminated = data, data.contaminated = df_aux)

  return(data.full)
}
