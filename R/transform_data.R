#' Transform Data by Applying a Function
#'
#' This function transforms genotype-environment data either into a summarized data frame
#' or into a matrix by applying a specified function (e.g., \code{mean} or \code{median})
#' on the yield values grouped by genotype and environment.
#'
#' @param dataframe A data frame containing genotype-environment data.
#' @param func A function to apply on the \code{yield} column (e.g., \code{mean} or \code{median}).
#' @param type A character string specifying the output type: either \code{"dataframe"} or \code{"matrix"}.
#'   Default is \code{"dataframe"}.
#'
#' @return Either a summarized data frame or a matrix of transformed yield values.
#'
#' @examples
#' sim_data <- sim.amb(seed = 123)
#' trans_df <- transform_usable_data(sim_data, mean, type = "dataframe")
#' str(trans_df)
#'
#' @export
transform_usable_data <- function(dataframe, func, type = c("dataframe", "matrix")) {
  type_aux <- c("dataframe", "matrix")
  match.arg(type, type_aux)

  if (type == "dataframe") {
    df <- dataframe %>%
      dplyr::group_by(gen, env) %>%
      dplyr::summarise(yield = func(yield)) %>%
      dplyr::arrange(env)
  } else if (type == "matrix") {
    df <- tapply(dataframe[, "yield"], dataframe[, c("gen", "env")], func)
  }

  return(df)
}
