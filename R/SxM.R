#' Steptoe x Morex Genotype–Environment Data
#'
#' A dataset containing genotype and environment observations from the
#' Steptoe x Morex mapping population. This classic dataset is widely used
#' in studies of genotype–environment interaction and plant breeding.
#'
#' @format A data frame with \eqn{2100} rows and \eqn{4} columns. Each row corresponds to a measurement for a genotype in a specific environment.
#'
#' \describe{
#'   \item{gen}{A factor indicating the genotype identifier.}
#'   \item{env}{A factor indicating the location where the measurement was taken.}
#'   \item{rep}{A factor indicating the replication number for each observation. It distinguishes repeated measurements or trials for the same genotype–environment combination.}
#'   \item{yield}{A numeric vector containing the observed trait values.}
#' }
#'
#' @details
#' The data originate from the Steptoe x Morex population, a well‐characterized mapping population used in
#' cereal genetics research. For more details on the experimental design and data collection, please visit:
#' \url{https://wheat.pw.usda.gov/ggpages/SxM/}
#'
#' @source \url{https://wheat.pw.usda.gov/ggpages/SxM/}
#'
"SxM"
