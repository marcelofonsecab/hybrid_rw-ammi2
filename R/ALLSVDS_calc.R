#' Aggregate SVD Results from Multiple Models
#'
#' This function fits several models (AMMI, W-AMMI, R-AMMI, and RW-AMMI) to the provided data,
#' and aggregates their singular value decomposition (SVD) results.
#'
#' @param data_not_contaminated The original (non-contaminated) data frame.
#' @param data_contaminated The contaminated data frame.
#' @param weights A list of weight vectors for the models.
#' @param Ncomp Number of SVD components to retain. Default is 2.
#' @param rSVD Logical; if \code{TRUE}, use robust SVD in RW-AMMI. Default is \code{FALSE}.
#'
#' @return A list containing the fitted models for:
#'   \describe{
#'     \item{RealSVD}{AMMI model fitted to the original data.}
#'     \item{AMMI}{AMMI model fitted to the contaminated data.}
#'     \item{WAMMI}{List of W-AMMI models for different weight schemes.}
#'     \item{RAMMI}{List of R-AMMI models for different weight schemes.}
#'     \item{RWAMMI}{List of RW-AMMI models for different weight schemes.}
#'   }
#'
#' @examples
#' sim_data <- sim_amb(seed = 123)
#' rep_data <- Create_Replications(sim_data, reps = 2, sig = 1)
#' contaminated <- Data_Contamination(rep_data, percentage = 5, seed = 1, type = "shift")
#' # Assume weights is a predefined list of weight vectors
#' models <- All_SVDS(data_not_contaminated = rep_data,
#'                    data_contaminated = contaminated$data.contaminated,
#'                    weights = list(LMM = list(w1 = rep(1, nrow(rep_data))),
#'                                   RLMM = list(w1 = rep(1, nrow(rep_data)))),
#'                    Ncomp = 2, rSVD = FALSE)
#' names(models)
#'
#' @export
All_SVDS <- function(data_not_contaminated, data_contaminated, weights, Ncomp = 2, rSVD = FALSE) {
  WAMMI <- RAMMI <- RWAMMI <- list()
  Estimated.real <- ammi_model(data_not_contaminated, Ncomp = Ncomp)
  AMMI.model <- ammi_model(data_contaminated, Ncomp = Ncomp)

  for (i in 1:length(weights$LMM)) {
    message(paste0(i, "/", length(weights$LMM)))
    WAMMI[[paste0("WAMMI_", names(weights$LMM)[i])]] <-
      wammi_model(data_contaminated, weight = weights$LMM[[i]], Ncomp = Ncomp)

    if (i == 5) {
      RAMMI[[paste0("RAMMI_", names(weights$RLMM)[i])]] <-
        rammi_model(data_contaminated, weight = NULL, Ncomp = Ncomp)
    } else {
      RAMMI[[paste0("RAMMI_", names(weights$RLMM)[i])]] <-
        rammi_model(data_contaminated, weight = weights$RLMM[[i]], Ncomp = Ncomp)
    }

    RWAMMI[[paste0("RWAMMI_", names(weights$RLMM)[i])]] <-
      rwammi_model(data_contaminated, weight = weights$RLMM[[i]], Ncomp = Ncomp, rSVD = rSVD)
  }

  ALLModels <- list(RealSVD = Estimated.real,
                    AMMI = AMMI.model,
                    WAMMI = WAMMI,
                    RAMMI = RAMMI,
                    RWAMMI = RWAMMI)
  return(ALLModels)
}
