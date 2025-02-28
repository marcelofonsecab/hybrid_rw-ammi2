#' Create Replications of Genotype-Environment Data
#'
#' This function generates replications of a genotype-environment dataset by adding random residuals.
#' It simulates experimental replications by perturbing the response variable (yield) with residual error,
#' allowing for the assessment of variability across replications.
#'
#' @param data A data frame containing the original genotype-environment data. The data frame must have the following columns:
#'   \describe{
#'     \item{gen}{A factor representing genotype identifiers.}
#'     \item{env}{A factor representing environment identifiers.}
#'     \item{yield}{A numeric vector representing the response variable.}
#'   }
#' @param reps An integer specifying the number of replications to create. Default is 2.
#' @param sig A numeric value specifying the standard deviation of the residuals to add. Default is 1.
#'
#' @return A data frame with the replicated data. Each replication is identified by the "rep" column.
#'
#' @examples
#' # Generate simulated data using sim.amb (assumes sim.amb is already defined and documented)
#' sim_data <- sim.amb(seed = 123)
#'
#' # Create 3 replications with a residual standard deviation of 1.5
#' rep_data <- Create.Replications(sim_data, reps = 3, sig = 1.5)
#' head(rep_data)
#'
#' @export
Create.Replications <- function(data, reps = 2, sig = 1) {

  # Determine the number of environments and genotypes from the data factors
  Nenv <- length(levels(data$env))
  Ngen <- length(levels(data$gen))
  N <- Nenv * Ngen

  # Initialize a container for the replicated data
  data.tmp <- NULL

  # Loop to create each replication
  for (i in 1:reps) {
    # Generate random residuals for the entire dataset
    resids <- rnorm(N, mean = 0, sd = sqrt(sig))

    # Create a copy of the original data and add residuals to the yield
    data.aux <- data
    data.aux$yield <- data.aux$yield + resids
    data.aux$rep <- as.factor(i)  # Mark the replication number as a factor

    # Append the replicated data to the overall dataset
    data.tmp <- rbind(data.tmp, data.aux)
  }

  return(data.tmp)
}
