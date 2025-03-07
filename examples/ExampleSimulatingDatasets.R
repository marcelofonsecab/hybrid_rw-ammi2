# Please read the README file before starting this code
# In particular you have indications on how to install the hybridRWAMMI package

# Example of Simulating Synthetic Data Using the hybridRWAMMI Package
# This script demonstrates how to:
#   - Generate synthetic genotype-environment data using the sim_amb function.
#   - Create replications of the simulated data.
#   - Introduce various types of contamination into the data (shift, variance inflation, point mass).
#   - Compute error variances and model weights based on the contaminated data.
#
# =============================================================================
# REQUIRED PACKAGES (if not already installed)
# install.packages("dplyr")
# install.packages("pbapply")
# install.packages("lme4")
# install.packages("robustlmm")
# install.packages("MASS")
# install.packages("reshape2")
# install.packages("rsvddpd")
# =============================================================================
#
library(hybridRWAMMI)

# Set contamination parameters:
shift_values <- c(4, 7, 10)          # Shift values (k) for 'shift' contamination
varinflated_values <- 0.2            # Variance inflation values (c) for 'varinflated' contamination
pointmass_values <- 10               # Point mass values (c) for 'pointmass' contamination
percent_values <- 0.02               # Contamination percentage as 2% (you can use a vector like c(0.02, 0.05, 0.10) for multiple percents of contamination)
seeds.values <- 1                    # Seed values; here using seed 1 (you can use a vector like 1:100 for multiple simulations)

# Optionally, if you want to speed up simulations using parallel processing,
# you can create a cluster using the 'parallel' package. For this example, we use no cluster.
cl <- NULL

# Initialize an empty list to store simulation results
results_list <- list()

for(i in 1:length(seeds.values)) {
  set.seed(seeds.values[i])
  cat("Seed:", seeds.values[i], "\n")

  # Generate synthetic data with the sim_amb function
  sim_data <- sim_amb(seed = seeds.values[i],
                      Ngen = 100, Nenv = 8, Ncomp = 2,
                      effectGlobal = c(mean = 15, sd = sqrt(3)),
                      effectGen = c(mean = 5, sd = 1),
                      effectEnv = c(mean = 8, sd = sqrt(2)),
                      k = c(28, 15))

  # Create replications of the data (e.g., 2 replications with residual standard deviation of 1)
  rep_data <- Create_Replications(sim_data, reps = 2, sig = 1)

  for(j in 1:length(percent_values)) {
    cat("Percentage:", percent_values[j], "\n")

    # ---------------------------
    # Introduce 'shift' type contamination
    for(k in 1:length(shift_values)) {
      cat("Shift:", shift_values[k], "\n")
      data_cont <- Data_Contamination(data = rep_data,
                                      percentage = percent_values[j],
                                      seed = seeds.values[i],
                                      type = "shift",
                                      k = shift_values[k])
      err_var <- suppressMessages(Error_Var(data_cont$data.cont, cluster = cl))
      weights_tmp <- R_Weights(data_cont$data.cont, err_var)

      results_list[[paste0("shift_perc", percent_values[j]*100, "_k", shift_values[k])]] <- list(
        data = data_cont,
        errorvariances = err_var,
        weights = weights_tmp
      )
    }

    # ---------------------------
    # Introduce 'varinflated' type contamination
    for(k in 1:length(varinflated_values)) {
      cat("High Variance (varinflated):", varinflated_values[k], "\n")
      data_cont <- Data_Contamination(data = rep_data,
                                      percentage = percent_values[j],
                                      seed = seeds.values[i],
                                      type = "varinflated",
                                      c = varinflated_values[k])
      err_var <- suppressMessages(Error_Var(data_cont$data.cont, cluster = cl))
      weights_tmp <- R_Weights(data_cont$data.cont, err_var)

      results_list[[paste0("varinflated_perc", percent_values[j]*100)]] <- list(
        data = data_cont,
        errorvariances = err_var,
        weights = weights_tmp
      )
    }

    # ---------------------------
    # Introduce 'pointmass' type contamination
    for(k in 1:length(pointmass_values)) {
      cat("PointMass:", pointmass_values[k], "\n")
      data_cont <- Data_Contamination(data = rep_data,
                                      percentage = percent_values[j],
                                      seed = seeds.values[i],
                                      type = "pointmass",
                                      c = pointmass_values[k])
      err_var <- suppressMessages(Error_Var(data_cont$data.cont, cluster = cl))
      weights_tmp <- R_Weights(data_cont$data.cont, err_var)

      results_list[[paste0("pointmass_perc", percent_values[j]*100)]] <- list(
        data = data_cont,
        errorvariances = err_var,
        weights = weights_tmp
      )
    }
  }

  cat("No contamination\n")
  err_var <- suppressMessages(Error_Var(rep_data, cluster = cl))
  weights_tmp <- R_Weights(rep_data, err_var)

  results_list[["nocontamination"]] <- list(
    data = rep_data,
    errorvariances = err_var,
    weights = weights_tmp
  )
}

# Example outputs:
results_list$shift_perc2_k4      # Data with shift contamination using k = 4
results_list$shift_perc2_k7      # Data with shift contamination using k = 7
results_list$shift_perc2_k10     # Data with shift contamination using k = 10
results_list$varinflated_perc2   # Data with variance-inflated contamination
results_list$pointmass_perc2     # Data with point mass contamination
results_list$nocontamination     # Data without contamination
