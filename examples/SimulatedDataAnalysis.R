# ============================================================
# Example Analysis of Simulated Data Using the hybridRWAMMI Package
#
# In this example, we load a simulated dataset with 2% contamination and a shift type
# contamination (shift k = 4). The package provides several simulated datasets:
#
#   SDATA_notcontaminated   - Not contaminated data
#   SDATA_p2_shift4         - 2% contamination with shift (k = 4)
#   SDATA_p5_shift4         - 5% contamination with shift (k = 4)
#   SDATA_p10_shift4        - 10% contamination with shift (k = 4)
#   SDATA_p2_shift7         - 2% contamination with shift (k = 7)
#   SDATA_p5_shift7         - 5% contamination with shift (k = 7)
#   SDATA_p10_shift7        - 10% contamination with shift (k = 7)
#   SDATA_p2_shift10        - 2% contamination with shift (k = 10)
#   SDATA_p5_shift10        - 5% contamination with shift (k = 10)
#   SDATA_p10_shift10       - 10% contamination with shift (k = 10)
#   SDATA_p2_pmass10        - 2% contamination with point mass (c = 10)
#   SDATA_p5_pmass10        - 5% contamination with point mass (c = 10)
#   SDATA_p10_pmass10       - 10% contamination with point mass (c = 10)
#   SDATA_p2_varinfl_2      - 2% contamination with variance inflation (c = 0.2)
#   SDATA_p5_varinfl_2      - 5% contamination with variance inflation (c = 0.2)
#   SDATA_p10_varinfl_2     - 10% contamination with variance inflation (c = 0.2)
#
# ============================================================

# Load required packages and the hybridRWAMMI package
library(hybridRWAMMI)
# (Make sure that the required packages are installed as noted in the package documentation)

# Load the simulated data; here we use the dataset with 2% contamination and shift (k = 4)
alldata <- SDATA_p2_shift4

# ----------------------------------------------------------------
# Step 1: Process the first dataset in alldata
# ----------------------------------------------------------------

attach(alldata[[1]])

# Compute SVD-based model fits using both the non-contaminated and contaminated data
allsvds <- All_SVDS(data$data.not_contaminated,
                    data$data.contaminated,
                    weights = weights)

# Calculate GEI metrics for the obtained SVD results
.allMetrics <- GEI_Metrics(allsvds)

# ----------------------------------------------------------------
# Step 2: Process the remaining datasets in alldata
# ----------------------------------------------------------------

# For example, iterate over datasets 2 to 5 in the list (adjust as needed)
for (i in 2:5) {
  attach(alldata[[i]])  # Attach the current dataset

  # Compute SVD-based models for the current dataset
  allsvds <- All_SVDS(data$data.not_contaminated,
                      data$data.contaminated,
                      weights = weights)

  # Calculate GEI metrics for the current dataset
  .tmpMetrics <- GEI_Metrics(allsvds)

  # Append the current metrics to the overall metrics list
  .allMetrics <- appendList(.allMetrics, .tmpMetrics)
}

# ----------------------------------------------------------------
# Step 3: Calculate and Combine Metrics for All Models
# ----------------------------------------------------------------

# Define the model names for the results (adjust or extend as needed)
model_names <- c("AMMI",
                 "W-AMMI{1}", "W-AMMI{2}", "W-AMMI{3}",
                 "W-AMMI{4}", "W-AMMI{5}", "W-AMMI{6}",
                 "W-AMMI{7}", "W-AMMI{8}", "W-AMMI{9}",
                 "R-AMMI{1}", "R-AMMI{2}", "R-AMMI{3}",
                 "R-AMMI{4}", "R-AMMI{5}", "R-AMMI{6}",
                 "R-AMMI{7}", "R-AMMI{8}", "R-AMMI{9}",
                 "RW-AMMI{1}", "RW-AMMI{2}", "RW-AMMI{3}",
                 "RW-AMMI{4}", "RW-AMMI{5}", "RW-AMMI{6}",
                 "RW-AMMI{7}", "RW-AMMI{8}", "RW-AMMI{9}")

# Calculate the "MPEV" metric (Mean Percentage of Explained Variability)
# with hypothesis testing (null hypothesis h0 = 100, alternative "greater")
Values <- Metricsapply(.allMetrics, "MPEV",
                       hypothesis.testing = TRUE,
                       h0 = 100, side = "greater")

# Assign the model names to the metrics and p-values
names(Values$Metric) <- names(Values$pvalues) <- model_names

# Note: To create the tables shown in the paper, repeat this procedure for each metric
# and for every contamination scenario (including the non-contaminated case),
# then bind all values into one matrix for presentation.
# Just change the "MPEV" to: "MSE.Singvals1"; "MSE.Singvals2"; "MtSPE"; or "Maxsub"

# ----------------------------------------------------------------
# Step 4: Generate Biplots
# ----------------------------------------------------------------

# For the biplots, we select a random seed to choose one simulated dataset.
# Here, we set the seed and sample one index from the available simulated datasets.
set.seed(768622)  # This seed was generated from a uniform distribution (range -999999 to 999999)
dataseed <- sample(1:100, 1)

# For example, use a dataset with 5% contamination and shift (k = 7)
alldata <- SDATA_p5_shift7
attach(alldata[[dataseed]])

# For Figure S5 (biplots without contamination), perform SVD on the non-contaminated data:
allsvds.notcontaminated <- All_SVDS(data$data.not_contaminated,
                                    data$data.not_contaminated,
                                    weights = weights)

# --- Generate biplots for non-contaminated data ---
## AMMI biplot
BiplotCreation(allsvds.notcontaminated$AMMI$SVD)

## W-AMMI biplots (e.g., model 1, model 5, model 6)
BiplotCreation(allsvds.notcontaminated$WAMMI$WAMMI_lmm.Weights.Env$SVD)  # W-AMMI{1}
BiplotCreation(allsvds.notcontaminated$WAMMI$WAMMI_Weights.Rlm$SVD)       # W-AMMI{5}
BiplotCreation(allsvds.notcontaminated$WAMMI$WAMMI_lmm.Weights.RlmdotEnv$SVD)  # W-AMMI{6}

## R-AMMI biplots (e.g., model 1, model 5, model 6)
BiplotCreation(allsvds.notcontaminated$RAMMI$RAMMI_rlmm.Weights.Env$SVD)    # R-AMMI{1}
BiplotCreation(allsvds.notcontaminated$RAMMI$RAMMI_Weights.Rlm$SVD)         # R-AMMI{5}
BiplotCreation(allsvds.notcontaminated$RAMMI$RAMMI_rlmm.Weights.RlmdotEnv$SVD)    # R-AMMI{6}

## RW-AMMI biplots (e.g., model 1, model 5, model 6)
BiplotCreation(allsvds.notcontaminated$RWAMMI$RWAMMI_rlmm.Weights.Env$SVD)  # RW-AMMI{1}
BiplotCreation(allsvds.notcontaminated$RWAMMI$RWAMMI_Weights.Rlm$SVD)       # RW-AMMI{5}
BiplotCreation(allsvds.notcontaminated$RWAMMI$RWAMMI_rlmm.Weights.RlmdotEnv$SVD)  # RW-AMMI{6}

# --- Generate biplots for contaminated data (Figures S6-S8) ---
# Now perform SVD on the contaminated data from the selected seed (e.g., 5% contamination with shift k = 7)
allsvds <- All_SVDS(data$data.not_contaminated,
                    data$data.contaminated,
                    weights = weights)

## AMMI biplots:
BiplotCreation(allsvds$RealSVD$SVD)  # Non-contaminated reconstruction
BiplotCreation(allsvds$AMMI$SVD)       # Contaminated AMMI model

## W-AMMI biplots:
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.Env$SVD)       # W-AMMI{1}
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.Gen$SVD)         # W-AMMI{2}
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.GendotEnv$SVD)     # W-AMMI{3}
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.GenplusEnv$SVD)    # W-AMMI{4}
BiplotCreation(allsvds$WAMMI$WAMMI_Weights.Rlm$SVD)             # W-AMMI{5}
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.RlmdotEnv$SVD)     # W-AMMI{6}
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.RlmdotGen$SVD)     # W-AMMI{7}
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.RlmdotGendotEnv$SVD) # W-AMMI{8}
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.RlmdotGenplusEnv$SVD)# W-AMMI{9}

## R-AMMI biplots:
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.Env$SVD)         # R-AMMI{1}
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.Gen$SVD)           # R-AMMI{2}
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.GendotEnv$SVD)       # R-AMMI{3}
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.GenplusEnv$SVD)      # R-AMMI{4}
BiplotCreation(allsvds$RAMMI$RAMMI_Weights.Rlm$SVD)                # R-AMMI{5}
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.RlmdotEnv$SVD)       # R-AMMI{6}
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.RlmdotGen$SVD)       # R-AMMI{7}
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.RlmdotGendotEnv$SVD)   # R-AMMI{8}
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.RlmdotGenplusEnv$SVD)  # R-AMMI{9}

## RW-AMMI biplots:
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.Env$SVD)         # RW-AMMI{1}
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.Gen$SVD)           # RW-AMMI{2}
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.GendotEnv$SVD)       # RW-AMMI{3}
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.GenplusEnv$SVD)      # RW-AMMI{4}
BiplotCreation(allsvds$RWAMMI$RWAMMI_Weights.Rlm$SVD)                # RW-AMMI{5}
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.RlmdotEnv$SVD)       # RW-AMMI{6}
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.RlmdotGen$SVD)       # RW-AMMI{7}
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.RlmdotGendotEnv$SVD)   # RW-AMMI{8}
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.RlmdotGenplusEnv$SVD)  # RW-AMMI{9}
