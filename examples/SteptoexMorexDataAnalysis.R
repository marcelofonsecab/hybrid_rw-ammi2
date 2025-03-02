# Example Application of the MyHybridAMMI Package
# This script demonstrates a step-by-step analysis using real genotype-environment data.
# All comments and variable names are in English to ensure global usability.

# Load the MyHybridAMMI package
library(hybrid_rw-ammi)

# Load your dataset.
# In this example, we are using the SteptoexMorex dataset.
# For more information about the dataset, you can use: ?SxM
data <- SxM

# IMPORTANT: Your dataset must be a data frame with exactly four columns in the following order:
#   1. gen  - A factor indicating the genotype identifiers.
#   2. env  - A factor representing the environments where the measurements were taken.
#   3. rep  - A factor indicating the replication number for each observation.
#           This distinguishes repeated measurements for the same genotype-environment combination.
#   4. yield - A numeric vector containing the observed trait values.
# It is essential that the column names are exactly: gen, env, rep, yield,
# otherwise the package functions may not work correctly.


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

# Step 1: Calculate the error variance for the dataset
errorVar <- Error_Var(data)

# Step 2: Calculate the weights for the dataset using the computed error variances
weights <- R_Weights(data, errorVar)

# Step 3: Compute all SVD-based model fits.
#         Here, we use the same dataset for both non-contaminated and contaminated data,
#         and set rSVD = TRUE to apply robust weighted SVD iteratively.
allSVDs <- All_SVDS(data, data, weights, rSVD = TRUE)

# Step 4: Perform diagnostic tests
# 4a: Test the normality of residuals from the AMMI model
shapiro.test(allSVDs$AMMI$residual)

# 4b: Test homogeneity of variances in yield grouped by genotype
fligner.test(yield ~ gen, data)

# 4c: Test homogeneity of variances in yield grouped by environment
fligner.test(yield ~ env, data)

# Step 5: Identify the best genotype for each environment
# Change bestRank to 2 or 3 to obtain the second or third best genotype, respectively.
bestRank <- 1
genotypeNames <- as.character(unique(data$gen))
environmentNames <- as.character(unique(data$env))

# Combine the best genotype selections from various model outputs
bestGenotypesMatrix <- cbind(
  # AMMI model
  best_envgen(allSVDs$AMMI, quant = bestRank, genotypeNames, environmentNames),

  # W-AMMI models
  best_envgen(allSVDs$WAMMI$WAMMI_lmm.Weights.Env, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$WAMMI$WAMMI_lmm.Weights.Gen, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$WAMMI$WAMMI_lmm.Weights.GendotEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$WAMMI$WAMMI_lmm.Weights.GenplusEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$WAMMI$WAMMI_Weights.Rlm, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$WAMMI$WAMMI_lmm.Weights.RlmdotEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$WAMMI$WAMMI_lmm.Weights.RlmdotGen, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$WAMMI$WAMMI_lmm.Weights.RlmdotEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$WAMMI$WAMMI_lmm.Weights.RlmdotGenplusEnv, quant = bestRank, genotypeNames, environmentNames),

  # R-AMMI models
  best_envgen(allSVDs$RAMMI$RAMMI_rlmm.Weights.Env, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RAMMI$RAMMI_rlmm.Weights.Gen, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RAMMI$RAMMI_rlmm.Weights.GendotEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RAMMI$RAMMI_rlmm.Weights.GenplusEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RAMMI$RAMMI_Weights.Rlm, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RAMMI$RAMMI_rlmm.Weights.RlmdotEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RAMMI$RAMMI_rlmm.Weights.RlmdotGen, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RAMMI$RAMMI_rlmm.Weights.RlmdotEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RAMMI$RAMMI_rlmm.Weights.RlmdotGenplusEnv, quant = bestRank, genotypeNames, environmentNames),

  # RW-AMMI models
  best_envgen(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.Env, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.Gen, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.GendotEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.GenplusEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RWAMMI$RWAMMI_Weights.Rlm, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.RlmdotEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.RlmdotGen, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.RlmdotEnv, quant = bestRank, genotypeNames, environmentNames),
  best_envgen(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.RlmdotGenplusEnv, quant = bestRank, genotypeNames, environmentNames)
)

# Display the final table of best genotypes for each environment
t(bestGenotypesMatrix)

# Step 6: Create Biplots to Visualize the Results

# Set environment names based on the AMMI estimated data
envNames <- colnames(allSVDs$AMMI$estimated.data)

## 6a: Biplot for the AMMI model
BiplotCreation(allSVDs$AMMI$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))

## 6b: Biplots for selected W-AMMI models ({1}, {5}, {6})
BiplotCreation(allSVDs$WAMMI$WAMMI_lmm.Weights.Env$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$WAMMI$WAMMI_Weights.Rlm$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$WAMMI$WAMMI_lmm.Weights.RlmdotEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))

## 6c: Biplots for selected R-AMMI models ({1}, {5}, {6})
BiplotCreation(allSVDs$RAMMI$RAMMI_rlmm.Weights.Env$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RAMMI$RAMMI_Weights.Rlm$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RAMMI$RAMMI_rlmm.Weights.RlmdotEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))

## 6d: Biplots for selected RW-AMMI models ({1}, {5}, {6})
BiplotCreation(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.Env$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RWAMMI$RWAMMI_Weights.Rlm$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.RlmdotEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))

# Step 7: (Supplementary) Create additional biplots for further model comparisons

## Biplots for additional W-AMMI models ({2}, {3}, {4}, {7}, {8}, {9})
BiplotCreation(allSVDs$WAMMI$WAMMI_lmm.Weights.Gen$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$WAMMI$WAMMI_lmm.Weights.GendotEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$WAMMI$WAMMI_lmm.Weights.GenplusEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$WAMMI$WAMMI_lmm.Weights.RlmdotGen$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$WAMMI$WAMMI_lmm.Weights.RlmdotGendotEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$WAMMI$WAMMI_lmm.Weights.GenplusEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))

## Biplots for additional R-AMMI models ({2}, {3}, {4}, {7}, {8}, {9})
BiplotCreation(allSVDs$RAMMI$RAMMI_rlmm.Weights.Gen$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RAMMI$RAMMI_rlmm.Weights.GendotEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RAMMI$RAMMI_rlmm.Weights.GenplusEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RAMMI$RAMMI_rlmm.Weights.RlmdotGen$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RAMMI$RAMMI_rlmm.Weights.RlmdotGendotEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RAMMI$RAMMI_rlmm.Weights.GenplusEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))

## Biplots for additional RW-AMMI models ({2}, {3}, {4}, {7}, {8}, {9})
BiplotCreation(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.Gen$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.GendotEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.GenplusEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.RlmdotGen$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.RlmdotGendotEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
BiplotCreation(allSVDs$RWAMMI$RWAMMI_rlmm.Weights.GenplusEnv$SVD,
               env.names = envNames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7),
               lims.V = c(-2.94, 3.2))
