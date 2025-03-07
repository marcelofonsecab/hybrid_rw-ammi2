# hybridRWAMMI

**hybridRWAMMI** is an R package developed for genotype-by-environment interaction (GEI) analysis in crop improvement breeding programs. It implements the classic AMMI model along with robust and weighted variations, and introduces a hybrid framework called **RW-AMMI** that integrates robust and weighted algorithms. This approach enhances the model’s ability to handle atypical data points—arising from measurement errors, genotype characteristics, diseases, or climatic extremes—thus maintaining reliable inferential results in the presence of data contamination or heterogeneous error variances.

This R package provides the implementation of the methods described in the paper *A Robust-Weighted AMMI Modeling Approach with Generalized Weighting Schemes*. It includes the code used for the analysis and modeling techniques presented in the paper, ensuring reproducibility and facilitating application to similar datasets. In addition, the package can also be used for the analysis of other two-way datasets, as illustrated in the examples described in the paper.

## Overview

The AMMI model and its variations have long been used to identify genotypes with specific adaptability and stability under diverse environmental conditions. However, contamination from measurement errors, diseases, or extreme climatic events may violate the assumptions of the model, leading to unreliable conclusions. **hybridRWAMMI** addresses these challenges by:

- **Robustifying the AMMI Model:** Incorporating robust algorithms (R-AMMI) to lessen the impact of outliers.
- **Weighted Approaches:** Introducing nine weighting schemes for weighted (W-AMMI), robust (R-AMMI), and hybrid (RW-AMMI) models.
- **Hybrid Framework (RW-AMMI):** Combining robust and weighted methodologies to enhance model performance, especially under contaminated conditions.
- **Comprehensive Evaluation:** Utilizing Monte Carlo simulations under both contaminated and uncontaminated scenarios, including heterogeneous error variance conditions, and validating the method with real crop data using ensemble strategies for improved genotype recommendation.

## Features

- **Data Transformation & Simulation:**  
  Easily transform GEI data and simulate synthetic datasets using functions like `sim_amb` and `Create_Replications`.

- **Model Fitting:**  
  Fit various models including the conventional AMMI, weighted AMMI (W-AMMI), robust AMMI (R-AMMI), and the hybrid robust weighted AMMI (RW-AMMI).

- **Error Variance & Weight Calculation:**  
  Compute error variances using linear mixed models (LMM) and robust LMM (RLMM), and generate weights with `Error_Var` and `R_Weights` using nine different weighting schemes.

- **SVD and Robust SVD:**  
  Perform standard and robust weighted singular value decomposition with functions such as `WeightedSvdRes`.

- **Model Evaluation:**  
  Calculate key GEI metrics (Mean Percentage of Explained Variability [MPEV], Mean Squared Error on singular values, Mean trimmed Squared Prediction Error [MtSPE], and Maxsub) using `GEI_Metrics` and `Metricsapply`.

- **Visualization:**  
  Generate informative biplots with `BiplotCreation` to visualize interaction patterns and model performance.

## Installation

You can install **hybridRWAMMI** from GitHub

```r
# install.packages("devtools")
devtools::install_github("marcelofonsecab/hybridRWAMMI")
```

Ensure you have the following required packages installed:

- dplyr
- pbapply
- lme4
- robustlmm
- MASS
- reshape2
- rsvddpd

## Usage

### Analyzing Real Data

Your genotype-environment dataset must be a data frame with exactly four columns in this order:

1. **gen**  - Factor indicating genotype identifiers.
2. **env**  - Factor representing the environments where measurements were taken.
3. **rep**  - Factor indicating the replication number (distinguishing repeated measurements).
4. **yield** - Numeric vector of observed trait values.

The analysis of the real Steptoe x Morex data presented in **Section 5** of the manuscript is available in the R script [`SteptoexMorexDataAnalysis.R`](examples/SteptoexMorexDataAnalysis.R). This script is designed to be self-explanatory and provides step-by-step instructions on how to:
- How to load your dataset (e.g., the *SteptoexMorex* dataset).
- How to compute error variances and weights.
- Fitting various AMMI-based models using SVD.
- Running diagnostic tests (normality, homogeneity of variances).
- Selecting the best genotype per environment.
- Generating biplots for visualization.

For any other dataset in the same format, the user can simply replace the Steptoe x Morex dataset with their new dataset to apply the same analysis.

### Simulating Data and Contamination

The package includes functions to generate synthetic GEI data and simulate various contamination scenarios presented in **Section 3** of the presented paper. See the examples:

- [`ExampleSimulatingDatasets.R`](examples/ExampleSimulatingDatasets.R)  
  *Demonstrates how to generate synthetic data, create replications, and apply different contamination types (shift, variance inflation, point mass).*

- [`SimulatedDataAnalysis.R`](examples/SimulatedDataAnalysis.R)  
  *Shows how to load simulated data, run SVD-based models, compute GEI metrics, and generate biplots for both contaminated and uncontaminated datasets.*

### Model Evaluation & Biplot Visualization

After fitting the models, you can evaluate them using `GEI_Metrics` and `Metricsapply`, which calculate metrics such as MPEV, MSE on singular values, MtSPE, and Maxsub. Biplots generated by `BiplotCreation` provide visual insights into the interaction patterns, supporting decisions in crop improvement.

## Documentation

All functions are documented with Roxygen comments. You can view documentation in R by running, for example:

```r
?ammi_model
?WeightedSvdRes
?GEI_Metrics
```
