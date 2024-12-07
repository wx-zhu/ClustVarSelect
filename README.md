# ClustVarSelect

ClustVarSelect implements feature selection for high-dimensional data using sparse manifold decomposition (SMD) with power k-means clustering. The package excels at analyzing single-cell RNA sequencing data through:

- Power k-means clustering with multiple Bregman divergences (Euclidean, KL, Itakura-Saito, logistic)

- Entropy-based (decision trees) and maximum margin (SVM) classification

- Dimension reduction via PCA or spectral methods

- Support for sparse matrices, SingleCellExperiment, and Seurat objects

## Installation

```R
devtools::install_github("wx-zhu/ClustVarSelect")
```

## Key Features

- Robust Clustering: Power k-means with adaptive weighting for outlier resilience

- Multiple Divergences: Support for various Bregman divergences to handle different data distributions

- Flexible Classification: Choice of decision trees or L1-penalized SVM for feature scoring

- Dimension Reduction: Integrated PCA and spectral clustering options

- Bioinformatics Integration: Direct support for Seurat and SingleCellExperiment objects


## Basic Usage

```R
library(ClustVarSelect)

# Basic feature selection with default parameters
result <- SMD(X, k_guess = 3)

# With custom parameters
result <- SMD(X, k_guess = 3,
              prop_algo = "kmeans",
              divergence = "kl",
              dim_reduction = "spectral")

# Print summary of results
print(result)
summary(result)

# Using KL divergence with power k-means
result <- power_kmeans_bregman(X, 
                              s = 2.0,
                              k = 3,
                              divergence = "kl",
                              dim_reduction = "pca")

# Feature selection with custom parameters
result <- SMD(X,
              k_guess = 3,
              trials = 100,
              n_sub = 1000,
              prop_algo = "kmeans",
              class_algo = "maxmargin",
              z_score = TRUE,
              divergence = "euclidean",
              dim_reduction = "spectral")
```


### Authors

- Zhiyuan Yu (zyyu@umich.edu)

- Wenxuan Zhu (zhuwx@umich.edu)

- Adrian Con Garcia (acong@umich.edu)


### Reference

Vellal, A., et al. (2022). Bregman Power k-Means for Clustering Exponential Family Data. ICML.

Melton, S., & Ramanathan, S. (2021). Discovering a sparse set of pairwise discriminating features in high-dimensional data. Bioinformatics.

