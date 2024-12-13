#' Power K-means Clustering with Bregman Divergence
#'
#' @description
#' Implements power k-means clustering algorithm with Bregman divergence for robust clustering,
#' supporting multiple data types and dimension reduction methods. The algorithm is particularly
#' effective for clustering exponential family data.
#'
#' @importFrom irlba irlba
#' @importFrom matrixStats colMins colSds  
#' @importFrom stats rnorm
#' 
#' @param X Input data: matrix, sparse matrix (dgCMatrix), SingleCellExperiment, or Seurat object
#' @param s Initial power parameter for weights calculation. Must be negative (e.g., -0.5)
#' @param k Number of clusters
#' @param npcs Number of principal components for dimension reduction. Default: min(50, ncol(X))
#' @param eta Learning rate for power parameter updates. Default: 1.05
#' @param scale Logical indicating whether to scale the data. Default: TRUE
#' @param initial_centers Initial cluster centers matrix or NULL for random initialization
#' @param divergence Type of Bregman divergence: "euclidean", "kl", "itakura-saito", or "logistic"
#' @param convergence_threshold Number of iterations with unchanged assignments for convergence. Default: 5
#' @param seed Random seed for reproducibility. Default: 1
#' @param max_iter Maximum number of iterations. Default: 1000
#' @param dim_reduction Dimension reduction method: "pca", "spectral", or "none". Default: "pca"
#' @param bw_spectral Bandwidth parameter for spectral clustering. Default: 1
#'
#' @return List containing:
#' \itemize{
#'   \item centers: Matrix of final cluster centers
#'   \item classes: Vector of cluster assignments
#'   \item reduced_dim: Dimension-reduced data matrix
#' }
#'
#' @details
#' The algorithm implements an adaptive power k-means approach with Bregman divergence,
#' specifically designed for clustering exponential family data. For KL and Itakura-Saito 
#' divergences, data is automatically shifted to ensure positivity. The power parameter s 
#' is updated during optimization to improve clustering robustness.
#'
#' @references
#' Vellal, A., Eldan, R., Zaslavsky, M., & Roeder, K. (2022). 
#' Bregman Power k-Means for Clustering Exponential Family Data. 
#' In International Conference on Machine Learning (pp. 22083-22101). PMLR.
#' https://proceedings.mlr.press/v162/vellal22a/vellal22a.pdf
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(1000), ncol=10)
#' result <- power_kmeans_bregman(X, s= -0.5, k=3)
#'
#' @export
power_kmeans_bregman <- function(
    X,
    s,
    k,
    npcs = min(50, ncol(X)),
    eta = 1.05,
    scale = TRUE,
    initial_centers = NULL,
    divergence = "euclidean",
    convergence_threshold = 5,
    seed = 1,
    max_iter = 1000,
    dim_reduction = "pca",
    bw_spectral = 1) {
  
  # Input validation
  if (k < 1 || k > nrow(X)) {
    stop("Invalid number of clusters")
  }
  if (!dim_reduction %in% c("pca", "spectral", "none")) {
    stop("Invalid dimension reduction method")
  }
  if (nrow(X) < 2 || ncol(X) < 1) {
    stop("Input matrix must have at least 2 rows and 1 column")
  }
  if (!is.numeric(s) || length(s) != 1 || is.na(s)) {
    stop("'s' must be a single numeric value")
  }
  if (!is.numeric(k) || k < 1 || k > nrow(X) || is.na(k)) {
    stop("'k' must be a positive integer less than nrow(X)")
  }
  if (!is.numeric(npcs) || npcs < 1 || is.na(npcs)) {
    stop("'npcs' must be a positive integer")
  }
  if (!is.numeric(s) || s >= 0) {
    stop("Parameter 's' must be negative")
  }
  
  # Data type using switch: matrix, sparse matrix, SingleCellExperiment, and Seurat
  data_type <- class(X)[1]
  do_scale <- FALSE

  X <- switch(data_type,
    "matrix" = {
      do_scale <- scale
      if (scale) {
        print("Data will be scaled later. Set scale = FALSE to use unscaled data.")
      } else {
        print("Data will be used as is.")
      }
      X
    },
    "dgCMatrix" = {
      do_scale <- scale
      if (scale) {
        print("Data will be scaled later. Set scale = FALSE to use unscaled data.")
      } else {
        print("Data will be used as is.")
      }
      as.matrix(X)
    },
    "SingleCellExperiment" = {
      do_scale <- scale
      if (scale) {
        print("Data will be scaled later. Set scale = FALSE to use unscaled data.")
      } else {
        print("Data will be used as is.")
      }
      t(SummarizedExperiment::assay(X))
    },
    "Seurat" = {
      layer <- if ("scale.data" %in% SeuratObject::Layers(X) && scale) {
        print("Scaled data exists, using scaled data")
        "scale.data"
      } else if ("data" %in% SeuratObject::Layers(X)) {
        do_scale <- scale
        if (scale) {
          print("Data is not scaled and will be scaled later. Set scale = FALSE to use unscaled data.")
        } else {
          print("Unscaled data will be centered and used")
        }
        "data"
      } else {
        stop("Neither data nor scale.data found in Seurat object")
      }
      t(Seurat::GetAssayData(X, layer = layer))
    },
    stop("Unsupported data type for X")
  )
  
  if (dim_reduction == "pca") {
    # Scale the data, remove features with zero variances
    sd_vec <- matrixStats::colSds(as.matrix(X)) # vectorized
    X <- X[, sd_vec > 0, drop = FALSE]
    X <- scale(X, center = TRUE, scale = do_scale)
    X_svd <- if (ncol(X) < 500) { # Dimension reduction for X
      svd(t(X), nu = npcs, nv = npcs)
    } else {
      irlba::irlba(t(X), nu = npcs, nv = npcs)
    }
    X_dim_reduced <- X_svd$v[, 1:npcs, drop = FALSE] %*% diag(X_svd$d[1:npcs])
  } 
  else if (dim_reduction == "spectral") {
    # Scale the data, remove features with zero variances
    sd_vec <- matrixStats::colSds(as.matrix(X)) # vectorized
    X <- X[, sd_vec > 0, drop = FALSE]
    X <- scale(X, center = TRUE, scale = do_scale)
    dist_mt <- as.matrix(dist(X, diag = TRUE, upper = TRUE))
    W <- exp(-dist_mt^2 / (2 * bw_spectral^2))
    D <- diag(rowSums(W))
    L <- D - W
    X_dim_reduced <- eigen(L, symmetric = TRUE, only.values = FALSE)$vectors[, nrow(L):(nrow(L) - npcs + 1)]
  } 
  else if (dim_reduction == "none") {
    print("no dimension reduction applied")
    X_dim_reduced <- X
  }
  # Validate divergence type
  if (!divergence %in% c("euclidean", "kl", "itakura-saito", "logistic")) {
    stop("Invalid divergence type")
  }

  # For KL and Itakura-Saito, ensure data is positive
  if (divergence %in% c("kl", "itakura-saito")) {
    shift_vec <- pmax(-matrixStats::colMins(X_dim_reduced), 0) + 1e-6
    X_dim_reduced <- sweep(X_dim_reduced, 2, shift_vec, "+")
    if (any(X_dim_reduced <= 0)) {
      stop("Data must be positive for KL or Itakura-Saito divergence")
    }
  }

  # Initialize centers if not provided
  if (is.null(initial_centers)) {
    set.seed(seed)
    centers <- X_dim_reduced[sample(nrow(X_dim_reduced), k), , drop = FALSE] +
      matrix(rnorm(k * ncol(X_dim_reduced), sd = 0.1), nrow = k)
  } 
  else {
    centers <- initial_centers
  }

  classes <- integer(nrow(X_dim_reduced))
  classes_old <- integer(nrow(X_dim_reduced))
  convergence_cnt <- 0

  # Main optimization loop with vectorized operations
  for (iter in 1:max_iter) {
    # Pairwise distances efficiently
    pairwise_dist <- pairwise_bregman(X_dim_reduced, centers)
    weights <- (pairwise_dist^(s - 1)) * (rowSums(pairwise_dist^s)^(1 / s - 1))

    # Update centers
    centers <- t(vapply(1:k, function(i) {
      colSums(X_dim_reduced * weights[, i]) / sum(weights[, i])
    }, FUN.VALUE = numeric(ncol(X_dim_reduced))))

    # Update power parameter s
    if (iter %% 2 == 0) {
      if (s > -1.0) {
        s <- s - 0.2
      } 
      else if (s > -120.0) {
        s <- s * eta
      }
    }

    # Update class assignments
    if (iter > 1) classes_old <- classes
    classes <- max.col(-pairwise_dist) # Faster than apply + which.min

    # Check convergence
    if (iter > 1 && identical(classes, classes_old)) {
      convergence_cnt <- convergence_cnt + 1
      if (convergence_cnt >= convergence_threshold) 
        break
    } 
    else {
      convergence_cnt <- 0
    }
  }

  # Create and return PowerKmeans object
  create_powerkmeans(
    centers = centers,
    classes = classes,
    reduced_dim = X_dim_reduced,
    params = list(
      s = s,
      k = k,
      eta = eta,
      divergence = divergence,
      dim_reduction = dim_reduction,
      convergence_threshold = convergence_threshold,
      max_iter = max_iter
    )
  )
}