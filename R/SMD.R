#' Sparse Manifold Decomposition (SMD) with Power K-means
#'
#' @description
#' Performs sparse manifold decomposition using power k-means clustering with Bregman
#' divergence for robust feature importance analysis in high-dimensional data, especially
#' for scRNA-seq data. Features are normalized before analysis and importance scores can
#' be standardized relative to shuffled data.
#'
#' @importFrom graphics hist
#' @importFrom stats median predict quantile rnorm sd
#' @importFrom utils head
#' @importFrom irlba irlba
#' @importFrom matrixStats colMins colSds
#'
#' @param X Numeric matrix with dimensions (number of data points) x (number of features)
#' @param k_guess Integer specifying best estimate of number of clusters present
#' @param trials Number of cluster proposals to average over in ensemble. Default: 2 * ncol(X)
#' @param n_sub Number of data points used in each proposal clustering. Default: 0.8 * nrow(X)
#' @param prop_algo Algorithm for constructing cluster proposals. Either "agglo" or "kmeans". If 'kmeans', uses power_kmeans_bregman. Default: "agglo"
#' @param class_algo Type of classifier used. Either "entropy" (decision trees) or "maxmargin" (L1-SVM). Default: "entropy"
#' @param cluster_prior_min Minimum number of clusters for proposals. Default: 3
#' @param z_score Logical indicating whether to return standardized scores. Default: TRUE
#' @param eta Power k-means learning rate parameter. Default: 1.05
#' @param initial_centers Matrix of initial cluster centers or NULL for random initialization
#' @param divergence Type of Bregman divergence: "euclidean", "kl", "itakura-saito", or "logistic"
#' @param convergence_threshold Number of iterations with unchanged assignments for convergence. Default: 5
#' @param seed Random seed for reproducibility. Default: 1
#' @param max_iter Maximum number of iterations for power k-means. Default: 1000
#' @param dim_reduction Dimension reduction method: "pca", "spectral", or "none". Default: "pca"
#'
#'
#' @details
#' The function implements sparse manifold decomposition (SMD) combining power k-means clustering
#' with Bregman divergence for robust clustering. This method is particularly effective for
#' discovering discriminating features in high-dimensional biological data. For entropy-based
#' classification, it uses decision trees (rpart). For maxmargin classification, it uses
#' L1-penalized SVM.
#'
#' @references
#' Melton, S., & Ramanathan, S. (2021). Discovering a sparse set of pairwise
#' discriminating features in high-dimensional data. Bioinformatics (Oxford, England),
#' 37(2), 202-212. https://doi.org/10.1093/bioinformatics/btaa690
#'
#' Vellal, A., Eldan, R., Zaslavsky, M., & Roeder, K. (2022).
#' Bregman Power k-Means for Clustering Exponential Family Data.
#' In International Conference on Machine Learning (pp. 22083-22101). PMLR.
#' https://proceedings.mlr.press/v162/vellal22a/vellal22a.pdf
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' X <- matrix(rnorm(1000), ncol = 20)
#'
#' # Basic usage
#' result <- SMD(X, k_guess = 3)
#'
#' # With custom parameters
#' result <- SMD(X,
#'   k_guess = 3,
#'   prop_algo = "kmeans",
#'   divergence = "kl",
#'   dim_reduction = "spectral"
#' )
#'
#' @export
SMD <- function(X,
                k_guess,
                trials = NULL,
                n_sub = NULL,
                prop_algo = "agglo",
                class_algo = "entropy",
                cluster_prior_min = NULL,
                z_score = TRUE,
                eta = 1.05, # for power_kmeans_bergman, adjusted if needed
                initial_centers = NULL,
                divergence = "euclidean",
                convergence_threshold = 5,
                seed = 1,
                max_iter = 1000,
                dim_reduction = "pca") {
  # Get dimensions
  N <- nrow(X)
  D <- ncol(X)

  # Input validation
  if (k_guess < 1) {
    stop("k_guess must be positive")
  }
  if (!is.null(trials) && trials < 1) {
    stop("Number of trials must be positive")
  }
  if (!is.null(n_sub) && (n_sub < 1 || n_sub > nrow(X))) {
    stop("Invalid number of samples for subsampling")
  }
  if (nrow(X) < 2 || ncol(X) < 1) {
    stop("Input matrix must have at least 2 rows and 1 column")
  }
  if (!is.numeric(k_guess) || length(k_guess) != 1 || k_guess < 1 || is.na(k_guess)) {
    stop("'k_guess' must be a positive integer")
  }
  if (!is.null(trials) && (!is.numeric(trials) || trials < 1 || is.na(trials))) {
    stop("'trials' must be NULL or a positive integer")
  }

  # Set default parameters
  if (is.null(n_sub)) n_sub <- as.integer(N * 0.8)
  if (is.null(trials)) trials <- as.integer(2 * D)
  if (is.null(cluster_prior_min)) cluster_prior_min <- 3

  # Normalize X
  X <- scale(X, center = FALSE, scale = apply(X, 2, sd))

  # Validate algorithms
  if (!prop_algo %in% c("agglo", "kmeans")) {
    stop(paste(prop_algo, "is not a supported clustering algorithm."))
  }
  if (!class_algo %in% c("entropy", "maxmargin")) {
    stop(paste(class_algo, "is not a supported classifier."))
  }

  # Set cluster range
  cluster_prior <- c(cluster_prior_min, as.integer(2 * k_guess))

  # Create power kmeans parameters list
  pkm_params <- list(
    eta = eta,
    initial_centers = initial_centers,
    divergence = divergence,
    convergence_threshold = convergence_threshold,
    seed = seed,
    max_iter = max_iter,
    dim_reduction = dim_reduction
  )

  # Choose classifier function based on class_algo
  classifier_func <- if (class_algo == "entropy") {
    function(X, n_sub, cluster_prior, prop_algo) {
      find_classifier_dims_entropy(X, n_sub, cluster_prior, prop_algo, pkm_params)
    }
  } else {
    function(X, n_sub, cluster_prior, prop_algo) {
      find_classifier_dims_maxmargin(X, n_sub, cluster_prior, prop_algo, pkm_params)
    }
  }

  # Collect counts
  counts <- unlist(lapply(1:trials, function(trial) {
    classifier_func(X, n_sub, cluster_prior, prop_algo)
  }))

  # Calculate histogram of counts
  gd <- hist(counts, breaks = seq(0, D), plot = FALSE)$density

  if (z_score) {
    X_shuffled <- shuffle_data(X)
    counts_shuffled <- unlist(lapply(1:trials, function(trial) {
      classifier_func(X_shuffled, n_sub, cluster_prior, prop_algo)
    }))
    g_shuffled <- hist(counts_shuffled, breaks = seq(0, D), plot = FALSE)$density
    scores <- (gd - mean(g_shuffled)) / sd(g_shuffled)
  } else {
    scores <- gd
  }

  # Create and return ClustVarSelect object
  create_clustvarselect(
    scores = scores,
    k_guess = k_guess,
    prop_algo = prop_algo,
    class_algo = class_algo,
    params = list(
      trials = trials,
      n_sub = n_sub,
      cluster_prior_min = cluster_prior_min,
      z_score = z_score,
      divergence = divergence,
      dim_reduction = dim_reduction
    )
  )
}




##### Helper functions

#' Randomly Shuffle Feature Values
#'
#' @description
#' Creates a copy of the input matrix where values within each feature (column)
#' are randomly shuffled while maintaining the marginal distribution of each feature.
#'
#' @param X Numeric matrix with dimensions (number of data points) x (number of features)
#'
#' @return A matrix of the same dimensions as X with randomly permuted values within
#'         each column
#'
#' @details
#' This function performs independent random permutations of values within each
#' feature column. It preserves the marginal distribution of each feature while
#' breaking any relationships between features. Useful for generating null
#' distributions in feature importance analyses.
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(20), ncol = 2)
#' X_shuffled <- shuffle_data(X)
#'
#' @export
shuffle_data <- function(X) {
  X_copy <- X
  for (feature in 1:ncol(X)) {
    X_copy[, feature] <- sample(X[, feature])
  }
  return(X_copy)
}


# incorporate power kmeans parameters

#' Find Important Dimensions Using Entropy-Based Classification
#'
#' @description
#' Identifies important dimensions in data using entropy-based classification with
#' either power k-means or hierarchical clustering.
#'
#' @param X Numeric matrix of data points
#' @param cell_sample Number of cells to sample
#' @param k_sub Vector of length 2 specifying range of clusters
#' @param cluster_algo Character specifying clustering algorithm ("kmeans" or "agglo")
#' @param pkm_params List of power k-means parameters including eta, initial_centers,
#'        divergence, convergence_threshold, seed, max_iter, and dim_reduction
#'
#' @return Vector of important feature indices
#'
#' @importFrom stats dist hclust cutree
#' @importFrom rpart rpart rpart.control
#'
#' @examples
#' X <- matrix(rnorm(1000), ncol = 20)
#' pkm_params <- list(
#'   eta = 1.05, initial_centers = NULL, divergence = "euclidean",
#'   convergence_threshold = 5, seed = 1, max_iter = 1000,
#'   dim_reduction = "pca"
#' )
#' dims <- find_classifier_dims_entropy(X, 30, c(3, 6), "kmeans", pkm_params)
#' @export
find_classifier_dims_entropy <- function(X, cell_sample, k_sub, cluster_algo, pkm_params) {
  set.seed(pkm_params$seed)
  cell_use <- sample(nrow(X), cell_sample)
  k_sub_i <- sample(k_sub[1]:k_sub[2], 1)
  Xit <- X[cell_use, ]

  k_guess <- if (cluster_algo == "kmeans") {
    power_kmeans_bregman(
      Xit,
      s = -0.5,
      k = k_sub_i,
      eta = pkm_params$eta,
      initial_centers = pkm_params$initial_centers,
      divergence = pkm_params$divergence,
      convergence_threshold = pkm_params$convergence_threshold,
      seed = pkm_params$seed,
      max_iter = pkm_params$max_iter,
      dim_reduction = pkm_params$dim_reduction
    )$classes
  } else {
    cutree(hclust(dist(Xit), method = "ward.D2"), k = k_sub_i)
  }
  k_guess <- as.integer(factor(k_guess, labels = 1:length(unique(k_guess))))

  # Input Validation
  if (is.null(k_guess) || any(is.na(k_guess))) {
    stop("Invalid cluster assignments")
  }
  if (length(unique(k_guess)) < 2) {
    stop("Clustering produced less than 2 distinct clusters")
  }

  gnout <- NULL
  for (i in 2:max(k_guess)) {
    for (j in 1:(i - 1)) {
      gnout <- c(gnout, one_tree(Xit, k_guess, i, j))
    }
  }
  return(gnout)
}


#' Find Important Dimensions Using Maximum Margin Classification
#'
#' @description
#' Identifies important dimensions in data using SVM-based maximum margin classification
#' with power k-means clustering.
#'
#' @param X Numeric matrix of data points
#' @param cell_sample Number of cells to sample
#' @param k_sub Vector of length 2 specifying range of clusters
#' @param cluster_algo Character specifying clustering algorithm (only "kmeans" supported)
#' @param pkm_params List of power k-means parameters including eta, initial_centers,
#'        divergence, convergence_threshold, seed, max_iter, and dim_reduction
#'
#' @return Vector of important feature indices
#'
#' @examples
#' X <- matrix(rnorm(1000), ncol = 20)
#' pkm_params <- list(
#'   eta = 1.05, initial_centers = NULL, divergence = "euclidean",
#'   convergence_threshold = 5, seed = 1, max_iter = 1000,
#'   dim_reduction = "pca"
#' )
#' dims <- find_classifier_dims_maxmargin(X, 30, c(3, 6), "kmeans", pkm_params)
#' @export
find_classifier_dims_maxmargin <- function(X, cell_sample, k_sub, cluster_algo, pkm_params) {
  set.seed(pkm_params$seed)
  cell_use <- sample(nrow(X), cell_sample)
  k_sub_i <- sample(k_sub[1]:k_sub[2], 1)
  Xit <- X[cell_use, ]

  k_guess <- power_kmeans_bregman(
    Xit,
    s = -0.5,
    k = k_sub_i,
    eta = pkm_params$eta,
    initial_centers = pkm_params$initial_centers,
    divergence = pkm_params$divergence,
    convergence_threshold = pkm_params$convergence_threshold,
    seed = pkm_params$seed,
    max_iter = pkm_params$max_iter,
    dim_reduction = pkm_params$dim_reduction
  )$classes
  k_guess <- as.integer(factor(k_guess, labels = 1:length(unique(k_guess))))

  # Input Validation
  if (is.null(k_guess) || any(is.na(k_guess))) {
    stop("Invalid cluster assignments")
  }
  if (length(unique(k_guess)) < 2) {
    stop("Clustering produced less than 2 distinct clusters")
  }

  gnout <- NULL
  for (i in 2:max(k_guess)) {
    for (j in 1:(i - 1)) {
      gnout <- c(gnout, one_plane(Xit, k_guess, i, j))
    }
  }
  return(gnout)
}


# decision tree based feature importance (entropy)

#' Identify Important Feature Using Decision Tree
#'
#' @description
#' Finds the most important feature distinguishing between two clusters using a
#' single-level decision tree.
#'
#' @param Xit Numeric matrix of data points
#' @param k_guess Vector of cluster assignments
#' @param ik1 First cluster index
#' @param ik2 Second cluster index
#'
#' @return Index of most important feature (or random feature if no split found)
#'
#' @details
#' Uses rpart to build a decision tree with maximum depth of 1 and returns the
#' feature with highest importance. If no split is found, returns a random feature.
#'
#' @importFrom rpart rpart rpart.control
#'
#' @examples
#' X <- matrix(rnorm(100), ncol = 5)
#' clusters <- sample(1:2, 20, replace = TRUE)
#' feature <- one_tree(X, clusters, 1, 2)
#' @export
one_tree <- function(Xit, k_guess, ik1, ik2) {
  # Subset data for the two clusters
  mask <- k_guess == ik1 | k_guess == ik2
  XTT <- Xit[mask, ]
  KTT <- k_guess[mask]

  # Train decision tree
  tree <- rpart(KTT ~ .,
    data = as.data.frame(XTT),
    method = "class",
    control = rpart.control(maxdepth = 1)
  )

  # Get feature importance
  importance <- tree$variable.importance
  if (length(importance) > 0) {
    return(which.max(importance))
  } else {
    return(sample(ncol(Xit), 1)) # Fallback if no split found
  }
}


# svm based feature importance (maxmargin)

#' Identify Important Feature Using L1-Penalized SVM
#'
#' @description
#' Finds the most important feature distinguishing between two clusters using
#' L1-penalized Support Vector Machine.
#'
#' @param Xit Numeric matrix of data points
#' @param k_guess Vector of cluster assignments
#' @param ik1 First cluster index (not used directly)
#' @param ik2 Second cluster index (not used directly)
#'
#' @return Index of feature with largest absolute coefficient
#'
#' @details
#' Uses penalizedSVM package to fit an L1-penalized SVM model and returns
#' the feature corresponding to the largest absolute coefficient.
#'
#' @importFrom penalizedSVM svmfs
#'
#' @examples
#' X <- matrix(rnorm(100), ncol = 5)
#' clusters <- sample(1:2, 20, replace = TRUE)
#' feature <- one_plane(X, clusters, 1, 2)
#' @export
one_plane <- function(Xit, k_guess, ik1, ik2) {
  # Subset data for the two clusters
  mask <- k_guess == ik1 | k_guess == ik2
  XTT <- Xit[mask, ]
  KTT <- factor(k_guess[mask], labels = c(-1, 1))
  colnames(XTT) <- 1:ncol(XTT)
  # Use penalizedSVM package to impose L1 penalty
  model <- penalizedSVM::svmfs(
    XTT,
    y = KTT,
    fs.method = "1norm",
    inner.val.method = "cv",
    cross.inner = 5,
    grid.search = "discrete",
    lambda1.set = c(0.01, 0.1, 0.5, 1, 10),
    parms.coding = "none",
    verbose = FALSE
  )

  # Get coefficients
  return(as.integer(names(which.max(abs(model$model$w)))))
}
