#' Sparse Manifold Decomposition (SMD) with Power K-means
#'
#' @description
#' Performs sparse manifold decomposition using power k-means clustering with Bregman
#' divergence for robust feature importance analysis in high-dimensional data, especially
#' for scRNA-seq data. Features are normalized before analysis and importance scores can
#' be standardized relative to shuffled data. Supports parallel processing for improved performance.
#' 
#' @importFrom parallel mclapply detectCores
#' @importFrom stats median predict quantile rnorm sd as.formula dist hclust cutree
#' @importFrom graphics hist
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
#' @param n_cores Number of CPU cores to use for parallel processing. Default: min(detectCores() - 1, 2)
#'
#' @details
#' The function implements sparse manifold decomposition (SMD) combining power k-means clustering
#' with Bregman divergence for robust clustering. This method is particularly effective for
#' discovering discriminating features in high-dimensional biological data. For entropy-based
#' classification, it uses decision trees (rpart). For maxmargin classification, it uses
#' L1-penalized SVM.
#'
#' The function leverages parallel processing through the parallel package to improve
#' computational efficiency. The number of cores can be specified via the n_cores parameter.
#' If not specified, it defaults to using one less than the available cores (up to a maximum
#' of 2 cores to be conservative). On Windows systems, parallel processing defaults to
#' sequential processing due to platform limitations.
#'
#' @return A ClustVarSelect object containing:
#' \itemize{
#'   \item scores: Vector of feature importance scores
#'   \item k_guess: Number of clusters used
#'   \item prop_algo: Clustering algorithm used
#'   \item class_algo: Classification algorithm used
#'   \item params: List of parameters used including trials, n_sub, and other settings
#' }
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
#' # Basic usage with default parallel processing
#' result <- SMD_parallel(X, k_guess = 3)
#'
#' # Specify number of cores for parallel processing
#' result <- SMD_parallel(X, k_guess = 3, n_cores = 4)
#'
#' # With custom parameters
#' result <- SMD_parallel(X,
#'   k_guess = 3,
#'   prop_algo = "kmeans",
#'   divergence = "euclidean",
#'   dim_reduction = "pca",
#'   n_cores = 2
#' )
#'
#' @export
SMD_parallel <- function(X,
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
                dim_reduction = "pca", 
                n_cores = NULL) {
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
  
  # Parallel processing function for trials
  run_parallel_trials <- function(trial_indices) {
    # Set seed for reproducibility within each parallel process
    set.seed(seed + as.integer(trial_indices[1]))
    unlist(lapply(trial_indices, function(i) {
      classifier_func(X, n_sub, cluster_prior, prop_algo)
    }))
  }
  
  # Set default parameters
  if (is.null(n_sub)) n_sub <- as.integer(N * 0.8)
  if (is.null(trials)) trials <- as.integer(2 * D)
  if (is.null(cluster_prior_min)) cluster_prior_min <- 3
  if (is.null(n_cores)) {
    n_cores <- min(parallel::detectCores() - 1, 2)  # Conservative default
  } else {
    # Ensure n_cores doesn't exceed system limits or reasonable bounds
    max_cores <- min(parallel::detectCores(), 2)
    if (n_cores > max_cores) {
      warning(sprintf("Requested %d cores exceeds limit. Using %d cores instead.", 
                      n_cores, max_cores))
      n_cores <- max_cores
    }
  }
  
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
  
  # Split trials into chunks for parallel processing
  trial_chunks <- split(1:trials, cut(1:trials, n_cores, labels = FALSE))
  
  # Collect counts
  counts <- unlist(parallel::mclapply(trial_chunks, 
                                      run_parallel_trials, 
                                      mc.cores = n_cores))
  
  counts <- unlist(lapply(1:trials, function(trial) {
    classifier_func(X, n_sub, cluster_prior, prop_algo)
  }))
  
  # Calculate histogram of counts
  gd <- hist(counts, breaks = seq(0, D), plot = FALSE)$density
  
  if (z_score) {
    X_shuffled <- shuffle_data(X)
    # Parallel processing for shuffled data
    counts_shuffled <- unlist(parallel::mclapply(trial_chunks, 
                                                 function(chunk) {
                                                   set.seed(seed + 1000 + as.integer(chunk[1]))
                                                   unlist(lapply(chunk, function(i) {
                                                     classifier_func(X_shuffled, n_sub, 
                                                                     cluster_prior, prop_algo)
                                                   }))
                                                 }, 
                                                 mc.cores = n_cores))
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


