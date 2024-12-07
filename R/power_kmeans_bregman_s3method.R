#' Power K-means Result Object
#'
#' @param centers Matrix of cluster centers
#' @param classes Vector of cluster assignments
#' @param reduced_dim Dimension-reduced data matrix
#' @param params List of parameters used
#'
#' @return A PowerKmeans object
#'
#' @keywords internal
create_powerkmeans <- function(centers, classes, reduced_dim, params = list()) {
  structure(
    list(
      centers = centers,
      classes = classes,
      reduced_dim = reduced_dim,
      params = params
    ),
    class = "PowerKmeans"
  )
}

#' Print Method for PowerKmeans Objects
#'
#' @param x A PowerKmeans object
#' @param ... Additional arguments passed to print
#'
#' @return Invisibly returns the input object
#'
#' @export
print.PowerKmeans <- function(x, ...) {
  cat("Power K-means Clustering Results\n\n")
  
  # Basic information
  cat(sprintf("Number of clusters: %d\n", nrow(x$centers)))
  cat(sprintf("Number of samples: %d\n", length(x$classes)))
  cat(sprintf("Number of features: %d\n", ncol(x$centers)))
  
  # Cluster sizes
  cluster_sizes <- table(x$classes)
  cat("\nCluster sizes:\n")
  print(cluster_sizes)
  
  # Parameters used
  cat("\nParameters:\n")
  if (!is.null(x$params)) {
    cat(sprintf("Dimension reduction: %s\n", x$params$dim_reduction))
    cat(sprintf("Divergence type: %s\n", x$params$divergence))
  }
  
  invisible(x)
}

#' Summary Method for PowerKmeans Objects
#'
#' @param object A PowerKmeans object
#' @param ... Additional arguments passed to summary
#'
#' @return A list containing summary statistics
#'
#' @export
summary.PowerKmeans <- function(object, ...) {
  # Calculate summary statistics
  cluster_sizes <- table(object$classes)
  within_ss <- sapply(1:nrow(object$centers), function(i) {
    cluster_points <- object$reduced_dim[object$classes == i, , drop = FALSE]
    if(nrow(cluster_points) > 0) {
      sum(rowSums((cluster_points - matrix(object$centers[i,], 
                                           nrow = nrow(cluster_points), 
                                           ncol = ncol(cluster_points), 
                                           byrow = TRUE))^2))
    } else {
      0
    }
  })
  
  stats <- list(
    n_clusters = nrow(object$centers),
    n_samples = length(object$classes),
    n_features = ncol(object$centers),
    cluster_sizes = cluster_sizes,
    min_cluster_size = min(cluster_sizes),
    max_cluster_size = max(cluster_sizes),
    within_ss = within_ss,
    total_within_ss = sum(within_ss)
  )
  
  # Print summary
  cat("Summary of Power K-means Clustering\n")
  cat("==================================\n\n")
  
  cat("Clustering Statistics:\n")
  cat(sprintf("Number of clusters: %d\n", stats$n_clusters))
  cat(sprintf("Number of samples: %d\n", stats$n_samples))
  cat(sprintf("Number of features: %d\n", stats$n_features))
  
  cat("\nCluster Sizes:\n")
  print(stats$cluster_sizes)
  cat(sprintf("\nSmallest cluster: %d samples\n", stats$min_cluster_size))
  cat(sprintf("Largest cluster: %d samples\n", stats$max_cluster_size))
  
  cat("\nWithin Cluster Sum of Squares:\n")
  for(i in 1:length(stats$within_ss)) {
    cat(sprintf("Cluster %d: %.4f\n", i, stats$within_ss[i]))
  }
  cat(sprintf("\nTotal Within SS: %.4f\n", stats$total_within_ss))
  
  invisible(stats)
}

#' Plot Method for PowerKmeans Objects
#'
#' @param x A PowerKmeans object
#' @param ... Additional arguments passed to plot
#'
#' @return Invisibly returns the input object
#'
#' @importFrom stats prcomp
#' @importFrom graphics plot points legend
#' @export
plot.PowerKmeans <- function(x, ...) {
  if(ncol(x$reduced_dim) > 2) {
    # If more than 2 dimensions, use PCA for visualization
    pca <- prcomp(x$reduced_dim, scale. = TRUE)
    plot_data <- pca$x[,1:2]
  } else {
    plot_data <- x$reduced_dim
  }
  
  # Create plot
  plot(plot_data, 
       col = x$classes, 
       pch = 16,
       xlab = "First component",
       ylab = "Second component",
       main = "Power K-means Clustering Results",
       ...)
  
  # Add cluster centers
  if(ncol(x$centers) > 2) {
    centers_proj <- predict(pca, x$centers)[,1:2]
  } else {
    centers_proj <- x$centers
  }
  points(centers_proj, col = 1:nrow(x$centers), 
         pch = 8, cex = 2, lwd = 2)
  
  legend("topright", 
         legend = paste("Cluster", 1:nrow(x$centers)),
         col = 1:nrow(x$centers),
         pch = 16)
  
  invisible(x)
}