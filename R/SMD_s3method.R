#' Print Method for ClustVarSelect Objects
#'
#' @param x A ClustVarSelect object returned by SMD function
#' @param ... Additional arguments passed to print
#'
#' @return Invisibly returns the input object
#'
#' @export
print.ClustVarSelect <- function(x, ...) {
  cat("Feature Selection Results via Power K-means Clustering\n\n")
  
  # Print basic information
  cat(sprintf("Number of features evaluated: %d\n", length(x$scores)))
  cat(sprintf("Number of clusters used: %d\n", x$k_guess))
  
  # Print top features
  cat("\nTop 10 features by importance score:\n")
  sorted_scores <- sort(x$scores, decreasing = TRUE)
  top_features <- head(sorted_scores, 10)
  for(i in seq_along(top_features)) {
    cat(sprintf("%d. Feature %d: %.4f\n", 
                i, 
                which(x$scores == top_features[i]), 
                top_features[i]))
  }
  
  cat("\nClustering method:", x$prop_algo)
  cat("\nClassification method:", x$class_algo)
  
  invisible(x)
}

#' Summary Method for ClustVarSelect Objects
#'
#' @param object A ClustVarSelect object returned by SMD function
#' @param ... Additional arguments passed to summary
#'
#' @return A list containing summary statistics
#'
#' @export
summary.ClustVarSelect <- function(object, ...) {
  # Calculate summary statistics
  scores <- object$scores
  stats <- list(
    mean_score = mean(scores),
    sd_score = sd(scores),
    median_score = median(scores),
    quantiles = quantile(scores, probs = c(0.25, 0.75)),
    n_features = length(scores),
    n_positive = sum(scores > 0),
    top_features = which(scores > quantile(scores, 0.9))
  )
  
  # Print summary
  cat("Summary of Feature Selection Results\n")
  cat("===================================\n\n")
  
  cat("Feature Score Statistics:\n")
  cat(sprintf("Mean score: %.4f\n", stats$mean_score))
  cat(sprintf("Standard deviation: %.4f\n", stats$sd_score))
  cat(sprintf("Median score: %.4f\n", stats$median_score))
  cat(sprintf("1st Quartile: %.4f\n", stats$quantiles[1]))
  cat(sprintf("3rd Quartile: %.4f\n", stats$quantiles[2]))
  
  cat("\nFeature Counts:\n")
  cat(sprintf("Total features: %d\n", stats$n_features))
  cat(sprintf("Features with positive scores: %d\n", stats$n_positive))
  
  cat("\nTop Features (90th percentile):", 
      paste(stats$top_features, collapse = ", "), "\n")
  
  # Return statistics invisibly
  invisible(stats)
}

#' Create ClustVarSelect Object
#'
#' @param scores Numeric vector of feature importance scores
#' @param k_guess Number of clusters used
#' @param prop_algo Clustering algorithm used
#' @param class_algo Classification algorithm used
#' @param params List of additional parameters used
#'
#' @return A ClustVarSelect object
#'
#' @keywords internal
create_clustvarselect <- function(scores, k_guess, prop_algo, class_algo, params = list()) {
  structure(
    list(
      scores = scores,
      k_guess = k_guess,
      prop_algo = prop_algo,
      class_algo = class_algo,
      params = params
    ),
    class = "ClustVarSelect"
  )
}