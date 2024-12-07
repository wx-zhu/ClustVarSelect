# Different Bregman divergence functions

#' Euclidean Bregman Divergence Computation
#'
#' @description
#' Computes Euclidean Bregman divergence components for a given matrix
#'
#' @param X Numeric matrix of data points
#' @param out Character specifying output type: "f" for function value or "df" for derivative
#'
#' @return
#' If out="f", returns vector of squared L2 norms.
#' If out="df", returns matrix of derivatives (2X).
#'
#' @examples
#' X <- matrix(1:6, nrow=2)
#' euclidean_bregman(X, "f")
#' euclidean_bregman(X, "df")
#'
#' @export
euclidean_bregman <- function(X, out = "f") {
  # Input validation
  if (is.null(X)) {
    stop("Input matrix X cannot be NULL")
  }
  if (!is.numeric(X)) {
    stop("Input matrix X must contain numeric values")
  }
  
  if (out == "f") {
    f <- rowSums(X^2)
    return(f)
  } else if (out == "df") {
    df <- 2 * X
    return(df)
  } else {
    stop("Invalid output type")
  }
}

#' Kullback-Leibler Bregman Divergence Computation
#'
#' @description
#' Computes KL Bregman divergence components for positive data
#'
#' @param X Numeric matrix of positive data points
#' @param out Character specifying output type: "f" for function value or "df" for derivative
#'
#' @return
#' If out="f", returns vector of x*log(x) sums.
#' If out="df", returns matrix of derivatives (log(x) + 1).
#'
#' @details
#' Values less than 1e-10 are replaced with 1e-10 to prevent log(0).
#'
#' @examples
#' X <- matrix(abs(rnorm(6)), nrow=2)
#' kl_bregman(X, "f")
#' kl_bregman(X, "df")
#'
#' @export
kl_bregman <- function(X, out = "f") {
  # Input validation
  if (is.null(X)) {
    stop("Input matrix X cannot be NULL")
  }
  if (!is.numeric(X)) {
    stop("Input matrix X must contain numeric values")
  }
  
  if (any(X < 0)) {
    warning("Negative values in input matrix will be truncated to 1e-10")
  }
  
  # KL divergence (for positive data)
  X <- pmax(X, 1e-10)  # Prevent log(0)
  if (out == "f") {
    f <- rowSums(X * log(X))
    return(f)
  } else if (out == "df") {
    df <- log(X) + 1
    return(df)
  } else {
    stop("Invalid output type")
  }
}

#' Itakura-Saito Bregman Divergence Computation
#'
#' @description
#' Computes Itakura-Saito Bregman divergence components for positive data
#'
#' @param X Numeric matrix of positive data points
#' @param out Character specifying output type: "f" for function value or "df" for derivative
#'
#' @return
#' If out="f", returns vector of -log(x) sums minus dimension.
#' If out="df", returns matrix of derivatives (-1/x).
#'
#' @details
#' Values less than 1e-10 are replaced with 1e-10 to prevent division by zero.
#'
#' @examples
#' X <- matrix(abs(rnorm(6)), nrow=2)
#' itakura_saito_bregman(X, "f")
#' itakura_saito_bregman(X, "df")
#'
#' @export
itakura_saito_bregman <- function(X, out = "f") {
  # Input validation
  if (is.null(X)) {
    stop("Input matrix X cannot be NULL")
  }
  if (!is.numeric(X)) {
    stop("Input matrix X must contain numeric values")
  }
  
  if (any(X < 0)) {
    warning("Negative values in input matrix will be truncated to 1e-10")
  }
  
  # Itakura-Saito divergence (for positive data)
  X <- pmax(X, 1e-10)
  if (out == "f") {
    f <- -rowSums(log(X)) - ncol(X)
    return(f)
  } else if (out == "df") {
    df <- -1/X
    return(df)
  } else {
    stop("Invalid output type")
  }
}

#' Logistic Bregman Divergence Computation
#'
#' @description
#' Computes logistic Bregman divergence components
#'
#' @param X Numeric matrix of data points
#' @param out Character specifying output type: "f" for function value or "df" for derivative
#'
#' @return
#' If out="f", returns vector of log(1 + exp(x)) sums.
#' If out="df", returns matrix of derivatives (exp(x)/(1 + exp(x))).
#'
#' @examples
#' X <- matrix(rnorm(6), nrow=2)
#' logistic_bregman(X, "f")
#' logistic_bregman(X, "df")
#'
#' @export
logistic_bregman <- function(X, out = "f") {
  # Input validation
  if (is.null(X)) {
    stop("Input matrix X cannot be NULL")
  }
  if (!is.numeric(X)) {
    stop("Input matrix X must contain numeric values")
  }
  
  # Logistic loss
  if (out == "f") {
    f <- rowSums(log(1 + exp(X)))
    return(f)
  } else if (out == "df") {
    df <- exp(X)/(1 + exp(X))
    return(df)
  } else {
    stop("Invalid output type")
  }
}


#' Compute Pairwise Bregman Divergences
#'
#' @description 
#' Calculates pairwise Bregman divergences between two sets of points using
#' specified divergence type
#'
#' @param X First numeric matrix of data points
#' @param Y Second numeric matrix of data points
#' @param divergence Character specifying divergence type: "euclidean", "kl", 
#'        "itakura-saito", or "logistic"
#'
#' @return Matrix of pairwise Bregman divergences between X and Y
#'
#' @details
#' Computes D(x,y) = phi(x) - phi(y) - <nabla phi(y), x-y> where phi is determined
#' by the divergence type.
#'
#' @examples
#' X <- matrix(rnorm(6), nrow=2)
#' Y <- matrix(rnorm(9), nrow=3)
#' dist <- pairwise_bregman(X, Y, "euclidean")
#'
#' @seealso
#' \code{\link{euclidean_bregman}}, \code{\link{kl_bregman}},
#' \code{\link{itakura_saito_bregman}}, \code{\link{logistic_bregman}}
#'
#' @export
pairwise_bregman <- function(X, Y, divergence = "euclidean") {
  # Input validation
  if (!is.matrix(X) || !is.matrix(Y)) {
    stop("X and Y must be matrices")
  }
  if (ncol(X) != ncol(Y)) {
    stop("X and Y must have the same number of columns")
  }
  
  # Select the appropriate Bregman divergence function
  bregman_func <- switch(divergence,
                         "euclidean" = euclidean_bregman,
                         "kl" = kl_bregman,
                         "itakura-saito" = itakura_saito_bregman,
                         "logistic" = logistic_bregman,
                         stop("Unknown divergence type")
  )
  
  phi_X <- bregman_func(X)
  phi_Y <- bregman_func(Y)
  grad_Y <- bregman_func(Y, out = "df")
  
  pairwise_distances <- t(-outer(phi_Y, phi_X, "-") - 
                            (tcrossprod(grad_Y, X) - diag(tcrossprod(Y, grad_Y))))
  return(pairwise_distances)
}

