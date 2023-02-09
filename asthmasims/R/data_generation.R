# %%%%% Functions used for data generation %%%%%


#' @title generate_data
#' @description Generate data where each participant is a single row.
#' @param alpha Intercept parameters
#' @param eta Linear predictor for each row
#' @return Simulated data y of length length(eta)
generate_data_eta <- function(alpha, eta) {
  n <- length(eta)
  p <- ordinal_probs_mats(alpha, eta)
  return(extraDistr::rcat(n, p))
}


#' @title generate_data_X
#' @description Generate data where each participant is a single row.
#' @param X Design matrix
#' @param alpha Intercept parameters
#' @param beta Design matrix parameters
#' @return Simulated data y of length nrow(X)
generate_data_X <- function(X, alpha, beta) {
  K <- length(alpha) + 1
  eta <- as.vector(X %*% beta)
  return(generate_data_eta(alpha, eta))
}


#' @title generate_data_x
#' @description Generate data where each participant assignment is x.
#' @param x Participant intervention indicators
#' @param Xdes Design matrix
#' @param alpha Intercept parameters
#' @param beta Coefficient parameters
#' @return Simulated dataset
generate_data_x <- function(x, Xdes, alpha, beta) {
  X <- Xdes[x, ]
  return(generate_data_X(X, alpha, beta))
}


#' @title generate_data_agg
#' @description Generate aggregated data where each row is counts across a covariate pattern.
#' @param n Sample size
#' @param alpha Intercept parameters
#' @param eta Group specific shift on alpha
#' @return Simulated data matrix y of dimension lenth(eta) x length(alpha) + y
generate_data_agg <- function(n, alpha, eta) {
  K <- length(alpha) + 1
  N <- length(eta)
  p <- ordinal_probs_mats(alpha, eta)
  y <- t(sapply(seq_len(N), function(z) stats::rmultinom(1, n[z], p[z, ])))
  return(y)
}


#' @title Generate ordinal data in aggregate
#' @description Simulate data using p instead of alpha and eta/beta
#' @param n The sample size per arm
#' @param p The outcome level probabilities (columns) per arm (rows)
#' @return A vector same dimension as p with outcome level counts by arm
generate_data_p <- function(n, p) {
  N <- nrow(p)
  y <- t(sapply(seq_len(N), function(z) stats::rmultinom(1, n[z], p[z, ])))
  return(y)
}
