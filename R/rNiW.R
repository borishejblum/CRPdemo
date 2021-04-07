#' Sample from a Normal inverse-Wishart distribution
#' whose parameter are given by the structure hyper
#'
#'@param hyper hyper parameter list
#'
#'@param diagVar whether to use a diagonal Covariance matrix.
#'Default is \code{TRUE}.
#'
#'@importFrom stats rgamma rnorm

rNiW <- function(hyper, diagVar = TRUE){
  
  mu0 = hyper[["mu"]]
  kappa0 = hyper[["kappa"]]
  nu0 = hyper[["nu"]]
  lambda0 = hyper[["lambda"]]
  
  # Sample S from an inverse Wishart distribution
  if(diagVar){
    betas <- diag(lambda0)
    S <- diag(1/stats::rgamma(n=length(mu0), shape=nu0,
                              rate=betas))
  }else{
    S = riW(n = nu0, lambda = lambda0)
  }
  
  # Sample mu from a normal distribution
  muSupp <- matrix(stats::rnorm(length(mu0)), nrow=1,
                   ncol=length(mu0))%*%chol(S/kappa0)
  mu = mu0 + as.vector(muSupp)
  
  
  return(list("S"=S, "mu"=mu))
}
