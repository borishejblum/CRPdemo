#' Sample from a inverse-Wishart distribution
#' 
#'@param n degrees of freedom
#'
#'@param lambda scale parameter
#'
riW <- function(n,lambda){
    
    p <- ncol(lambda)
    
    S <- try(solve(rW(n = n, Sigma = solve(lambda))), silent=TRUE)
    
    if(inherits(S, "try-error")){
        S <- solve(rW(n = n, Sigma = solve((lambda + diag(p)))))
    }
    
    return(S)
}
