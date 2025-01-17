#'multivariate  Normal  probability density function
#'
#'@param x a p x n data matrix with n the number of observations and
#'p the number of dimensions
#'
#'@param mean mean vector
#'
#'@param varcovM variance-covariance matrix o
#'

mvnpdf <- function(x, mean, varcovM){
    if(!is.matrix(x)){
        stop("x should be a matrix")
    }
    
    n <- dim(x)[1]
    p <- dim(x)[2]
    
    if(is.null(mean)){
        mean <- rep(0, p)
        x0 <- x
    } else if(is.vector(mean) & length(mean)==p){
        x0 <- x-mean
    } else{
        stop("wrong input for mean")
    }
    
    
    if(is.null(varcovM)){
        varcovM <- diag(1, nrow=p, ncol=p)
    }
    if(dim(varcovM)[1]!=dim(varcovM)[2] | dim(varcovM)[1]!=p){
        stop("varcovM is not a square matrix or is of the wrong size")
    }
    
    R <- chol(varcovM)
    xRinv <- x0%*%solve(R)
    logSqrtDetvarcovM <- sum(log(diag(R)))
    
    quadform <- xRinv%*%t(xRinv)
    y <- exp(-0.5*quadform - logSqrtDetvarcovM -p*log(2*pi)/2)
    
    return(y)
    
}