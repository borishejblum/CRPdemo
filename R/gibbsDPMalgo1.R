#'Gibbs Sampling with Alghorithm 1 from Neal
#'
#'Polya Urn representation
#'
#'@param z data matrix \code{p x n}
#'
#'@param hyperG0 a list of hyper parameters
#'
#'@param alpha concentration parameter
#'
#'@param niter number of MCMC iterations
#'
#'@param doPlot whether to print plots. Default is \code{TRUE}.
#'
#'@param pause how long to pause between graphical outputs. 
#'Necessary for proper display. Default is \code{0.1}s.
#'
#'@author Francois Caron, Boris Hejblum
#'
#'@importFrom viridis plasma
#'
#'@importFrom stats runif
#'
#'@importFrom graphics lines
#'
#'@export gibbsDPMalgo1
#'
#'@examples
#' rm(list=ls())
#' #Number of data
#' n <- 100
#' set.seed(1231)
#' 
#' # Sample data
#' m <- matrix(nrow=2, ncol=4, c(-1, 1, 0, 2, 1, -2, -1, -2))
#' p <- c(0.2, 0.1, 0.4, 0.3) # frequence des clusters
#' 
#' library(expm)
#' s <- array(dim=c(2,2,4))
#' s[, ,1] <- matrix(nrow=2, ncol=2, c(0.1, 0, 0, 0.1))
#' s[, ,2] <- matrix(nrow=2, ncol=2, c(0.01, 0, 0, 0.1))
#' s[, ,3] <- matrix(nrow=2, ncol=2, c(0.1, 0.08, 0.08, 0.1))
#' s[, ,4] <- .1*diag(2)
#' c <- rep(0,n)
#' z <- matrix(0, nrow=2, ncol=n)
#' for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  z[,k] <- m[, c[k]] + expm::sqrtm(s[, , c[k]])%*%matrix(rnorm(2, mean = 0, sd = 1), nrow=2, ncol=1)
#' }
#'  
#'  # Set parameters of G0
#'  hyperG0 <- list()
#'  hyperG0[["mu"]] <- c(0,0)
#'  hyperG0[["kappa"]] <- 1
#'  hyperG0[["nu"]] <- 4
#'  hyperG0[["lambda"]] <- diag(2)
#'  # Scale parameter of DPM
#'  alpha <- 3
#'  # Number of iterations
#'  N <- 30 
#'  # do some plots
#'  doPlot <- TRUE 
#'  
#'  plot(t(z), pch=16, xlim=c(-3,3), ylim=c(-3,3), xlab="Dim 1", ylab="Dim 2")
#'  # Gibbs sampler for Dirichlet Process Mixtures
#'  if(interactive()){
#'    GibSample <- gibbsDPMalgo1(z, hyperG0, alpha, niter = N, doPlot)
#'  gifski::save_gif(expr = gibbsDPMalgo1(z, hyperG0, alpha = 2, niter=15, doPlot=TRUE, pause = 0),
#'           delay = 0.02, width=600, height=500, gif_file = "BMgibbs_alpha2.gif")
#'  }
#'
#'
gibbsDPMalgo1 <- function (z, hyperG0, alpha, niter, doPlot=TRUE, pause = 0.1){

    
    p <- dim(z)[1]
    n <- dim(z)[2]
    theta_mu <- matrix(0, nrow=p, ncol=n)
    theta_Sigma = array(0, dim=c(p, p, n))
    
    
    # Initialization----
    i=1
    hyper <- update_SS(z[, 1], hyperG0)
    NiW <- rNiW(hyper)
    theta_mu[, 1] <- NiW[["mu"]]
    theta_Sigma[, , 1] <- NiW[["S"]]
    for (k in 2:n){
        sampTheta <- sample_theta(alpha, z=z[, k], hyperG0=hyperG0, 
                                  theta_mu_notk=theta_mu[, 1:(k-1)], 
                                  theta_Sigma_notk=theta_Sigma[, , 1:(k-1)])
        theta_mu[, k] <- sampTheta[["mu"]]
        theta_Sigma[, , k] = sampTheta[["Sigma"]]
    }
    if(!doPlot){
        cat(i, "/", niter, " samplings\n", sep="")
    }
    
    for(i in 2:niter){
        for (k in 1:n){
            # Sample theta_k | theta_{-k}, z_k
            #ind_notk <- c(1:n)[-k]
            sampTheta <- sample_theta(alpha, z=z[, k], hyperG0, 
                                      theta_mu_notk=theta_mu[, -k], 
                                      theta_Sigma_notk=theta_Sigma[, , -k])
            theta_mu[, k] <- sampTheta[["mu"]]
            theta_Sigma[, , k] = sampTheta[["Sigma"]]
        
            if(doPlot){
                plot_DPM1(z, theta_mu, k, n, i, pause)
            }
        }
        
        if(!doPlot){
            cat(i, "/", niter, " samplings\n", sep="")
        }
    }
    return(list("mu"=theta_mu, "Sigma"=theta_Sigma))
}








# Subfunctions ----

sample_theta <- function(alpha, z, hyperG0, theta_mu_notk, theta_Sigma_notk){

    if(is.null(dim(theta_mu_notk))){
        p <- 1
        n <-mvnpdf(x = matrix(z, nrow= 1, ncol=length(z)) , mean = theta_mu_notk, varcovM = theta_Sigma_notk)
    }else{
        p <- dim(theta_mu_notk)[2]
        n <- numeric(p)
        for (i in 1:p){
            n[i] <- mvnpdf(x = matrix(z, nrow= 1, ncol=length(z)) , mean = theta_mu_notk[,i], varcovM = theta_Sigma_notk[, , i])  
        }
    }
    
    n0 <- pred(z, hyperG0)
    const <- alpha*n0+sum(n)
    p0 <- alpha*n0/const

    u <- runif(n=1, min = 0, max = 1)
    if (u<p0){
        # Accept: sample new value
        # cat("acceptation:", u, "<", p0, "\n")
        hyper <- update_SS(z, hyperG0)
        NiW <- rNiW(hyper)
        theta_mu_k <- NiW[["mu"]]
        theta_Sigma_k <- NiW[["S"]]
    } 
    else{
        # Reject: sample old value
        # cat("rejection:", u, ">=", p0, "\n")
        u1  <-  u - p0
        ind <- match(TRUE, cumsum(n/const) > u1)
        
        if(is.null(dim(theta_mu_notk))){
            theta_mu_k = theta_mu_notk  
            theta_Sigma_k = theta_Sigma_notk          
        } else{
            theta_mu_k = theta_mu_notk[,ind]
            theta_Sigma_k = theta_Sigma_notk[, , ind] 
        }  
    }
    return(list("mu"=theta_mu_k, "Sigma"= theta_Sigma_k))
}


plot_DPM1 <- function(z, theta_mu, k, n, i, pause = 0.1){
    ind <- unique(theta_mu[1,])
    ncl <- length(ind)
    cl <- numeric(n)
    U_mu <- matrix(0, nrow=2, ncol=length(ind))
    
    for(j in 1:ncl){
        ind2 <- which(theta_mu[1, ]==ind[j])
        cl[ind2] <- j
        U_mu[, j] <- theta_mu[, ind2[1]]
    }
    mycols <- viridis::plasma(n=ncl)
    plot(t(z[,-k]), pch=16, col=mycols[cl[-k]],
         xlab = "Dim 1", ylab = "Dim 2", 
         main = paste0('Iteration ', i, ',  Obs. ', k, '\nNb. of non-empty clusters: ', ncl,
                       "\nPolya urn Gibbs sampler"), 
         xlim=c(-3, 3), ylim=c(-3,3)
    )
    lines(t(z[,k]), pch=16, col=mycols[cl[k]], type="p", cex=2)
    lines(t(U_mu), pch=4, col=mycols, type="p", cex=1.4)
    if(interactive()){
        Sys.sleep(pause)
    }
}
