#'Gibbs Sampling with Alghorithm 2
#'
#'Chinese Restaurant Process representation
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
#'Necessary for proper display. Default is \code{0.08}s.
#'
#'@author Francois Caron, Boris Hejblum
#'
#'@importFrom viridis plasma
#'
#'@export gibbsDPMalgo2
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
#'  
#'  if(interactive()){
#'  # Gibbs sampler for Dirichlet Process Mixtures
#'  gifski::save_gif(expr = gibbsDPMalgo2(z, hyperG0, alpha = 0.0001, niter=15, 
#'                   doPlot=TRUE, pause = 0),
#'           delay = 0.02, width=600, height=500, gif_file = "CRPgibbs, alpha00001.gif")
#'  gifski::save_gif(expr = gibbsDPMalgo2(z, hyperG0, alpha = 2, niter=15, 
#'                   doPlot=TRUE, pause = 0),
#'           delay = 0.02, width=600, height=500, gif_file = "CRPgibbs, alpha2.gif")
#'  gifski::save_gif(expr = gibbsDPMalgo2(z, hyperG0, alpha = 30, niter=15, 
#'                   doPlot=TRUE, pause = 0),
#'           delay = 0.02, width=600, height=500, gif_file = "CRPgibbs, alpha30.gif")
#'           
#'           
#'  hyperG0[["lambda"]] <- 100*diag(2)
#'  gifski::save_gif(expr = gibbsDPMalgo2(z, hyperG0, alpha = 2, niter=15, 
#'                                        doPlot=TRUE, pause = 0),
#'           delay = 0.02, width=600, height=500, gif_file = "CRPgibbs, alpha2_lambda10.gif")
#' }
#'
gibbsDPMalgo2 <- function(z, hyperG0, alpha, niter, doPlot = TRUE, pause = 0.08){
    
    p <- dim(z)[1]
    n <- dim(z)[2]
    U_mu <- matrix(0, nrow=p, ncol=n)
    U_Sigma = array(0, dim=c(p, p, n))
    
    # U_SS is a list where each U_SS[k] contains the sufficient
    # statistics associated to cluster k
    U_SS <- list()
    
    m <- numeric(n)
    c <-numeric(n)
    
    # Initialisation: each observation is assigned to a different cluster----
    
    i <- 1
    for (k in 1:n){
        c[k] <- k
        U_SS[[c[k]]] <- update_SS(z=z[, k], S=hyperG0)
        NiW <- rNiW(U_SS[[c[k]]])
        U_mu[, c[k]] <- NiW[["mu"]]
        U_Sigma[, , c[k]] <- NiW[["S"]]
        m[c[k]] <- m[c[k]]+1
    }
    
    if(!doPlot){
        cat(i, "/", niter, " samplings\n", sep="")
    }
    
    
    for(i in 2:niter){
        for (k in 1:n){
            # Update cluster assignments c
            
            m[c[k]] <- m[c[k]] - 1
            U_SS[[c[k]]] <- downdate_SS(z[, k], U_SS[[c[k]]])
            c[k] <- sample_c(m, alpha, z[, k], hyperG0, U_mu, U_Sigma)
            m[c[k]] <- m[c[k]] + 1
            
            if(m[c[k]]>1){
                U_SS[[c[k]]] <- update_SS(z[, k], U_SS[[c[k]]])
            } else {
                U_SS[[c[k]]] <- update_SS(z[, k], hyperG0)
                NiW <- rNiW(U_SS[[c[k]]])
                U_mu[, c[k]] <- NiW[["mu"]]
                U_Sigma[, , c[k]] <- NiW[["S"]]
            }
            
            if(doPlot){
                plot_DPM2(z, U_mu, m, k, c, i, pause)
            }
        }
        
        # Update cluster locations U
        ind <- which(m!=0)
        
        for(j in 1:length(ind)){
            NiW <- rNiW(U_SS[[ind[j]]])
            U_mu[, ind[j]] <- NiW[["mu"]]
            U_Sigma[, , ind[j]] <- NiW[["S"]]
        }
        
        if(!doPlot){
            cat(i, "/", niter, " samplings\n", sep="")
        }
    }
    
    return(list("clusters" = c, "U_mu" = U_mu, "U_Sigma" = U_Sigma, 
                "partition" = m))
}








# Subfunctions ----
sample_c <- function(m, alpha, z, hyperG0, U_mu, U_Sigma){
    
    fullCl <- which(m!=0) # indexes of non empty clusters
    r <- sum(m)
    n <- numeric(length(fullCl))
    for (i in 1:length(fullCl)){
        n[i] <- mvnpdf(x = matrix(z, nrow= 1, ncol=length(z)) , 
                       mean = U_mu[, fullCl[i]], 
                       varcovM = U_Sigma[, , fullCl[i]])*m[fullCl[i]]  
    }
    
    n0 <- pred(z, hyperG0)
    const <- sum(n)/(alpha + r) + alpha/(alpha + r)*n0
    p0 <- alpha/(alpha + r)*n0/const # probability of sampling a new item
    
    u <- runif(n=1, min = 0, max = 1)
    if (u<p0){
        # Accept: allocate to a new cluster
        # cat("acceptation:", u, "<", p0, "\n")
        K <- which(m==0)[1]
    } 
    else{
        # Reject: allocate to a previously non empty cluster
        # cat("rejection:", u, ">=", p0, "\n")
        u1  <-  u - p0
        ind <- which(cumsum(n/const/(alpha+r))>rep(u1,length(n)))[1]
        K <- fullCl[ind]
    }
    return(K)
}


plot_DPM2 <- function(z, U_mu, m, k, c, i, pause = 0.08){
    n <- ncol(z)
    fullCl <- which(m!=0)
    U_mu2plot <- U_mu[, fullCl]    
    zClusters <- as.factor(c)
    levels(zClusters) <- as.character(1:length(levels(zClusters)))
    ncl <- length(levels(zClusters))
    
    mycols <- viridis::plasma(n=ncl)
    plot(t(z[,-k]), pch=16, col=mycols[zClusters[-k]],
         xlab = "Dim 1", ylab = "Dim 2", 
         main = paste0('Iteration ', i, ',  Obs. ', k, '\nNb. of non-empty clusters: ', ncl, 
                       "\nCRP Gibbs sampler"), 
         xlim=c(-3, 3), ylim=c(-3,3)
    )
    lines(t(z[,k]), pch=16, col=mycols[zClusters[k]], type="p", cex=2)
    lines(t(U_mu2plot), pch=4, col=mycols, type="p", cex=1.4)
    
    if(interactive()){
        Sys.sleep(pause)
    }
}
