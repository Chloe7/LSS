#reference:: http://www.pfaffikus.de/files/books/springer/rex2/Rcode-2-1.R

#source("C:/Users/chloe/Documents/Research-2015/simulateFactorModel.r")

#install.packages("mvtnorm")
#install.packages("vars")
#install.packages("dse")

library(mvtnorm)
library(vars)
library(dse)

simulateFactorModel = function(n=10, r=2, T=100, fac.err.mean = rep(0, r), fac.err.sd = diag(r), Phi = list(matrix(c(0.7, 0, 0, -0.5),ncol=r)), Y.start = c(5, 1), obs.err.mean = rep(0, n), obs.err.sd = diag(n), Lambda) {    
	# the order of VAR
	s = length(Phi)

	fac.err = rmvnorm(T, mean=fac.err.mean, sigma=fac.err.sd) # T x r
	
	fac = matrix(0, T, r)
	fac[1,] = Phi[[1]] %*% Y.start + fac.err[1,]
	for(i in 2:(T)) {
	    fac[i,] = Phi[[1]] %*% fac[(i-1),] + fac.err[i,]
	}		
	colnames(fac) <- paste("factor-", seq(1:r), sep="")

    #--- Simulate observes X ---
    #--- Parameters for observation function ---

    #- True Lambda -
    obs.err = rmvnorm(T, mean=obs.err.mean, sigma=obs.err.sd)  # T x n
  
    obs = fac %*% t(Lambda) + obs.err
  
    return(list(X=obs, Y=fac))
}

