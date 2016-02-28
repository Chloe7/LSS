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

use.simulateFactorModel() {
  n=6
  r=2
  T=100 
  fac.err.mean = rep(0, r)
  fac.err.sd = diag(r)
  Phi.true = list(matrix(c(0.9, 0, 0, 0.6),ncol=r))
  Y.start = c(0,0)
  obs.err.mean = rep(0, n)
  obs.err.sd = diag(n)
  set.seed(20141221)
  Lambda.true = matrix( runif(n*r, min=0, max=1) , nrow=n, ncol=r) 
  fac.mod = simulateFactorModel(n=n, r=r, T=T, fac.err.mean = fac.err.mean, fac.err.sd = fac.err.sd, Phi = Phi.true, Y.start = Y.start, obs.err.mean = obs.err.mean, obs.err.sd = obs.err.sd, Lambda.true)
  
  #--- Get factor Y and obs X --- 
  X = fac.mod$X
  Y = fac.mod$Y
  
  #--- Add column names ---
  Y.namestrings = paste0("Y", seq(1:r))
  dimnames(Y) = list(NULL, Y.namestrings)
  colMeans(Y)

  X.namestrings = paste0("X", seq(1:n))
  dimnames(X) = list(NULL, X.namestrings)
  
  #--- Plot ---
  par(mfrow=c(2,1))
  ts.plot(Y, lty=1:r, ylab="Factors") 
  legend(locator(1), legend=Y.namestrings, lty=1:r)

  ts.plot(X, col=1:n, ylab="Observations")
  legend(locator(1), legend=colnames(X), col=1:n)
}

