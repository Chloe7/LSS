#source("C:/Users/chloe/Documents/Research-2015/simXandY.orthon.r")

#install.packages("mvtnorm")
library(mvtnorm)

simXandY = function(n=10, r=2, T=100, fac.err.mean = rep(0, r), fac.err.sd = diag(r), Phi1 = matrix(c(0.7, 0, 0, -0.5), r, r), Y.start = c(5, 1), obs.err.mean = rep(0, n), obs.err.sd = diag(n), Lambda) {
  #--- Simulate factors Y ---
  #--- parameters for factor dynamics ---
  #- ft = phi1 * ft-1 + et

  #n=10; r=2; T=100
  e = rmvnorm(T, mean=fac.err.mean, sigma=fac.err.sd)
  e = t(e); #e

  Y = matrix(0, r, T)
  Y[,1] = Phi1 %*% Y.start + e[,1]

  for(t in 2:T) {
    Y[,t] = Phi1 %*% Y[,(t-1)] + e[,t]
  }
  Y = t(Y)
#--- Simulate observes X ---
#--- Parameters for observation function ---

#- True Lambda -
  eps = rmvnorm(T, mean=obs.err.mean, sigma=obs.err.sd); 
  
  X = Y %*% t(Lambda) + eps
  
  return(list(X=X, Y=Y, Lambda=Lambda, Phi1=Phi1, fac.err.mean=fac.err.mean, fac.err.sd=fac.err.sd, obs.err.mean=obs.err.mean, obs.err.sd=obs.err.sd, Y.start=Y.start))
}

