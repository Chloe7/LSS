#updated on 2016-01-29 on zero limit of tau.lo and tau.hi

#source("C:/Users/chloe/Documents/Research-2015/stiefel1.h.r")

library(MASS)

#TEST = TRUE

#--- Matrix Operations ---
is.square <- function(M) {
  return(is.matrix(M) && (nrow(M) == ncol(M)))
 }
		   
matTrace <- function(M) {
  return(ifelse(is.square(M), sum(diag(M)), NA))
}

svd_inv <- function(M) {
    M.SVD <- svd(M)
    return( with(M.SVD, v %*% diag(1/d) %*% t(u)) )   # calculate inverse with SVD
}

### loss as a function of Y
loss_func <- function(X, Y, w=1) {  
  n = ncol(X);n
  r = ncol(Y);r
  T = nrow(X);T
  if(T != nrow(Y)) {stop("function: loss_func(): X and Y not compatible!")}
  
  #err = (diag(T)-Y %*% t(Y)) %*% X
  #matTrace(t(X) %*% t(diag(T)-Y %*% t(Y)) %*% (diag(T)-Y %*% t(Y)) %*% X)
  #matTrace(t(X) %*% (diag(T)-Y %*% t(Y)) %*% X)
  #matTrace( t(X) %*% X - t(X) %*% Y %*% t(Y) %*% X )
  
  tXY = t(X) %*% Y  
  L1 = sum(X * X) - sum(tXY * tXY) #matTrace( t(X) %*% X - t(X) %*% Y %*% t(Y) %*% X )
  G1 = - 2 * X %*% tXY
  
  #L1 = matTrace( t(X) %*% X - t(X) %*% Y %*% t(Y) %*% X )
  #G1 = - 2 * X %*% t(X) %*% Y
  
  A = rbind( matrix(0, 1, T),
           cbind( matrix(0, T-1, 1), diag(T-1)) )
		   
  B = rbind( cbind( matrix(0, T-1, 1), diag(T-1) ),
           matrix(0, 1, T) )
		   
  C = rbind( cbind( diag(T-1), matrix(0, T-1, 1) ), 
           matrix(0, 1, T) )
  E = diag(r)
  S11 = t(Y) %*% A %*% Y
  S01 = t(Y) %*% B %*% Y
  S00 = t(Y) %*% C %*% Y
  
  L2 = matTrace( t(Y) %*% A %*% Y )
  G2 = (A + t(A)) %*% Y 
  for(i in 1:r) {
    Ei = E[,i, drop=FALSE]
	#L2 = L2 - as.numeric(t(Ei)%*%t(S01)%*%Ei) * as.numeric(t(Ei)%*%S00%*%Ei)^(-1) * as.numeric(t(Ei)%*%S01%*%Ei)
	L2 = L2 - S01[i,i] * S00[i,i]^(-1) * S01[i,i]
	
	#G2 = G2 - 2 * ( as.numeric(t(Ei)%*%t(S01)%*%Ei) * as.numeric(t(Ei)%*%S00%*%Ei)^(-1) * (B + t(B)) %*% Y %*% Ei %*% t(Ei) - as.numeric(t(Ei)%*%t(S01)%*%Ei)^2 * as.numeric(t(Ei)%*%S00%*%Ei)^(-2) * C %*% Y %*% Ei %*% t(Ei) )
	G2 = G2 - 2 * ( S01[i,i] * S00[i,i]^(-1) * (B + t(B)) %*% Y %*% Ei %*% t(Ei) - S01[i,i]^2 * S00[i,i]^(-2) * C %*% Y %*% Ei %*% t(Ei) )	
  }
  L = L1 + w * L2
  G = G1 + w * G2
  
  return(list(L=L, G=G, L1=L1, L2=L2))
}


curve_template = function(Y, U, V) {
  function(tau) {
    return( Y - tau * U %*% svd_inv(diag(ncol(V)) + 0.5 * tau * t(V) %*% U) %*% t(V) %*% Y )
  }
}

#v_curve_template <- Vectorize(curve_template, c("tau"))
#v_curve_func <- Vectorize(curve_func, c("tau"))

zoom <- function(tau.lo, tau.hi, X, Y, G, U, V, w=1, rho1=1e-4, rho2=0.9) {
  ## Interpolate (using quadratic, cubic, or bisection) to find a trial step length between tau.lo and tau.hi
  ## use bisection here  

}

# Linear search part for step size is based on Algorithm 3.5 (Line Search Algorithm)
# rho1=1e-4; rho2=0.9; tau.max=0.01
curvilinear_search <- function(X, Y, G, U, V, w=1, k, rho1=1e-4, rho2=0.9, tau.max=0.01) {
  # curve_func fixed through the following iterations
  # curve_func = curve_template(Y, U, V) 
  
  AA <- G %*% t(Y) - Y %*% t(G) 
  
  curve.val.at.0 <- Y   
  curve.diff.at.0 <- - AA %*% Y 
  loss.val.at.0 <- loss_func(X, curve.val.at.0, w=w)$L
  loss.diff.at.0 <- sum(G * curve.diff.at.0)    
  
  tau.hist <- c(0, 0.5 * tau.max)
  loss.val.hist <- c(loss.val.at.0)
  i <- 2
  while(TRUE) {
    inv.at.tau <- svd_inv(diag(ncol(V)) + 0.5 * tau.hist[i] * t(V) %*% U)
    val.at.tau <- U %*% inv.at.tau %*% t(V)
  
    curve.val.at.tau <- Y - tau.hist[i] * val.at.tau %*% Y	
    curve.diff.at.tau <- -(diag(nrow(U)) - 0.5 * tau.hist[i] * val.at.tau ) %*% AA %*% (Y + curve.val.at.tau)  # reformulate to avoid re-calculation
    loss.val.hist <- c(loss.val.hist, loss_func(X, curve.val.at.tau, w=w)$L)
    loss.diff.at.tau <- sum(G * curve.diff.at.tau)  # matTrace(t(G) %*% curve.diff.at.tau)
	
	if(loss.val.hist[i] > loss.val.at.0 + rho1 * tau.hist[i] * loss.diff.at.0 | (loss.val.hist[i] >= loss.val.hist[i-1] & i>2)) {
	  #tau.opt <- zoom(tau.hist[i-1], tau.hist[i])
	  { ## zoom(tau.hist[i-1], tau.hist[i])
	    tau.lo = tau.hist[i-1]
        tau.hi = tau.hist[i]
	    if(TEST) {
	      cat(paste("tau history at Y search step", k, "-", i, ":: enter condition 1: zoom(tau.hist[i-1], tau.hist[i]), (tau.lo =", tau.lo, "tau.hi =", tau.hi,")\n")); flush.console()
        }        
		while(TRUE) {
		  ## trail tau
          tau.tr = 0.5 * (tau.lo + tau.hi) 
        
          inv.at.tau.tr <- svd_inv(diag(ncol(V)) + 0.5 * tau.tr * t(V) %*% U)
          val.at.tau.tr <- U %*% inv.at.tau.tr %*% t(V)
  
          curve.val.at.tau.tr <- Y - tau.tr * val.at.tau.tr %*% Y	
          curve.diff.at.tau.tr <- - (diag(nrow(U)) - 0.5 * tau.tr * val.at.tau.tr ) %*% AA %*% (Y + curve.val.at.tau.tr)  # reformulate to avoid re-calculation
          loss.val.at.tau.tr <- loss_func(X, curve.val.at.tau.tr, w=w)$L
          loss.diff.at.tau.tr <- sum(G * curve.diff.at.tau.tr)  # matTrace(t(G) %*% curve.diff.at.tau)      
        
		  ## values at tau.lo
		  curve.val.at.tau.lo <- Y - tau.lo * U %*% svd_inv(diag(ncol(V)) + 0.5 * tau.lo * t(V) %*% U) %*% t(V) %*% Y	
		  loss.val.at.tau.lo <- loss_func(X, curve.val.at.tau.lo, w=w)$L
		
          if(loss.val.at.tau.tr > loss.val.at.0 + rho1 * tau.tr * loss.diff.at.0 | (loss.val.at.tau.tr >= loss.val.at.tau.lo)) {
		    tau.hi <- tau.tr			
		    if(TEST) {
			  cat(paste("--- zoom(tau.hist[i-1], tau.hist[i]) condition 1 :: tau.hi <- tau.tr, (tau.lo =", tau.lo, "tau.hi =", tau.hi, ")\n")); flush.console()
            }		    
		  } else {
		    if ( abs(loss.diff.at.tau.tr) <= - rho2 * loss.diff.at.0 ) {
			  if(TEST) {
			    cat(paste("--- zoom(tau.hist[i-1], tau.hist[i]) condition 2 :: tau.opt <- tau.tr = ", tau.tr, "\n")); flush.console()
              }
		      tau.opt <- tau.tr
			  curve.val.at.tau.opt <- curve.val.at.tau.tr
			  break
		    } else if ( loss.diff.at.tau.tr * (tau.hi - tau.lo) >= 0 ) {
			  tau.hi <- tau.lo 
			  if(TEST) {
			    cat(paste("--- zoom(tau.hist[i-1], tau.hist[i]) condition 3 :: tau.hi <- tau.lo & tau.lo <- tau.tr, (tau.lo =", tau.lo, "tau.hi =", tau.hi, ")\n")); flush.console()
              }		      
		    }
			tau.lo <- tau.tr
		  }
		  
		  if(abs(tau.hi - tau.lo) < 1e-15) {
			  tau.opt <- tau.lo
			  curve.val.at.tau.opt <- curve.val.at.tau.lo
			  break
		  } ## end if(abs(tau.hi - tau.lo) < 1e-5)
		} ## end while(TRUE)		
	  } ## end zoom(tau.hist[i-1], tau.hist[i])
	  break
	} else if( abs(loss.diff.at.tau) <= -rho2 * loss.diff.at.0) {
	  if(TEST) {
	      cat(paste("tau history at Y search step", k, "-", i, ":: enter condition 2: tau.opt <- tau.hist[i] <- ", tau.hist[i],"\n")); flush.console()
        }
	  tau.opt <- tau.hist[i]
	  curve.val.at.tau.opt <- curve.val.at.tau
	  break
	} else if(loss.diff.at.tau >= 0) {
	  ##tau.opt <- zoom(tau.hist[i], tau.hist[i-1])
	  { ## zoom(tau.hist[i], tau.hist[i-1])
	    tau.lo = tau.hist[i]
        tau.hi = tau.hist[i-1]
		
	    if(TEST) {
	      cat(paste("tau history at Y search step", k, "-", i, ":: enter condition 3: zoom(tau.hist[i], tau.hist[i-1]), (tau.lo =", tau.lo, "tau.hi =", tau.hi, ")\n")); flush.console()
        }
        
		while(TRUE) {
		  ## trail tau
          tau.tr = 0.5 * (tau.lo + tau.hi) 
        
          inv.at.tau.tr <- svd_inv(diag(ncol(V)) + 0.5 * tau.tr * t(V) %*% U)
          val.at.tau.tr <- U %*% inv.at.tau.tr %*% t(V)
  
          curve.val.at.tau.tr <- Y - tau.tr * val.at.tau.tr %*% Y	
          curve.diff.at.tau.tr <- - (diag(nrow(U)) - 0.5 * tau.tr * val.at.tau.tr ) %*% AA %*% (Y + curve.val.at.tau.tr)  # reformulate to avoid re-calculation
          loss.val.at.tau.tr <- loss_func(X, curve.val.at.tau.tr, w=w)$L
          loss.diff.at.tau.tr <- sum(G * curve.diff.at.tau.tr)  # matTrace(t(G) %*% curve.diff.at.tau)      
        
		  ## values at tau.lo
		  curve.val.at.tau.lo <- Y - tau.lo * U %*% svd_inv(diag(ncol(V)) + 0.5 * tau.lo * t(V) %*% U) %*% t(V) %*% Y	
		  loss.val.at.tau.lo <- loss_func(X, curve.val.at.tau.lo, w=w)$L
		
          if(loss.val.at.tau.tr > loss.val.at.0 + rho1 * tau.tr * loss.diff.at.0 | (loss.val.at.tau.tr >= loss.val.at.tau.lo)) {
		    tau.hi <- tau.tr   		    
		    if(TEST) {
			  cat(paste("--- zoom(tau.hist[i], tau.hist[i-1]) condition 1 :: tau.hi <- tau.tr, (tau.lo =", tau.lo, "tau.hi =", tau.hi, ")\n")); flush.console()
            }		    
		  } else {
		    if ( abs(loss.diff.at.tau.tr) <= -rho2 * loss.diff.at.0 ) {
			  tau.opt <- tau.tr
			  curve.val.at.tau.opt <- curve.val.at.tau.tr
			  
			  if(TEST) {
			    cat(paste("--- zoom(tau.hist[i], tau.hist[i-1]) condition 2 :: tau.opt <- tau.tr = ", tau.tr," \n")); flush.console()
              }			  
			  break
		    } else if ( loss.diff.at.tau.tr * (tau.hi - tau.lo) >= 0 ) {			  
			  tau.hi <- tau.lo			  			  
			  if(TEST) {
			    cat(paste("--- zoom(tau.hist[i], tau.hist[i-1]) condition 3 :: tau.hi <- tau.lo & tau.lo <- tau.tr, (tau.lo =", tau.lo, "tau.hi =", tau.hi, ")\n")); flush.console()
              }		      
		    }
			tau.lo <- tau.tr
		  }
		  if(abs(tau.hi - tau.lo) < 1e-15) {
			  tau.opt <- tau.lo
			  curve.val.at.tau.opt <- curve.val.at.tau.lo
			  break
		  } ## end if(abs(tau.hi - tau.lo) < 1e-5)
		} ## end while(TRUE)			
	  }  ## end zoom(tau.hist[i], tau.hist[i-1])
	  break
	}
	
	tau.hist <- c(tau.hist, 0.5 * (tau.hist[i]+tau.max))
	
	if(TEST) {
	  cat(paste("tau history at Y search step", k, "-", i, "::"), paste0(tau.hist," ,"), "\n"); flush.console()
    }
	i <- i + 1	
  } ## end outer while(TRUE)
  
  #---------   test  -----------
  if(TEST.PLOT) { 
    xx = seq(0,tau.max, length.out=100)
    yy = NULL
    for( x in xx ) {
      #yy = c(yy, loss_func(X, Y - x * U %*% solve(diag(ncol(V)) + 0.5 * x * t(V) %*% U) %*% t(V) %*% Y, w=w)$L)
	  aaSVD = svd(diag(ncol(V)) + 0.5 * x * t(V) %*% U)
	  aa_inv = with(aaSVD, v %*% diag(1/d) %*% t(u))
      yy = c(yy, loss_func(X, Y - x * U %*% aa_inv %*% t(V) %*% Y, w=w)$L)
    }  
	jpeg(file=paste(getwd(), "/test/iter_Y_", k,".jpg",sep=""))
    plot(xx, yy, lty=2, xlab="tau", ylab="loss value") 
	abline(v=tau.opt)
	title(paste("Loss function v.s. tau at Y_", k,sep=""))
	dev.off()
  }
  #-----------------------------
  
  return(curve.val.at.tau.opt) # cbind(curve.val.at.tau, Y)
}
