#setwd("C:/Users/chloe/Documents/Research-2015/current")
#source("main_data_v6.r")

#install.packages("Matrix")
library(Matrix)
library(vars)
library(MASS)


#source file in the same directory
source(paste0(getwd(),"/stiefel.v5.h.r"))
source(paste0(getwd(),"/ls.stiefel.h.r"))
source(paste0(getwd(),"/dfm.h.r"))
source(paste0(getwd(),"/genPosDefMat.r"))
source(paste0(getwd(),"/simulateFactorModel.r"))
#source(paste0(getwd(),"/simXandY.orthon.r"))
source(paste0(getwd(),"/sink.reset.r"))
source(paste0(getwd(),"/genData.r"))

TEST <- TRUE
TEST.PLOT <- FALSE

#--- Xialu paper P8
#--- Orthonormal matrix H1 and H2
spaceDist <- function(H1, H2) {
  div = max(ncol(H1), ncol(H2))
  return(sqrt(1- 1/div * matTrace(H1 %*% t(H1) %*% H2 %*% t(H2))))
}

#--- Orthonormal basis 
orthonormal_basis <- function(M) {
  return( qr.Q(qr(M)) )
}

directory <- function ( d ) {
  path <- file.path(getwd(), d)
  dir.create(file.path(getwd(), d), showWarnings = FALSE)
  return(path)  
}

generateMCData <- function( T.vec=c(100, 300, 500), 
                            n.vec=c(50, 100, 150, 200, 500, 700, 1000), 
							r=2, 
							n.mc=100, 
							fac.obs.noise.ratio = 1,
							Lambda.true = matrix( runif(n*r, min=0, max=1) , nrow=n, ncol=r),
							Phi.true = list(matrix(c(0.7, 0, 0, -0.2), r, r)),
                            Y.start = c(5, 1),
                            out.dir	= directory("data")					
						   ) {
   
  #T.vec <- c(100, 300) 
  #n.vec <- c(50, 100, 150, 200, 500, 1000)
  #n.vec <- c(200, 500, 1000)
  #r <- 2
  #n.mc <- 100  
  # fac.obs.noise.ratio = 0.1   # with w = 10 in ls.stiefel()
  # fac.obs.noise.ratio = 1     # with w = 1 in ls.stiefel()
  
  #--- generate data ---
  for(T in T.vec) {
    for(n in n.vec) {	 	  
	  # Factor Model Parameters
	  set.seed(20141221)
      Lambda.true = matrix( runif(n*r, min=0, max=1) , nrow=n, ncol=r) 
      #Phi.true = list(matrix(c(0.7, 0, 0, -0.2), r, r)) 	  
	  fac.err.sd.true = diag(r) * fac.obs.noise.ratio
	  obs.err.sd.true = diag(n) 	  
	  sim.data.list <- replicate(n = n.mc,
		                         expr = simulateFactorModel(n=n, r=r, T=T+50, Lambda=Lambda.true, Phi=Phi.true, Y.start = Y.start, fac.err.sd=fac.err.sd.true, obs.err.sd=obs.err.sd.true),
                                 simplify = FALSE)	
      name.string = paste0("XY.r",r, "T",T, "n",n, ".FONR",fac.obs.noise.ratio)						 
	  assign(name.string, 
	         list(Lambda.true=Lambda.true, 
			      Phi.true=Phi.true,
				  fac.err.sd.true=fac.err.sd.true,
				  obs.err.sd.true=obs.err.sd.true,
				  sim.data.list=sim.data.list))	  
	  #str(get(name.string))	  
	  save(list=name.string, file=file.path(out.dir, paste0(name.string, ".RData")))
	  rm(list=c(name.string))
	}  ## end for( n in n.vec)
  }  ## end for (T in T.vec)
}

calculatePerformance <- function( T.vec=c(100, 300), 
                           n.vec=c(50, 100, 150, 200, 500, 700,1000), 
						   r=2, 
						   n.mc=100, 
						   fac.obs.noise.ratio = 1,
						   h = 1,
						   w = 1,
						   in.dir = directory("data"),
						   out.dir	= directory("results"),
						   seed = 20141221
						 ) {  
  
  #--- MC simulation --- 
  #sink(file.path(out.dir, paste0("MC.simulation.r",r, ".FONR",fac.obs.noise.ratio, ".log")))
  
  res.summary = NULL
  for(T in T.vec) {
    # T = 300
    for(n in n.vec) {
	  # n = 1000; T = 300
	  ptm.lss = NULL
      ptm.pca = NULL	   
       	
      PFF.lss <- rep(NA, n.mc)
      PFF.pca <- rep(NA, n.mc)
      FF.lss <- rep(NA, n.mc)
      FF.pca <- rep(NA, n.mc)
		 
	  Y.lss.space.dist <- rep(NA, n.mc)
	  Y.pca.space.dist <- rep(NA, n.mc)
	  factor.lss.space.closer <- rep(NA, n.mc)
		 
	  Lambda.oracle.space.dist <- rep(NA, n.mc)
	  Lambda.lss.space.dist <- rep(NA, n.mc)
	  Lambda.pca.space.dist <- rep(NA, n.mc)
	  Lambda.lss.space.closer <- rep(NA, n.mc)
		 
      obs.mse.lss <- rep(NA, n.mc)
      fac.mse.lss <- rep(NA, n.mc)
      obs.mse.pca <- rep(NA, n.mc) 
      fac.mse.pca <- rep(NA, n.mc) 
	  obs.mse.oracle <- rep(NA, n.mc)
      fac.mse.oracle <- rep(NA, n.mc) 		 
		 
      Lambda.lss <- array(NA, c(n, r, n.mc))
      Phi.lss <- array(NA, c(r, r, n.mc))
      Y.pred.lss <- array(NA, c(h, r, n.mc))
      Xh.pred.lss <- array(NA, c(h, n, n.mc))
  
      Xh.pred.err.Fnorm.lss <- rep(NA, n.mc)
      Xh.pred.err.Fnorm.pca <- rep(NA, n.mc)
         
	  Xh.pred.err.Fnorm.oracle <- rep(NA, n.mc)

      Xh.pred.mse.lss <- rep(NA, n.mc)
      Xh.pred.mse.pca <- rep(NA, n.mc)
      Xh.pred.mse.oracle <- rep(NA, n.mc)

	  # n = 150; T = 300	   
	  # Load simulation data
	  name.string = paste0("XY.r",r, "T",T, "n",n, ".FONR",fac.obs.noise.ratio)
	  load.name = load(file=file.path(in.dir, paste0(name.string, ".RData")))	 
      sim.dfm = get(load.name)
	  rm(list = c(load.name))
	  sim.data.total = length(sim.dfm$sim.data.list)
	  
	  Lambda.true = t(sim.dfm$Lambda.true)
	  Phi.true = sim.dfm$Phi.true
	  Phi1.true = Phi.true[[1]]
	  var.order = length(Phi.true)
	  
	  load.name = load(file=file.path(out.dir, paste0("res.",name.string,".w",w, ".RData")))	 
      sim.result.list = get(load.name)
	  rm(list = c(load.name))
	  
	  # i = 1
	  set.seed(seed)
	  for(i in 1:sim.data.total) { #
		 dfm.obs = sim.dfm$sim.data.list[[i]]$X
		 dfm.fac = sim.dfm$sim.data.list[[i]]$Y
         X.train = dfm.obs[1:T,]
         X.test = dfm.obs[(T+1):nrow(dfm.obs),]
		 Y.true = dfm.fac[1:T,]
		 Y.test = dfm.fac[(T+1):nrow(dfm.fac),]
		 
		 X.train = scale(X.train, center = TRUE, scale = TRUE)		 
		 
		 cat("Estimation for r = ",r, ", T = ",T, ", n = ",n , ", MC = ", i, "\n"); flush.console()
		 
		 # LS Stiefel 
		 #proc.start = proc.time()
		 #lss.res = ls.stiefel(X.train, r=r, w=w, i.mc=i)            
         #proc.time() - proc.start			 
         #ptm.lss = rbind(ptm.lss, (proc.time() - proc.start)["elapsed"])
	     #Y.lss = lss.res$Y[,,lss.res$k.opt]
			
		 # PCA 
		 #proc.start = proc.time()
         #pca.res = prcomp(X.train, retx=TRUE, center=TRUE, scale=TRUE)
         #Y.pca = pca.res$x[,1:r]; #dim(Y.pca)		 
         #ptm.pca = rbind(ptm.pca, (proc.time() - proc.start)["elapsed"])
		 ## save simulation results
		 #sim.result.list[[i]] <- list(lss.res=lss.res, Y.pca=Y.pca)
		 
		 lss.res = sim.result.list[[i]]$lss.res
		 Y.lss = lss.res$Y[,,lss.res$k.opt]
         Y.pca = sim.result.list[[i]]$Y.pca		 
		 
		 ## Calculate performance criterions
		 #--- 0. Factor Space Distance 
		 Y.lss.space.dist[i] <- spaceDist(orthonormal_basis(Y.true), orthonormal_basis(Y.lss))
		 Y.pca.space.dist[i] <- spaceDist(orthonormal_basis(Y.true), orthonormal_basis(Y.pca))
		 factor.lss.space.closer[i] <- (Y.pca.space.dist[i]-Y.lss.space.dist[i])/Y.pca.space.dist[i]
		 		 
		 #--- 1. Trace R2 of the multivariate regression of estimated factors onto true factors 
         #Pf = Y.true %*% ginv(t(Y.true) %*% Y.true) %*% t(Y.true)
         Pf = Y.true %*% ginv(t(Y.true) %*% Y.true) %*% t(Y.true)
         PFF.lss[i] = matTrace(t(Y.lss) %*% Pf %*% Y.lss)
         FF.lss[i] = matTrace(t(Y.lss) %*% Y.lss)
         PFF.pca[i] = matTrace(t(Y.pca) %*% Pf %*% Y.pca)
         FF.pca[i] = matTrace(t(Y.pca) %*% Y.pca)
  
         #--- 2. How close the forecast based on the estimated factors is to the infeasible forecast based on the true factors
         dfm.lss = dfm(X.train, Y.lss, h=1)
         dfm.pca = dfm(X.train, Y.pca, h=1)
         dfm.oracle = dfm(X.train, Y.true, h=1)
  
         obs.mse.lss[i] = dfm.lss$obs.mse
         fac.mse.lss[i] = dfm.lss$fac.mse		 
         Lambda.lss = dfm.lss$Lambda.est
         Phi.lss = dfm.lss$Phi.est
         Y.pred.lss = dfm.lss$fac.pred
         Xh.pred.lss = dfm.lss$obs.pred
  
         obs.mse.pca[i] = dfm.pca$obs.mse 
         fac.mse.pca[i] = dfm.pca$fac.mse 		 
         Lambda.pca = dfm.pca$Lambda.est
         Phi.pca = dfm.pca$Phi.est
         Y.pred.pca = dfm.pca$fac.pred
         Xh.pred.pca = dfm.pca$obs.pred
  
         obs.mse.oracle[i] = dfm.oracle$obs.mse
         fac.mse.oracle[i] = dfm.oracle$fac.mse 
         Lambda.oracle = dfm.oracle$Lambda.est
         Phi.oracle = dfm.oracle$Phi.est
         Y.pred.oracle = dfm.oracle$fac.pred
         Xh.pred.oracle = dfm.oracle$obs.pred      
		 
		 Lambda.oracle.space.dist[i] <- spaceDist(orthonormal_basis(t(Lambda.true)), orthonormal_basis(t(Lambda.oracle)))
		 Lambda.lss.space.dist[i] <- spaceDist(orthonormal_basis(t(Lambda.true)), orthonormal_basis(t(Lambda.lss)))
		 Lambda.pca.space.dist[i] <- spaceDist(orthonormal_basis(t(Lambda.true)), orthonormal_basis(t(Lambda.pca)))
		 Lambda.lss.space.closer[i] <- (Lambda.pca.space.dist[i]-Lambda.lss.space.dist[i])/Lambda.pca.space.dist[i]
  
         #--- Estimate X using Lambda ---
         Xh.true = X.test[1:h,]

         Xh.pred.err.lss = Xh.pred.lss - Xh.true
         Xh.pred.mse.lss[i] = mean(Xh.pred.err.lss^2)  
  
         Xh.pred.err.pca = Xh.pred.pca - Xh.true
         Xh.pred.mse.pca[i] = mean(Xh.pred.err.pca^2)
  
         Xh.pred.err.oracle = Xh.pred.oracle - Xh.true         
         Xh.pred.mse.oracle[i] = mean(Xh.pred.err.oracle^2)		 
		 
	  } # end MC
		 
	  #assign(paste0("res.",name.string,".w",w), sim.result.list)	
      #save(list=paste0("res.",name.string,".w",w), file=file.path(out.dir, paste0("res.",name.string,".w",w, ".RData")))	
      #rm(list=c(paste0("res.",name.string,".w",w)))	  
		 
	  R2FF.lss <- mean(na.omit(PFF.lss)) / mean(na.omit(FF.lss)); R2FF.lss
      R2FF.pca <- mean(na.omit(PFF.pca)) / mean(na.omit(FF.pca)); R2FF.pca	  
	  
	  pred.mse.lss = mean(na.omit(Xh.pred.mse.lss)); pred.mse.lss
      pred.mse.pca <- mean(na.omit(Xh.pred.mse.pca)); pred.mse.pca
      pred.mse.oracle <- mean(na.omit(Xh.pred.mse.oracle)); pred.mse.oracle
	  
	  # T.vec <- c(100, 200); n.vec <- c(10, 50, 100, 150, 200, 500, 1000)
	  res.summary = rbind( res.summary, 
	                   data.frame(T = c(T, T), 
					              N = c(n, n), 
								  Describ = c("Mean", "Sd"),
	                              Factor.R2FF.LSS = c(R2FF.lss,NA), 
					              Factor.R2FF.PCA = c(R2FF.pca, NA), 
					              Factor.Space.Dist.LSS = c(mean(Y.lss.space.dist), sd(Y.lss.space.dist)), 
					              Fact.Space.Dist.PCA = c(mean(Y.pca.space.dist), sd(Y.pca.space.dist)), 
								  LSS.Space.Closer = c(mean(factor.lss.space.closer), sd(factor.lss.space.closer)),					 
					              Loading.Space.Dist.LSS = c(mean(Lambda.lss.space.dist), sd(Lambda.lss.space.dist)),
						          Loading.Space.Dist.PCA = c(mean(Lambda.pca.space.dist), sd(Lambda.pca.space.dist)), 
								  Loading.lss.Space.Closer = c(mean(Lambda.lss.space.closer), sd(Lambda.lss.space.closer)),
					              PredX.MSE.Oracle = c(mean(Xh.pred.mse.oracle), sd(Xh.pred.mse.oracle)), 
					              PredX.MSE.LSS = c(mean(Xh.pred.mse.lss), sd(Xh.pred.mse.lss)), 
					              PredX.MSE.PCA = c(mean(Xh.pred.mse.pca), sd(Xh.pred.mse.pca)), 
						          ptm.lss = c(mean(ptm.lss), sd(ptm.lss)),
						          ptm.pca = c(mean(ptm.pca), sd(ptm.pca)) ) )		 
	  } # end n
	} # end T
	
	write.table(res.summary, file=file.path(out.dir, paste0("summary.all.","r",r,".FONR",fac.obs.noise.ratio,".w",w, ".csv")), sep="|", row.names=FALSE, append=TRUE)
  
    write.table(res.summary[which(res.summary$Describ == "Mean"), ], sep="|", file=file.path(out.dir, paste0("summary.all.","r",r,".FONR",fac.obs.noise.ratio,".w",w, ".mean.csv")), row.names=FALSE, append=TRUE)
	   
    write.table(res.summary[which(res.summary$Describ == "Sd"), ], sep="|", file=file.path(out.dir, paste0("summary.all.","r",r,".FONR",fac.obs.noise.ratio,".w",w, ".sd.csv")), row.names=FALSE, append=TRUE)
		
	sink.reset()
}


main.r2 <- function() {
  
  h <- 1  #- forecasting length 
  T.vec <- c(100, 300) 
  n.vec <- c(50, 100, 150, 200, 500, 700, 1000)
  #n.vec <- c(200, 500, 1000)
  r <- 2
  n.mc <- 100
  
  #fac.obs.noise.ratio = 0.1   # with w = 10 in ls.stiefel()
  fac.obs.noise.ratio = 1
  
  generateMCData( T.vec=T.vec, n.vec=n.vec, r=r, n.mc=n.mc, fac.obs.noise.ratio = fac.obs.noise.ratio, Lambda.true = matrix( runif(n*r, min=0, max=1) , nrow=n, ncol=r), Phi.true = list(matrix(c(0.7, 0, 0, -0.2), r, r)), Y.start = c(5, 1), out.dir	= directory("data") )
  
  DFMSimulation( T.vec=T.vec, n.vec=n.vec, r=r, n.mc=n.mc, fac.obs.noise.ratio = fac.obs.noise.ratio, h = h, w=1/fac.obs.noise.ratio, in.dir = directory("data"), out.dir	= directory("results"), seed = 20141221 ) 
  
  calculatePerformance( T.vec=T.vec, n.vec=n.vec, r=r, n.mc=n.mc, fac.obs.noise.ratio = fac.obs.noise.ratio, h = h, w=1/fac.obs.noise.ratio, in.dir = directory("data"), out.dir	= directory("results"), seed = 20141221 ) 

  
}

main.r3 <- function() {
  
  h <- 1  #- forecasting length 
  T.vec <- c(100, 300) 
  n.vec <- c(50, 100, 150, 200, 500, 1000)
  #n.vec <- c(200, 500, 1000)
  r <- 3
  n.mc <- 100  
  
  # fac.obs.noise.ratio = 1
  fac.obs.noise.ratio = 0.1   # with w = 10 in ls.stiefel()
  
  generateMCData( T.vec=T.vec, n.vec=n.vec, r=r, n.mc=n.mc, fac.obs.noise.ratio = fac.obs.noise.ratio, Phi.true = list(matrix(c(0.7, 0, 0, 0, -0.4, 0, 0, 0, 0.1), r, r)), Y.start = c(5, 1, 3) )
  
  DFMSimulation( T.vec=T.vec, n.vec=n.vec, r=r, n.mc=n.mc, fac.obs.noise.ratio = fac.obs.noise.ratio, h = h, w=1/fac.obs.noise.ratio) 

  
}

main.r5 <- function() {
  
  h <- 1  #- forecasting length 
  T.vec <- c(100, 300) 
  n.vec <- c(50, 100, 150, 200, 500, 1000)
  #n.vec <- c(200, 500, 1000)
  r <- 5
  n.mc <- 100  
  
  fac.obs.noise.ratio = 1
  # fac.obs.noise.ratio = 0.1   # with w = 10 in ls.stiefel()
  
  generateMCData( T.vec=T.vec, n.vec=n.vec, r=r, n.mc=n.mc, fac.obs.noise.ratio = fac.obs.noise.ratio, Phi.true = list(matrix(c(0.7, 0, 0, 0, -0.4, 0, 0, 0, 0.1), r, r)), Y.start = c(5, 1, 3) )
  
  DFMSimulation( T.vec=T.vec, n.vec=n.vec, r=r, n.mc=n.mc, fac.obs.noise.ratio = fac.obs.noise.ratio, h = h, w=1/fac.obs.noise.ratio) 

  
}


	 