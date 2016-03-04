# X = X.train; r = r; w = 1; K = 1000; H = 5; i.mc = 1;
ls.stiefel = function( X, r=2, w=1, K=1000, H=5, i.mc, out.dir ,out.file ) {
  #w = 1
  #K = 1000  # maximum number of iteration
  n = ncol(X)
  T = nrow(X)
  Y = array(NA, c(T, r, K))
  L = rep(NA, K)
  G = array(NA, c(T, r, K))
  DF = array(NA, c(T, r, K))
  toly = rep(NA, K)
  toll = rep(NA, K)
  
  ###
  #H = 5   # window length
  eps <- epsy <- epsl <- 1e-5

  ### 
  #set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
  #set.seed(20141221)
  
  Y[,,1] = genOrthogonal(T)[, 1:r] 
  tmp = loss_func(X, Y[,,1], w) #--- Y = Y[,,1]
  L[1] = tmp$L
  G[,,1] = tmp$G
  DF[,,1] = G[,,1] - Y[,,1] %*% t(G[,,1]) %*% Y[,,1]
  flag.conv = 0
  conv.gradient.norm = FALSE
  conv.relative.change.yl = FALSE
  conv.relative.change.ylmean = FALSE
  proc.start = proc.time()
  k = 1
  sink(file.path(out.dir, out.file), append=TRUE)
  cat("@@@@@ Monte Carlo Step ", i.mc, " @@@@@\n"); flush.console()
  while(k < K) {  
    U = cbind(G[,,k], Y[,,k])
    V = cbind(Y[,,k], -G[,,k])   
    
    if(TRUE) {
	  cat("# Searching Y Step ", k, "\n"); flush.console()
    }
    ###
  
    #Y.saved = Y; Y = Y[,,k]; G.saved = G; G = G[,,k]
    Y[,,k+1] = curvilinear_search(X, Y[,,k], G[,,k], U, V, w=1, k, rho1=1e-4, rho2=0.9, tau=0.01) 
	# Y = Y.saved; G = G.saved
	# Y[,,k+1] = curve.val.at.tau.opt
    
    tmp = loss_func(X, Y[,,k+1], w)  
    L[k+1] = tmp$L
    G[,,k+1] = tmp$G  
  
    DF[,,k+1] = G[,,k+1] - Y[,,k+1] %*% t(G[,,k+1]) %*% Y[,,k+1]
  
    toly[k] = norm( Y[,,k]-Y[,,k+1], "F" ) / sqrt(T)
    toll[k] = (L[k] - L[k+1]) / (abs(L[k])+1)
  
    conv.gradient.norm = (norm(DF[,,k],"M") <= eps)
    conv.relative.change.yl = (toly[k] <= epsy & toll[k] <= epsl)
    conv.relative.change.ylmean = (mean(toly[I(k-min(k,H)+1):k]) <= epsy & mean(toll[I(k-min(k,H)+1):k]) <= epsl)
    ###

    ###
    # | conv.relative.change1 | conv.relative.change2
    if( conv.gradient.norm | conv.relative.change.yl | conv.relative.change.ylmean ) {
      flag.conv = 1
	  break	# break before k+1, thus k is the optimal
    }
    k = k + 1
  } ## end while(TRUE)
  
  
  if(TEST) {
    cat("-conv.gradient.norm:", conv.gradient.norm, "\n"); flush.console()
    cat("-conv.relative.change.yl:", conv.relative.change.yl, "\n"); flush.console()
    cat("-conv.relative.change.ylmean:", conv.relative.change.ylmean, "\n"); flush.console()
	}
    sink()

  if(conv.gradient.norm | conv.relative.change.yl | conv.relative.change.ylmean) {
    k.opt = k
  } else {k.opt = k - 1}

  proc.end = proc.time()

  proc.duration = proc.end - proc.start
  proc.duration

  k.opt

  #orthon.res.01 = list(k.opt=k.opt, Y=Y, L=L, G=G, DF=DF, toly=toly, toll=toll)
  #save(orthon.res.01, file=paste(getwd(),"/test/orthon.res.01.RData", sep=""))

  #load(file=paste(getwd(),"/XY.n10.r2.T100/res.01.RData",sep=""))

  #feasi = norm( t(Y[,,k]) %*% Y[,,k] - diag(ncol(Y[,,k])), "F"); feasi
  #norm(DF[,,k],"M")
  #L[k]

  #install.packages("vars")
  return(list(k.opt=k.opt, Y=Y, L=L, G=G, DF=DF, toly=toly, toll=toll))
}

