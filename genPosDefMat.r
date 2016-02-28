#source("genPosDefMat.r")
#--- ref: genPositiveDefMat in Package ¡®clusterGeneration¡¯
#--- getCluster.R , line 2105
#--- https://cran.r-project.org/web/packages/clusterGeneration/clusterGeneration.pdf

#--- generate orthogonal matrix
#--- dim : dimenstion
genOrthogonal <- function(dim) {
  Q <- MOrthogonal(runif(dim))
  return(Q)
}

#--- Construct an orthogonal matrix whose first few columns are standardized 'M', where columns of 'M' are orthogonal, Here "standardized 'M'" means each its columns has length 1
MOrthogonal <- function(M) {
  # can set the parameter "tol" of "qr" to decide how small value shoule be 0
  # M <- runif(dim)
  tmp <- qr(M)
  Q <- qr.Q(tmp, complete=TRUE)
  if(is.vector(M)) { 
    if( Q[1]*M[1]<0 ) Q <- -Q 
  } else { 
    if( Q[1,1]*M[1,1]<0 ) Q <- -Q 
  }
  return(Q)
}

genPosDefMat <- function( dim, lambdaLow=1, ratioLambda=10 ) {
  low <- lambdaLow
  upp <- lambdaLow * ratioLambda
  D <- matrix(0, dim, dim) 
  eigenvalues <- runif(dim, min=low, max=upp)
  D <- diag(eigenvalues) 
  if(dim > 1) {
    Q <- genOrthogonal(dim)
	Sigma <- Q %*% D %*% t(Q)
  }  # ... dim = 1 ?????
  return(list(eigenvalues=eigenvalues, Sigma=Sigma))
}