source(paste(getwd(),"/simXandY.r", sep=""))

#n.mc = 100
#n <- 1000; r <- 2; T <- 150

TEST <- FALSE

genData <- function(n, r, T, i.mc){  
    assign(paste("XY.n",n, "r",r, "T",T, ".",mc, sep=""), simXandY(n=n, r=r, T=T))
    save(list=paste("XY.n",n, "r",r, "T",T, ".",mc, sep=""), file=paste(getwd(),"/data/XY.n",n, "r",r, "T",T, ".", i.mc, ".RData",sep="")) 
}