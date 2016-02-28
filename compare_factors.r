# compare ture Y's and estimated Y's
par(mfrow=c(2,1))
ts.plot(Y.true, lty=1:2, ylab="y1, y2")
legend(locator(1), legend=c("y1","y2"), lty=1:2)
title("True Y's")
ts.plot(Y[,,k], lty=1:2, ylab="y1, y2")
title("Estimated Y's")
#ts.plot(cbind(Y.true[,1],Y[,1,k]), col=1:2, ylab="")
legend(locator(1), legend=colnames(transY), lty=1:2)


y.model = lm(Y.true ~ Y[,,k] - 1)
y.model = lm(Y.true ~ Y.lss - 1)
summary(y.model)

y1.model = lm(Y.true[,1] ~ Y[,2,k] - 1)
summary(y1.model)

y2.model = lm(Y.true[,2] ~ Y[,2,k] - 1)
summary(y2.model)




