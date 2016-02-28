
dfm <- function(X, Y, h=1) {
  ## estimate measurement equation
  obs.model <- lm(X ~ Y - 1)
  obs.mse <- mean(resid(obs.model)^2)
  Lambda.est <- coef(obs.model)
  
  ## estimate factor equation
  fac.var1 <- VAR(Y, p=1, type="none")
  fac.mse <- mean(resid(fac.var1)^2) 
  Phi.est = do.call(rbind, lapply(coef(fac.var1), function(x){t(subset(x, select="Estimate"))}))
  
  ## predict
  fac.pred <- do.call(cbind, lapply(predict(fac.var1, n.ahead=h)$fcst, function(x){subset(x, select="fcst", drop=F)}))
  
  #obs.pred <- predict(obs.model, as.data.frame(fac.pred))
  obs.pred <- fac.pred %*% Lambda.est
  
  return(list(obs.model=obs.model, obs.mse=obs.mse, Lambda.est=Lambda.est, fac.var1=fac.var1, fac.mse=fac.mse, Phi.est=Phi.est, fac.pred=fac.pred, obs.pred=obs.pred, h=h))
}