###############################################################################################
# Tunning Random Forests for Causal Inference Under Cluster-level Unmeasured Confounding
# Youmi Suk and Hyunseung Kang 

# :::: Simulation codes ::::

# :: load packages and DGP codes
source("DGP_clusterOVB.R") # R codes for data generating models 
library(grf)
library(drgee)
library(lme4)
library(dplyr)

li.clusterATE <- function(ps, Y, Z, id) {
  dataf <- data.frame(Y, Z, id)
  dataf$temp.wt <- with(dataf, Z/ps + (1-Z)/(1-ps))
  temp.dat <- data.frame(clusterATE= dataf %>% group_by(id) %>% summarize(clusterATE=lm(Y ~ Z, weights=temp.wt)$coef[2], .groups = 'drop') %>% pull(clusterATE),
                         wt=tapply(dataf$temp.wt, dataf$id, sum))
  
  return(list(ATE=sum(temp.dat$clusterATE * temp.dat$wt, na.rm = TRUE) / sum(temp.dat$wt[!is.na(temp.dat$clusterATE)]), clusterATE=temp.dat[,1])) # Li's method
}

IPWfun <- function(ps, Z) Z/ps + (1-Z)/(1-ps)

# :: set simulation parameters
smpl.size = "200.20" # sample size (the number of clusters, the mean cluster size)
o.val=2  # the main effect of U_j in the outcome model; (2, 4) for linear outcome model and (0.3, 0.9) for binary outcome model
m.val=0  # the cross-level interaction with U_j in the outcome model; (0, 2, 4) for linear outcome model and (0, 0.3, 0.9) for binary outcome model

reps <- 1000 # the number of replications

# save simulation results
simu.rlst <- data.frame(matrix(NA, nrow=reps, ncol=21))

for (i in 1:reps) {
  
# Please choose DGP
  dat <- twolevel.normalU(Smpl.size = smpl.size, o.val=o.val, m.val=m.val) # Design 1
#  dat <- twolevel.uniformU(Smpl.size = smpl.size, o.val=o.val, m.val=m.val) # Designs 2 and 3
#  dat <- twolevel.binary(Smpl.size = smpl.size, o.val=o.val, m.val=m.val) # Design 4
#  dat <- twolevel.experror(Smpl.size = smpl.size, o.val=o.val, m.val=m.val) # Design 5
  
  # group deviation
  dat$X1grpcent <- dat$X1 - ave(dat$X1, dat$id)
  dat$X2grpcent <- dat$X2 - ave(dat$X2, dat$id)
  dat$X3grpcent <- dat$X3 - ave(dat$X3, dat$id)
  dat$Ygrpcent <- dat$Y -ave(dat$Y, dat$id)
  dat$Zgrpcent <- dat$Z -ave(dat$Z, dat$id)
  dat$X1sqgrpcent <- dat$X1^2 - ave(dat$X1^2, dat$id)
  dat$X2W1grpcent <- dat$X2*dat$W1 - ave(dat$X2*dat$W1, dat$id)
  dat$X1X3indgrpcent <- dat$X1*(dat$X3 < 0.3) - ave(dat$X1*( dat$X3 < 0.3), dat$id)
  
  # group mean
  dat$meanX1 <- ave(dat$X1, dat$id)
  dat$meanX2 <- ave(dat$X2, dat$id)
  dat$meanX3 <- ave(dat$X3, dat$id)
  dat$meanZ <- ave(dat$Z, dat$id)
  
  # ::: comparison: DRCGEE
  both.cor <- drgee(oformula=formula(Y ~ X1 + X2 + X3 + I(X1^2) + I(X2*W1) + I(X1*( X3 < 0.3))), eformula=(Z ~ X1 + X2 + X3 + I(X1^2) + I(X2*W1) + I(X1*( X3 < 0.3))), 
                    iaformula = formula(~ X3 + W1), olink="identity", elink="identity", 
                    data=dat, estimation.method= "dr", clusterid = dat$id, cond=T)
  drgee.est <- mean(as.matrix(cbind(1, dat[, c("X3", "W1")]))%*%summary(both.cor)$coef[,1]) 
  
  # when the propensity score model is misspecified.
  mis.sel <- drgee(oformula=formula(Y~ X1 + X2 + X3 + I(X1^2) + I(X2*W1) + I(X1*( X3 < 0.3))), eformula=(Z~ X1 + X2  + X3), 
                   iaformula = formula(~ X3 + W1), olink="identity", elink="identity", 
                   data=dat, estimation.method= "dr", clusterid = dat$id, cond=T)
  
  drgee.mis.sel.est <- mean(as.matrix(cbind(1, dat[, c("X3", "W1")]))%*%summary(mis.sel)$coef[,1])
  
  # when both the propensity score model and the outcome model are misspecified.
  strong.man <- drgee(oformula=formula(Y~ X1 + X2 + X3), eformula=(Z~ X1 + X2  + X3), 
                      iaformula = formula(~ X3 + W1), olink="identity", elink="identity", 
                      data=dat, estimation.method= "dr", clusterid = dat$id, cond=T)
  
  drgee.mis.est <- mean(as.matrix(cbind(1, dat[, c("X3", "W1")]))%*%summary(strong.man)$coef[,1])
  
  # ::: comparison: IPW estiamtors 
  sel.re <- glmer(Z~ X1 + X2 + X3 + W1 + I(X2*W1) + I(X1^2) + I(X1*( X3 < 0.3)) +  (1|id), data = dat, family = binomial) # random effects logistic 
  ps.est.re <- predict(sel.re, type = 'response')   # propensity score  
  
  sel.re.mis <- glmer(Z~ X1 + X2 + X3 + W1 + (1|id), data = dat, family = binomial)
  ps.est.re.mis <- predict(sel.re.mis, type = 'response')   # mis-specified propensity score  
  
  sel.fe <- glm(Z ~ X1 + X2  + X3 + I(X2*W1) + I(X1^2) + I(X1*( X3 < 0.3)) + factor(id), data = dat, family = binomial) # fixed effects logistic 
  ps.est.fe <- predict(sel.fe, type = 'response')   # propensity score  
  
  sel.fe.mis <- glm(Z ~ X1 + X2  + X3 + factor(id), data = dat, family = binomial)
  ps.est.fe.mis <- predict(sel.fe.mis, type = 'response')   # mis-specified propensity score   
  
  # marignal IPW estimator
  marginal.rlst <- rbind(marginal.re=summary(lm(Y~Z, data=dat, weights=IPWfun(ps.est.re, dat$Z)))$coef[2,],
                     marginal.re.mis=summary(lm(Y~Z, data=dat, weights=IPWfun(ps.est.re.mis, dat$Z)))$coef[2,1],
                     marginal.fe=summary(lm(Y~Z, data=dat, weights=IPWfun(ps.est.fe, dat$Z)))$coef[2,],
                     marginal.fe.mis=summary(lm(Y~Z, data=dat, weights=IPWfun(ps.est.fe.mis, dat$Z)))$coef[2,])[,1]
  
  # clustered IPW estimator
  cluster.rlst <- c(cluster.re=li.clusterATE(ps.est.re,  dat$Y, dat$Z, dat$id)$ATE,
                    cluster.re.mis=li.clusterATE(ps.est.re.mis,  dat$Y, dat$Z, dat$id)$ATE,
                    cluster.fe=li.clusterATE(ps.est.fe,  dat$Y, dat$Z, dat$id)$ATE,
                    cluster.fe.mis=li.clusterATE(ps.est.fe.mis,  dat$Y, dat$Z, dat$id)$ATE)
  
  
  # ::: our modifications: Causal Forests (CF)
  # the default CF
  out.cf <- causal_forest(X=dat[,c("X1","X2","X3", "W1")], Y=dat$Y, W=dat$Z)
  pred.cf = predict(out.cf, type="vector", estimate.variance = TRUE)
  
  # CF+RePS
  out.cf.ps.re <- causal_forest(X=dat[,c("X1","X2","X3", "W1")], Y=dat$Y, W=dat$Z, W.hat=ps.est.re) # ps.est from random-effects logistic regression
  pred.cf.ps.re = predict(out.cf.ps.re, type="vector", estimate.variance = TRUE)
  
  # CF+RePS with mis-specified propensity scores
  out.cf.ps.re.mis <- causal_forest(X=dat[,c("X1","X2","X3", "W1")], Y=dat$Y, W=dat$Z, W.hat=ps.est.re.mis) # ps.est from random-effects logistic regression
  pred.cf.ps.re.mis = predict(out.cf.ps.re.mis, type="vector", estimate.variance = TRUE)
  
  # CF+FePS
  out.cf.ps.fe <- causal_forest(X=dat[,c("X1","X2","X3", "W1")], Y=dat$Y, W=dat$Z, W.hat=ps.est.fe) # ps.est from fixed-effects logistic regression
  pred.cf.ps.fe = predict(out.cf.ps.fe, type="vector", estimate.variance = TRUE)
  
  # CF+FePS with mis-specified propensity scores
  out.cf.ps.fe.mis <- causal_forest(X=dat[,c("X1","X2","X3", "W1")], Y=dat$Y, W=dat$Z, W.hat=ps.est.fe.mis) # ps.est from fixed-effects logistic regression
  pred.cf.ps.fe.mis = predict(out.cf.ps.fe.mis, type="vector", estimate.variance = TRUE)
  
  # CF+CLusterID
  id.dummes <- model.matrix( ~ id + 0, dat)
  out.cf.id.dummies <- causal_forest(X=cbind(dat[,c("X1","X2","X3", "W1")], id.dummes), Y=dat$Y, W=dat$Z)
  pred.cf.id.dummies = predict(out.cf.id.dummies, type="vector", estimate.variance = TRUE)
  
  # CF+Demean
  out.cf.demean <- causal_forest(X=dat[, c("X1grpcent", "X2grpcent", "X3grpcent", "meanX1", "meanX2", "meanX3", "W1")], Y=dat$Ygrpcent, W=dat$Zgrpcent)
  pred.cf.demean = predict(out.cf.demean, type="vector", estimate.variance = TRUE)
  
  # CF+Demean+PS 
  DD.pred <- predict(lm(Zgrpcent ~ X1grpcent + X2grpcent + X3grpcent + X2W1grpcent + X1sqgrpcent + X1X3indgrpcent , data=dat))
  out.cf.demean_psDD <- causal_forest(X=dat[, c("X1grpcent", "X2grpcent", "X3grpcent", "meanX1", "meanX2", "meanX3", "W1")], Y=dat$Ygrpcent, W=dat$Zgrpcent, W.hat=DD.pred)
  pred.cf.demean_psDD = predict(out.cf.demean_psDD, type="vector", estimate.variance = TRUE)
  
  # CF+Demean+PS with mis-specified propensity scores
  DD.pred.mis <- predict(lm(Zgrpcent ~ X1grpcent + X2grpcent + X3grpcent, data=dat))
  out.cf.demean_psDD.mis <- causal_forest(X=dat[, c("X1grpcent", "X2grpcent", "X3grpcent", "meanX1", "meanX2", "meanX3", "W1")], Y=dat$Ygrpcent, W=dat$Zgrpcent, W.hat=DD.pred.mis)
  pred.cf.demean_psDD.mis = predict(out.cf.demean_psDD.mis, type="vector", estimate.variance = TRUE)
  
  rlst.blp <- c(CF=best_linear_projection(out.cf, A=NULL)[,1],
                CF.RePS=best_linear_projection(out.cf.ps.re, A=NULL)[,1],
                CF.RePS.mis=best_linear_projection(out.cf.ps.re.mis, A=NULL)[,1],                    
                CF.FePS=best_linear_projection(out.cf.ps.fe, A=NULL)[,1],
                CF.FePS.mis=best_linear_projection(out.cf.ps.fe.mis, A=NULL)[,1],                    
                CF.ClusterID=best_linear_projection(out.cf.id.dummies, A=NULL)[,1],
                CF.Demean=best_linear_projection(out.cf.demean, A=NULL)[,1],
                CF.Demean.PS=best_linear_projection(out.cf.demean_psDD, A=NULL)[,1],
                CF.Demean.PS.mis=best_linear_projection(out.cf.demean_psDD.mis, A=NULL)[,1])

  simu.rlst[i, ] <- c(sampleATE=mean(dat$Y1-dat$Y0), drgee.est=drgee.est, drgee.mis.sel.est=drgee.mis.sel.est, drgee.mis.est=drgee.mis.est, marginal.rlst, cluster.rlst, rlst.blp)
  
}
