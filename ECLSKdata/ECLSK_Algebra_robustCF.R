####################################################################################################
# Tuning Random Forests for Causal Inference Under Cluster-Level Unmeasured Confounding
# : The ECLS-K data analysis 
# by Youmi Suk & Hyunseung Kang

# This is our complete data. 
# The original ECLSK 1998-99 dataset is available at https://nces.ed.gov/ecls/dataproducts.asp. 
# For more information, see Tourangeau et al. (2009) from https://nces.ed.gov/pubs2009/2009003.pdf.

# The variables in the dataset include:

# :: ID
# S7_ID    : school ID

# :: treatment
# C7DESMTH : whether students took the algebra or a higher-level math course (= 1) or not (= 0)

# :: outcome
# C7R4MSCL : math IRT scale score

# :: student-level covariates
# C6R4MSCL : prior math IRT scale score
# P6EXPECT : expected degree of child
# GENDER : male or female
# RACE : race/ethnicity
# WKSESL : socio-economic status measure
# WKPOV_R : poverty level
# WKMOMED : mother's education level
# P7HFAMIL : family type

# :: school-level covariates
# S7PUPRI : public or private school
# R7REGION : census region - WEST, MIDWEST, NORTHEAST, SOUTH
# R7URBAN : school location - urban, suburb, small town
####################################################################################################

# load packages/sources
library(grf)
library(lme4)
library(drgee)
library(dplyr)

demean_mean_covs <- function(X, ID) {
  IDnum <- as.numeric(factor(ID))
  IDfactor <- as.factor(ID)
  X_demeaned <- apply(X, 2, function(x) x - ave(x, IDnum))
  
  # Data structures for h2o
  lvl2var <-  colSums(apply(as.matrix(X), 2, function(x) tapply(x, IDnum, function(x) length(unique(x)))) == 1) == length(levels(IDfactor))
  X_lvl1 <- as.matrix(X)[, !lvl2var]
  X_demeaned_lvl1 <-  as.matrix(X_demeaned)[, !lvl2var]
  X_demeaned_lvl1_num <- as.matrix(as.matrix(X_demeaned_lvl1)[, apply(as.matrix(X_lvl1), 2, function(x) sum(!unique(x) %in% c(0, 1))) > 0])
  
  X_lvl1_mean <- apply(as.matrix(X)[, !lvl2var], 2, function(x) ave(x, IDnum))
  W_lvl2 <- X[, lvl2var]
  
  Data.Covs = data.frame(X_demeaned_lvl1, X_lvl1_mean, W_lvl2)
  
  return(Data.Covs)
  
}

IPWfun <- function(ps, Z) Z/ps + (1-Z)/(1-ps)

li.clusterATE.SE <- function(ps, Y, Z, id) {
  dataf <- data.frame(Y, Z, id)
  dataf$temp.wt <- with(dataf, Z/ps + (1-Z)/(1-ps))
  temp.dat <- data.frame(clusterATE= dataf %>% group_by(id) %>% summarize(clusterATE=lm(Y~Z, weights=temp.wt)$coef[2], .groups = 'drop') %>% pull(clusterATE),
                         clusterSE= dataf %>% group_by(id) %>% summarize(clusterSE=summary(lm(Y~Z, weights=temp.wt))$coef[2,2], .groups = 'drop') %>% pull(clusterSE),
                         wt=tapply(dataf$temp.wt, dataf$id, sum))
  temp.ate = sum(temp.dat$clusterATE * temp.dat$wt, na.rm = TRUE) / sum(temp.dat$wt[!is.na(temp.dat$clusterATE)])
  temp.se = sqrt(sum(temp.dat$wt^2*temp.dat$clusterSE^2, na.rm = TRUE)) / sum(temp.dat$wt, na.rm = TRUE)
  
  return(list(ATE=c(temp.ate, temp.se), clusterATE=temp.dat[,1])) # Li's method
}


# :: load data
dat <- read.csv("ECLSK_Algebra_complete.csv") 

dat$WKSESL_s <- as.numeric(scale(dat$WKSESL, center=T, scale=T))
dat$C6R4MSCL_s <- as.numeric(scale(dat$C6R4MSCL, center=T, scale=T))

covs.dum_c <- model.matrix(~ C6R4MSCL_s + WKSESL_s + P6EXPECT + GENDER +  RACE + WKPOV_R + WKMOMED + P7HFAMIL+ S7PUPRI + R7REGION + R7URBAN + 0,  data =dat) # create dummy variables for categorical confounders

covs_demean_c <- demean_mean_covs(X=covs.dum_c, ID=dat$S7_ID) # create cluster-demeaned and cluster-constant components of confounders

# :: implement the DRCGEE estimator
drgee_ATE <- summary(drgee(oformula=formula(C7R4MSCL ~ C6R4MSCL_s + P6EXPECT + GENDER +  RACE + WKSESL_s + WKPOV_R + WKMOMED + P7HFAMIL), 
                           eformula=formula(C7DESMTH ~  C6R4MSCL_s + I(C6R4MSCL_s^2) +  P6EXPECT  +  C6R4MSCL_s*GENDER + C6R4MSCL_s*RACE + C6R4MSCL_s*WKMOMED +  P7HFAMIL),  
                           iaformula=(~  C6R4MSCL_s), olink="identity", elink="identity", data=dat, estimation.method="dr", clusterid="S7_ID", cond=T))$coef[1, 1:2]

# :: implement the IPW estimators
sel.re <- glmer(C7DESMTH ~ C6R4MSCL_s +  I(C6R4MSCL_s^2) +  P6EXPECT  +  C6R4MSCL_s*GENDER + C6R4MSCL_s*RACE  + C6R4MSCL_s*WKMOMED +  P7HFAMIL + S7PUPRI + R7REGION + R7URBAN +  (1 |S7_ID), data = dat, 
                family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) # random effects logistic regression
ps.est.re <- predict(sel.re, type = 'response')   # random effects propensity score  

sel.fe <- glm(C7DESMTH ~ C6R4MSCL_s + I(C6R4MSCL_s^2) + P6EXPECT  +  C6R4MSCL_s*GENDER + C6R4MSCL_s*RACE + C6R4MSCL_s*WKMOMED +  P7HFAMIL + S7_ID, data = dat, family = binomial) # fixed effects logistic regression
ps.est.fe <- predict(sel.fe, type = 'response')   # fixed effects propensity score  

# marginal IPW estimator
lm.marginal.re = lm(C7R4MSCL ~ C7DESMTH, data=dat, weights=IPWfun(ps.est.re, dat$C7DESMTH))
lm.marginal.fe = lm(C7R4MSCL ~ C7DESMTH, data=dat, weights=IPWfun(ps.est.fe, dat$C7DESMTH))

marginal.ATEs <- rbind(Marginal.RePS=summary(lm.marginal.re)$coef[2,1:2], 
                       Marginal.FePS=summary(lm.marginal.fe)$coef[2,1:2])

# cluster IPW estimator
cluster.ATEs <- rbind(Clustered.RePS=li.clusterATE.SE(ps.est.re, dat$C7R4MSCL, dat$C7DESMTH, dat$S7_ID)$ATE,
                      Clustered.FePS=li.clusterATE.SE(ps.est.fe, dat$C7R4MSCL, dat$C7DESMTH, dat$S7_ID)$ATE)


# :: implement Causal Forests (CF) with modifications

# 1) default CF
out.cf <- causal_forest(X=covs.dum_c, Y=dat$C7R4MSCL, W=dat$C7DESMTH)
pred.cf <- predict(out.cf, type="vector", estimate.variance = TRUE)

# 2) CF+RePS
# RePS computed in the Marginal Estimator
out.cf.psRe <- causal_forest(X=covs.dum_c, Y=dat$C7R4MSCL, W=dat$C7DESMTH, W.hat=ps.est.re) # ps.est from random-effects logistic regression
pred.cf.psRe = predict(out.cf.psRe, type="vector", estimate.variance = TRUE)

# 3) CF+FePS
# FePS computed in the Marginal Estimator
out.cf.psFe <- causal_forest(X=covs.dum_c, Y=dat$C7R4MSCL, W=dat$C7DESMTH, W.hat=ps.est.fe) # ps.est from fixed-effects logistic regression
pred.cf.psFe = predict(out.cf.psFe, type="vector", estimate.variance = TRUE)

# 4) CF+IDdum
id.dummes <- model.matrix( ~ S7_ID + 0, dat)  # cluster dummies

out.cf.id.dummies <- causal_forest(X=cbind(covs.dum_c, id.dummes), Y=dat$C7R4MSCL, W=dat$C7DESMTH)
pred.cf.id.dummies <- predict(out.cf.id.dummies, type="vector", estimate.variance = TRUE)

# 5) CF+Demean 
dat$C7R4MSCL_dm <- as.matrix(dat$C7R4MSCL - ave(dat$C7R4MSCL, dat$S7_ID))
dat$C7DESMTH_dm <- as.matrix(dat$C7DESMTH - ave(dat$C7DESMTH, dat$S7_ID))

out.cf.demean <- causal_forest(X=covs_demean_c, Y=dat$C7R4MSCL_dm, W=dat$C7DESMTH_dm)
pred.cf.demean <- predict(out.cf.demean, type="vector", estimate.variance = TRUE)

# 6) CF+Demean+PS 
DmCovMat <- model.matrix(~ C6R4MSCL_s +  I(C6R4MSCL_s^2) +  P6EXPECT  +  C6R4MSCL_s*GENDER + C6R4MSCL_s*RACE + C6R4MSCL_s*WKMOMED +  P7HFAMIL + 0, data=dat)
DmCovMat <- apply(DmCovMat, 2, function(x) x - ave(x, dat$S7_ID))

DD.pred <- predict(lm(dat$C7DESMTH_dm ~ DmCovMat)) # extract demeaned PS

out.cf.demean_psDD <- causal_forest(X=covs_demean_c, Y=dat$C7R4MSCL_dm, W=dat$C7DESMTH_dm, W.hat=DD.pred)
pred.cf.demean_psDD <- predict(out.cf.demean_psDD, type="vector", estimate.variance = TRUE)


CF.ATEs <- rbind(best_linear_projection(out.cf)[1:2], best_linear_projection(out.cf.psRe)[1:2], best_linear_projection(out.cf.psFe)[1:2], 
                 best_linear_projection(out.cf.id.dummies)[1:2], best_linear_projection(out.cf.demean)[1:2], best_linear_projection(out.cf.demean_psDD)[1:2])

rownames(CF.ATEs) <- c("CF", "CF+RePS", "CF+FePS", "CF+IDdum", "CF+Demean", "CF+Demean+PS")


# :: summarize results
ATE_summary <- rbind(primafacie=summary(lm(dat$C7R4MSCL ~ dat$C7DESMTH))$coef[2,1:2],
                     DRCGEE=drgee_ATE, marginal.ATEs, cluster.ATEs, CF.ATEs) 
# in this manuscript, except for the DRCGEE estimator (and the praima facie, unadjusted estimator), standard errors were estimated using cluster bootstrap sampling with 5000 replicates. 
# you may get slightly different estimates from causal forests because each run involves random sample splitting for leaf splits and random sub-sampling to obtain multiple trees.

round(ATE_summary, 5)
