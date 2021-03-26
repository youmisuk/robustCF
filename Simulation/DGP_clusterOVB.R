###############################################################################################
# Tunning Random Forests for Causal Inference Under Cluster-level Unmeasured Confounding
# Youmi Suk and Hyunseung Kang 

# ::: Data Generating Models :::
# twolevel.normalU
# twolevel.uniformU
# twolevel.binary 
# twolevel.experror

# :: Design 1: Normally Distributed U_j and Linear Outcome Model ####

twolevel.normalU <- function(ids=1, Smpl.size = "150.30", s.val=0.3, o.val=2, m.val=2) {
  
  # ::::: 1) generate the number of observations in each cluster :::::
  
  size <- strsplit(Smpl.size, "\\.")[[1]]  # the number of cluster and cluster sizes
  
  J <- as.numeric(size[1]) # level-2 unit, the number of cluster
  n.clus <- as.numeric(size[2]) # level-1 unit, cluster sizes
  
  if (n.clus < 20) {  # variation of cluster sizes 
    sd.clus <- 0.5 
  } else if (n.clus >= 20 & n.clus < 30) { 
    sd.clus <- 1 
  } else if (n.clus >= 30) {
    sd.clus <- 2
  }
  
  N <- round(rnorm(J, n.clus, sd.clus))                   # normally distributed cluster sizes 
  id <- as.factor(rep(ids:(J + ids - 1), N))              # cluster id
  
  # ::::: 2) generate level-1 covariates, Xs with cluster-specific means :::::
  
  totalN <- length(id)
  X1 <- runif(totalN, -1, 1)
  X2 <- rnorm(totalN)
  X3 <- runif(totalN, 0, 1)
  
  # ::::: 3) generate level-2 covariates, W :::::
  
  W1 <- runif(J, -1, 1)
  W2 <- rnorm(J, 0, 1) # W2 is an unmeasured cluster-level confounder.	 
  names(W1) <- names(W2)  <- levels(id) 
  
  pop <- data.frame(id, X1, X2, X3, W1=W1[id], W2=W2[id])
  
  # ::::: 4) generate selection probabilities and potential outcome ::::::::
  E <- rnorm(sum(N), 0, 1)   # error terms for pot.  
  pop$lps <- -0.6 + 0.3* pop$X1 + 0.3* pop$X2 + 0.3*pop$X3 + 0.3* pop$W1 + s.val* (pop$W2) + 0.4* pop$X2*pop$W1 + 0.4 * pop$X1^2 +  0.4 * pop$X1*I(pop$X3 < 0.3)   # ps logit
  pop$Y0 <- 70  + 2*(pop$X1 + pop$X2 + pop$X3 + pop$W1 + pop$X2*pop$W1 + pop$X1^2 + pop$X1*I( pop$X3  < 0.3)) + o.val*(pop$W2)  + E
  pop$Y1 <- pop$Y0 + 2 + 2*pop$X3 + 2*pop$W1 + m.val*(pop$W2)
  
  pop$ps <- 1 / (1 + exp(-pop$lps))   # ps
  pop$Z <- rbinom(nrow(pop), 1, pop$ps) # treatment indicator  
  pop$Y <- with(pop, ifelse(Z == 1, Y1, Y0))
  
  pop
} 


# :: Designs 2 & 3: Uniformly Distributed U_j and Linear Outcome Model ####

twolevel.uniformU <- function(ids=1, Smpl.size = "150.30", s.val=0.3, o.val=2, m.val=2) {
  
  # ::::: 1) generate the number of observations in each cluster :::::
  
  size <- strsplit(Smpl.size, "\\.")[[1]]  # the number of cluster and cluster sizes
  
  J <- as.numeric(size[1]) # level-2 unit, the number of cluster
  n.clus <- as.numeric(size[2]) # level-1 unit, cluster sizes
  
  if (n.clus < 20) {  # variation of cluster sizes 
    sd.clus <- 0.5 
  } else if (n.clus >= 20 & n.clus < 30) { 
    sd.clus <- 1 
  } else if (n.clus >= 30) {
    sd.clus <- 2
  }
  
  N <- round(rnorm(J, n.clus, sd.clus))                   # normally distributed cluster sizes 
  id <- as.factor(rep(ids:(J + ids - 1), N))              # cluster id
  
  # ::::: 2) generate level-1 covariates, Xs with cluster-specific means :::::
  
  totalN <- length(id)
  X1 <- runif(totalN, -1, 1)
  X2 <- rnorm(totalN)
  X3 <- runif(totalN, 0, 1)
  
  # ::::: 3) generate level-2 covariates, W :::::
  
  W1 <- runif(J, -1, 1)
  W2 <- runif(J, -2, 2) # W2 is an unmeasured cluster-level confounder.	  
  names(W1) <- names(W2)  <- levels(id) 
  
  pop <- data.frame(id, X1, X2, X3, W1=W1[id], W2=W2[id])
  
  # ::::: 4) generate selection probabilities and potential outcome ::::::::
  
  E <- rnorm(sum(N), 0, 1)   # error terms for pot.  
  
  pop$lps <- -0.6 + 0.3* pop$X1 + 0.3* pop$X2 + 0.3*pop$X3 + 0.3* pop$W1 + s.val* (pop$W2) + 0.4* pop$X2*pop$W1 + 0.4 * pop$X1^2 +  0.4 * pop$X1*I(pop$X3 < 0.3)   # ps logit
  pop$Y0 <- 70  + 2*(pop$X1 + pop$X2 + pop$X3 + pop$W1 + pop$X2*pop$W1 + pop$X1^2 + pop$X1*I( pop$X3  < 0.3)) + o.val*(pop$W2) + E  
  pop$Y1 <- pop$Y0 + 2 + 2*pop$X3 + 2*pop$W1 + m.val*(pop$W2)
  
  
  pop$ps <- 1 / (1 + exp(-pop$lps))   # ps
  pop$Z <- rbinom(nrow(pop), 1, pop$ps) # treatment indicator  
  pop$Y <- with(pop, ifelse(Z == 1, Y1, Y0))
  
  pop
} 


# :: Design 4: Uniformly Distributed U_j and Nonlinear Outcome Model ####
twolevel.binary <- function(ids=1, Smpl.size = "150.30", s.val=0.3, o.val=0.3, m.val=0.3) {
  
  # ::::: 1) generate the number of observations in each cluster :::::
  
  size <- strsplit(Smpl.size, "\\.")[[1]]  # the number of cluster and cluster sizes
  
  J <- as.numeric(size[1]) # level-2 unit, the number of cluster
  n.clus <- as.numeric(size[2]) # level-1 unit, cluster sizes
  
  if (n.clus < 20) {  # variation of cluster sizes 
    sd.clus <- 0.5 
  } else if (n.clus >= 20 & n.clus < 30) { 
    sd.clus <- 1 
  } else if (n.clus >= 30) {
    sd.clus <- 2
  }
  
  N <- round(rnorm(J, n.clus, sd.clus))                   # normally distributed cluster sizes 
  id <- as.factor(rep(ids:(J + ids - 1), N))              # cluster id
  
  # ::::: 2) generate level-1 covariates, Xs with cluster-specific means :::::
  
  totalN <- length(id)
  X1 <- runif(totalN, -1, 1)
  X2 <- rnorm(totalN)
  X3 <- runif(totalN, 0, 1)
  
  # ::::: 3) generate level-2 covariates, W :::::
  
  W1 <- runif(J, -1, 1)
  W2 <- runif(J, -2, 2)  # W2 is an unmeasured cluster-level confounder.	  
  names(W1) <- names(W2)  <- levels(id) 
  
  pop <- data.frame(id, X1, X2, X3, W1=W1[id], W2=W2[id])
  
  # ::::: 4) generate selection probabilities and potential outcome ::::::::
  pop$lps <- -0.6 + 0.3* pop$X1 + 0.3* pop$X2 + 0.3*pop$X3 + 0.3* pop$W1 + s.val* (pop$W2) + 0.4* pop$X2*pop$W1 + 0.4 * pop$X1^2 +  0.4 * pop$X1*I(pop$X3 < 0.3)   # ps logit
  pop$lpsY0 <- -0.6  + 0.3*(pop$X1 + pop$X2 + pop$X3 + pop$W1 + pop$X2*pop$W1 + pop$X1^2 + pop$X1*I( pop$X3  < 0.3)) + o.val*(pop$W2)
  pop$lpsY1 <- pop$lpsY0 + 0.5 + 0.3*pop$X3 + 0.3*pop$W1 + m.val*(pop$W2)
  
  pop$ps <- 1 / (1 + exp(-pop$lps))   # ps
  pop$psY0 <- 1 / (1 + exp(-pop$lpsY0))    
  pop$psY1 <- 1 / (1 + exp(-pop$lpsY1))      
  pop$Z <- rbinom(nrow(pop), 1, pop$ps) # treatment indicator  
  pop$Y0 <- rbinom(nrow(pop), 1, pop$psY0)   # potential control outcome
  pop$Y1 <- rbinom(nrow(pop), 1, pop$psY1)   # potential treatment outcome
  pop$Y <- with(pop, ifelse(Z == 1, Y1, Y0))
  
  pop
} 

# :: Design 5: Uniformly Distributed $U_j$ and Exponential Errors in the Outcome Model ####

twolevel.experror <- function(ids=1, Smpl.size = "150.30", s.val=0.3, o.val=2, m.val=2) {
  
  # ::::: 1) generate the number of observations in each cluster :::::
  
  size <- strsplit(Smpl.size, "\\.")[[1]]  # the number of cluster and cluster sizes
  
  J <- as.numeric(size[1]) # level-2 unit, the number of cluster
  n.clus <- as.numeric(size[2]) # level-1 unit, cluster sizes
  
  if (n.clus < 20) {  # variation of cluster sizes 
    sd.clus <- 0.5 
  } else if (n.clus >= 20 & n.clus < 30) { 
    sd.clus <- 1 
  } else if (n.clus >= 30) {
    sd.clus <- 2
  }
  
  N <- round(rnorm(J, n.clus, sd.clus))                   # normally distributed cluster sizes 
  id <- as.factor(rep(ids:(J + ids - 1), N))              # cluster id
  
  # ::::: 2) generate level-1 covariates, Xs with cluster-specific means :::::
  
  totalN <- length(id)
  X1 <- runif(totalN, -1, 1)
  X2 <- rnorm(totalN)
  X3 <- runif(totalN, 0, 1)
  
  # ::::: 3) generate level-2 covariates, W :::::
  
  W1 <- runif(J, -1, 1)
  W2 <- runif(J, -2, 2)  # W2 is an unmeasured cluster-level confounder.	 	 
  names(W1) <- names(W2)  <- levels(id) 
  
  pop <- data.frame(id, X1, X2, X3, W1=W1[id], W2=W2[id])
  
  # ::::: 4) generate selection probabilities and potential outcome ::::::::
  
  E <- rexp(sum(N), rate = 1/5) # follow exponential dist.
  
  pop$lps <- -0.6 + 0.3* pop$X1 + 0.3* pop$X2 + 0.3*pop$X3 + 0.3* pop$W1 + s.val* (pop$W2) + 0.4* pop$X2*pop$W1 + 0.4 * pop$X1^2 +  0.4 * pop$X1*I(pop$X3 < 0.3)   # ps logit
  pop$Y0 <- 70  + 2*(pop$X1 + pop$X2 + pop$X3 + pop$W1 + pop$X2*pop$W1 + pop$X1^2 + pop$X1*I( pop$X3  < 0.3)) + o.val*(pop$W2) + E  
  pop$Y1 <- pop$Y0 + 2 + 2*pop$X3 + 2*pop$W1 + m.val*(pop$W2)
  
  
  pop$ps <- 1 / (1 + exp(-pop$lps))   # ps
  pop$Z <- rbinom(nrow(pop), 1, pop$ps) # treatment indicator  
  pop$Y <- with(pop, ifelse(Z == 1, Y1, Y0))
  
  pop
} 




