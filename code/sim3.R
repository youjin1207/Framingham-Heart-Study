library(igraph)
library(MASS)
library(netdep)
####
network_c2 = read.table("phs000153_c2.txt", sep = "\t", header = TRUE)
pheno.data = read.csv("pheno_c2_ex1_5.csv",sep = ",", header = TRUE)
pheno.info = pheno.data
focus.data = cbind(pheno.info$shareid, pheno.info$E487)
focus.data = as.data.frame(na.omit(focus.data))
names(focus.data) = c("shareid", "pheno")
focus.data = as.data.frame(na.omit(focus.data))
focus.data$pheno = ifelse(focus.data$pheno > 1, 1, focus.data$pheno)
network_c2_1.5 = network_c2
Adj = matrix(0, nrow = nrow(focus.data), ncol = nrow(focus.data))
ids = focus.data$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network_c2_1.5[which(network_c2_1.5$shareid == ids[i]), 34])
  friends2 = which(ids %in% network_c2_1.5[which(network_c2_1.5$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i ] = 1
}
sample.size = nrow(Adj)

peer.process.two = function(A, max.time = 3, mprob = 0.5, epsilon = 0.3){
  popn = nrow(A)
  # Generate the initial popn's outcomes (at t = 0)
  outcome = matrix(0, ncol = popn, nrow = (max.time+1) )
  outcome[1,] = rnorm(popn, 0, 1)
  # outcome t=1,2,...,max.time
  for (t in 2:(max.time + 1)){
    for (i in 1:popn){
      p = runif(1, 0, mprob) # susceptibility probability
      if (rowSums(A)[i] != 0){
        z = mean(outcome[t-1, which(rowSums(A)>0)]) + rnorm(1,0,epsilon)
      } else{
        z = mean(outcome[t-1,which(rowSums(A)==0)]) + rnorm(1,0,epsilon)
      }
      outcome[t,i] =  ifelse(rbinom(1, 1, p) == 1, z, outcome[t-1,i]) 
    }
  }
  outcomes = list()
  for(t in 1:(max.time + 1)){
    outcomes[[t]] = as.numeric(outcome[t,])
  }
  names(outcomes) = as.character(paste("time" ,c(0:max.time), sep = ""))
  return(outcomes)
}

epsilons = c(0.1, 0.05, 0.01)
# k = 1, 2, 3
k = 1
results = list()
for(r in 1:1000){
  set.seed(r)
  covariate = peer.process.two(Adj, max.time = 50, mprob = 0.5, epsilon = epsilons[k])$time50
  outcome = peer.process.two(Adj, max.time = 50, mprob = 0.5, epsilon = epsilons[k])$time50
  corr = cor(covariate, outcome)
  all.result = lm(outcome ~ covariate)
  all.p = summary(all.result)$coefficients[2,4]
  lower = summary(all.result)$coefficients[2,1] - 1.96*summary(all.result)$coefficients[2,2]
  upper = summary(all.result)$coefficients[2,1] + 1.96*summary(all.result)$coefficients[2,2]
  pval = as.numeric(all.p)  
  se = summary(all.result)$coefficients[2,2]
  est = summary(all.result)$coefficients[2,1]
  moran.cov = make.permute.moran(Adj, covariate, 500)
  moran.resi = make.permute.moran(Adj, all.result$residuals, 500)
  moran.outcome = make.permute.moran(Adj, outcome, 500)

  results[[r]] = list(corr = corr, moran.cov = moran.cov, moran.outcome = moran.outcome, 
                      moran.resi = moran.resi, glm.result = c(lower, upper, pval, est, se))
} 