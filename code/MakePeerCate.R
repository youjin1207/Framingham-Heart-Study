########################################
library(igraph)
library(MASS)
library(parallel)
library(doParallel)
##########################################
source("generatePeerConti.R")
source("generateLatentConti.R")
source("generatePeerCate.R")
source("snowball.R")
source("Phi.R")
#######################################
popn = 200
np = 500
########## implement!
nCores <- as.numeric(Sys.getenv('NSLOTS'))
RNGkind("L'Ecuyer-CMRG")
cl <- makeCluster(nCores) # create a cluster 
registerDoParallel(cl) # register the cluster


PeerCate = foreach(mm = 1:500, .packages = c('igraph', 'parallel')) %dopar% {

  set.seed(mm)
  
  # First of all, create a dependence space
  A = latent_space_peer(popn)
  
  outcome = nominal_peer_influence(A, c(0,1,2,3), 0.40, c(0.1, 0.2, 0.3, 0.25, 0.15))
  
  # save the outcomes for each time point
  result.t0 = make.permute.Phi(A, outcome[1,], np)
  result.t1 = make.permute.Phi(A, outcome[2,], np)
  result.t2 = make.permute.Phi(A, outcome[3,], np)
  result.t3 = make.permute.Phi(A, outcome[4,], np)

  result = rbind(result.t0, result.t1, result.t2, result.t3)

  
  return(list(outcome, result)) 
} 

save(PeerCate, file = "../data/PeerCate.RData")

