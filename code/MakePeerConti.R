########################################
library(igraph)
library(MASS)
library(parallel)
library(doParallel)
##########################################
source("generatePeerConti.R")
source("generateLatentConti.R")
source("MoranI.R")
#######################################
popn = 200
np = 500

times = c(0,1,2,3)
outcome = matrix(0, nrow = 4, ncol = popn); result = matrix(0, nrow = 4, ncol = 3)
########## implement!
nCores <- as.numeric(Sys.getenv('NSLOTS'))
RNGkind("L'Ecuyer-CMRG")
cl <- makeCluster(nCores) # create a cluster 
registerDoParallel(cl) # register the cluster


PeerConti200 = foreach(mm = 1:500, .packages = c('igraph', 'parallel', 'MASS')) %dopar% {

  set.seed(mm)
  
  # First of all, create a dependence space
  A = latent_space_peer(popn)
  outcome = peer_influence_dependence(A, times, mprob = 0.50)
    
  for(j in 1:4){
    result[j,] = make.permute.moran(A, outcome[j,], np)
  }

  # save the outcomes for each time point 
	return(list(outcome, result)) 
} 

save(PeerConti200, file = "../data/PeerConti.RData")
