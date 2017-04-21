########################################
library(igraph)
library(MASS)
library(parallel)
library(doParallel)
##########################################
source("generatePeerConti.R")
source("generateLatentConti.R")
source("snowball.R")
source("MoranI.R")
#######################################
popn = 200
np = 500

rhos = c(0.00, 0.20, 0.30, 0.40)
outcome = matrix(0, nrow = 4, ncol = popn); result = matrix(0, nrow = 4, ncol = 3)
########## implement!
nCores <- as.numeric(Sys.getenv('NSLOTS'))
RNGkind("L'Ecuyer-CMRG")
cl <- makeCluster(nCores) # create a cluster 
registerDoParallel(cl) # register the cluster


LatentConti = foreach(mm = 1:500, .packages = c('igraph', 'parallel', 'MASS')) %dopar% {

  set.seed(mm)
  
  # First of all, create a dependence space
  for(j in 1:4){
    set.seed(mm*j)
    tmpG = latent_variable_dependence(2*popn, rhos[j])
    subg = snowball_sampling(tmpG, popn)
    G = induced.subgraph(tmpG, as.integer(subg))
    outcome[j,] =  V(G)$outcome
    A = as.matrix(get.adjacency(G))
    result[j,] = make.permute.moran(A, outcome[j,], np)
  }
  # save the outcomes for each time point
  
  return(list(outcome, result)) 
} 

save(LatentConti, file = "../data/LatentConti.RData")

