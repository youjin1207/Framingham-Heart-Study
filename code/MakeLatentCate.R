########################################
library(igraph)
library(MASS)
library(parallel)
library(doParallel)
##########################################
source("generatePeerConti.R")
source("generateLatentConti.R")
source("generateLatentCate.R")
source("snowball.R")
source("Phi.R")
#######################################
popn = 200
np = 500

rhos = c(0.00, 0.30, 0.40, 0.50)
outcome = matrix(0, nrow = 4, ncol = popn); result = matrix(0, nrow = 4, ncol = 3)
########## implement!
nCores <- as.numeric(Sys.getenv('NSLOTS'))
RNGkind("L'Ecuyer-CMRG")
cl <- makeCluster(nCores) # create a cluster 
registerDoParallel(cl) # register the cluster


LatentCate = foreach(mm = 1:500, .packages = c('igraph', 'parallel', 'MASS')) %dopar% {

  set.seed(mm)
  
  # First of all, create a dependence space
  for(j in 1:4){
    set.seed(mm*j)
    tmpG = nominal_latent_dependence(2*popn, rho = rhos[j], multip = c(0.1, 0.2, 0.3, 0.25, 0.15))
    subg = snowball_sampling(tmpG, popn)
    G = induced.subgraph(tmpG, as.integer(subg))
    outcome[j,] =  V(G)$outcome
    A = as.matrix(get.adjacency(G))
    result[j,] = make.permute.Phi(A, outcome[j,], np)
  }
  # save the outcomes for each time point
  
  return(list(outcome, result)) 
} 

save(LatentCate, file = "../data/LatentCate.RData")

