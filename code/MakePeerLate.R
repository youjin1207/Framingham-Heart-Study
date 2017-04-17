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
popn = 1000
np = 500

########## implement!
nCores <- as.numeric(Sys.getenv('NSLOTS'))
RNGkind("L'Ecuyer-CMRG")
cl <- makeCluster(nCores) # create a cluster 
registerDoParallel(cl) # register the cluster

LatentConti00 = foreach(mm = 1:500, .packages = c('igraph', 'parallel', 'MASS')) %dopar% {

  set.seed(mm)
  # First of all, create a dependence space
  tmpG = latent_variable_dependence(2*popn, 0.00)
  subg = snowball_sampling(tmpG, popn)
  G = induced.subgraph(tmpG, as.integer(subg))
  outcome =  V(G)$outcome
  A = as.matrix(get.adjacency(G))
  result= make.permute.moran(A, outcome, np)
  # save the outcomes for each time point
  
	return(list(outcome, result)) 
} 

save(LatentConti00, file = "LatentConti00.RData")

###
LatentConti10 = foreach(mm = 1:500, .packages = c('igraph', 'parallel', 'MASS')) %dopar% {

  set.seed(mm)
  # First of all, create a dependence space
  tmpG = latent_variable_dependence(2*popn, 0.10)
  subg = snowball_sampling(tmpG, popn)
  G = induced.subgraph(tmpG, as.integer(subg))
  outcome =  V(G)$outcome
  A = as.matrix(get.adjacency(G))
  result= make.permute.moran(A, outcome, np)
  # save the outcomes for each time point
  
	return(list(outcome, result)) 
} 

save(LatentConti10, file = "LatentConti10.RData")

###
LatentConti15 = foreach(mm = 1:500, .packages = c('igraph', 'parallel', 'MASS')) %dopar% {

  set.seed(mm)
  # First of all, create a dependence space
  tmpG = latent_variable_dependence(2*popn, 0.15)
  subg = snowball_sampling(tmpG, popn)
  G = induced.subgraph(tmpG, as.integer(subg))
  outcome =  V(G)$outcome
  A = as.matrix(get.adjacency(G))
  result= make.permute.moran(A, outcome, np)
  # save the outcomes for each time point
  
	return(list(outcome, result)) 
} 

save(LatentConti15, file = "LatentConti15.RData")


###
LatentConti20 = foreach(mm = 1:500, .packages = c('igraph', 'parallel', 'MASS')) %dopar% {

  set.seed(mm)
  # First of all, create a dependence space
  tmpG = latent_variable_dependence(2*popn, 0.20)
  subg = snowball_sampling(tmpG, popn)
  G = induced.subgraph(tmpG, as.integer(subg))
  outcome =  V(G)$outcome
  A = as.matrix(get.adjacency(G))
  result= make.permute.moran(A, outcome, np)
  # save the outcomes for each time point
  
	return(list(outcome, result)) 
} 

save(LatentConti20, file = "LatentConti20.RData")

