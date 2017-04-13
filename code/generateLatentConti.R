library(igraph)
library(MASS)
library(mvrtn)

latent_variable_dependence = function(popn, rho){
  ### input 
  # popn : the total number of population (N)
  # rho : the correlation between a latent variable X and Y.
  ### outout
  # igraph G : having X and Y as nodal attributes
  
  Sigma = matrix(c(1,rho,rho,1),2,2) # Covariance matrix
  sample = mvrnorm(n = popn, mu = c(0,0), Sigma ) # Generate popn's samples (X,Y)
  X  = sample[,1]; Y = sample[,2]

  prob = matrix(0, popn, popn) # initiate a probability (of sharing a tie) matrix 
  A = matrix(0, popn, popn) # initiate an adjacency matrix
  
  for (i in 1:popn) {
    for (j in i:popn) {
            if (i == j) { 
              prob[i,j] = 0.0
            }
          else if (i != j) {
            dif = abs(X[i] - X[j])
            if (dif < 0.01) dif = 0.01
            if (dif < 0.10) {
              val = 2 / (dif)
            }else if(dif < 0.50){
              val = dif
            }else if(dif < 1.00){
              val = - dif
            }else if(dif < 1.50){
              val = - 2*dif
            }else{
              val = -10*dif 
            }
          prob[i,j] <- exp(val) / (1 + exp(val))
          }
      A[i,j] = rbinom(1, 1 , prob[i,j]) 
      A[j,i] = A[i,j]
    } 
  }
  
  G = graph.adjacency(A, "undirected")
  V(G)$outcome = Y
  V(G)$latent = X
  
  return(G)  
}
