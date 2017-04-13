library(igraph)
library(mvrtn)

nominal_latent_dependence = function(popn, rho, multip){
  ### input 
  # popn : the total number of population (N) (before sampling)
  # rho : the correlation between a latent variable X and tmpY (continuous version).
  # multip : multinomial probability each connectes tmpY and categorical Y
  ### outout
  # igraph G : having X and Y as nodal attributes
  
  Sigma = matrix(c(1,rho,rho,1),2,2) # Covariance matrix
  sample = mvrnorm(n = popn, mu = c(0,0), Sigma ) # Generate popn's samples (X,Y)
  X  = sample[,1]; tmpY = sample[,2]

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
      A[i,j] <- rbinom(1, 1 , prob[i,j]) 
      A[j,i] <- A[i,j]
    } 
  }
  
  
  multip = multip / sum(multip)
  cump = c()
  for(i in 1:length(multip)){
    cump[i] = sum(multip[1:i])
  }

  setpoint = qnorm(cump, 0, 1)
  
  Y = c()
  for(i in 1:popn){
  	if(tmpY[i] <= setpoint[1]){
  		Y[i] = 1
  	}
  	for(k in 2:length(setpoint)){
  		if(tmpY[i] > setpoint[(k-1)] & tmpY[i] <= setpoint[k]){
  			Y[i] = k
  		} 
  	}
  }
  
  G = graph.adjacency(A, "undirected")
 
  V(G)$outcome = Y
  V(G)$latent = X
  
  return(G)  
}
