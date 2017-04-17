## generate latent space model for continuous observations in Peer Influence dependence
latent_space_peer = function(popn){ 
  ### input
  # popn : the number of subjects.
  ### output
  # A : (popn x popn) adjacency matrix.
  
  X <- rnorm(n = popn, 0, 1) # benerate popn's iid X's
   
  prob = matrix(0, popn, popn) # initiate a probability (of sharing a tie) matrix
  A = matrix(0, popn, popn) # initiate an adjacency matrix
  
  for (i in 1:popn){
    for (j in i:popn){
      if (i == j) { 
        prob[i,j] = 0.0
      }
      else if (i != j) {
        dif = abs(X[i] - X[j])
        
        if (dif < 0.05) dif = 0.05
        
        if (dif < 0.10) val = 1 / dif
        else val = -10 * dif
        
        prob[i,j] <- exp(val) / (1 + exp(val)) 
      }
      
      A[i,j] = rbinom(1, 1 , prob[i,j]) 
      A[j,i] = A[i,j]     
    } 
  }
  
  return(A)   
}


peer_influence_dependence = function(A, time_point, mprob){
  ### input
  # A : an adjacency matrix.
  # time_point : a vector of time point you want to observe the outcomes.
  # mprob : a maximum susceptibability probability
  ### output
  # an array of observations at time_point.


  popn = nrow(A)
  
  max_t = time_point[length(time_point)] + 1  # maximum number of time for observation
  max_t = as.integer(max_t)
  
  outcome <- matrix(0, ncol = popn, nrow = max_t) # initialize the outcome data frame
  
  # Generate the initial popn's outcomes (at t = 0)
  outcome[1,] = rnorm(popn, 0, 1)
  
  for (t in 2:max_t){
    for (i in 1:popn){
        
        p = runif(1,0, mprob) # Generate a new susceptibility probability
        deg = rowSums(A)[i] # the number of friends of subject i
        if (deg != 0){
          # the average of all of subject i's neighbor's outcomes at time (t-1)
          z = sum(A[1:popn,i]*outcome[t-1,1:popn]) / deg 
        } else
          z = outcome[t-1, i]
       
      outcome[t,i] = (1 - p)*outcome[t-1,i] + p*z + rnorm(1, 0, 0.3)
    }
  }
  
  return(outcome[time_point + 1,])
}

