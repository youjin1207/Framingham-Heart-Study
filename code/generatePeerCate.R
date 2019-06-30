library(igraph)
library(mvrtn)


nominal_peer_influence = function(A, time_point, mprob, multip){
  ### input
  # A : an adjacency matrix.
  # time_point : a vector of time point you want to observe the outcomes.
  # mprob : a maximum susceptibability probability.
  # multip : multinomial probability.
  ### output
  #  an array of observations at time_point.
  
  popn = nrow(A)
  # popn : the number of subjects (n)
  
  max_t = time_point[length(time_point)] + 1  # maximum number of time for observation
  max_t = as.integer(max_t)

  outcome = matrix(0, ncol = popn, nrow = max_t) # initialize the outcome data frame
  
  # Generate observations at initial stage
  k = length(multip)
  for(i in 1:popn){
    dumb = rmultinom(1, size = 1, multip)
    outcome[1,i] = which(dumb == 1) 
  }
  
  for (t in 2:max_t){
    for (i in 1:popn){  
      p = runif(1,0,mprob) # Generate a new susceptibility probability
      dummy = rbinom(1,1,p)
      
      if(dummy == 0){
        outcome[t,i] = outcome[(t-1),i]
      }else{
        temp = table(A[1:popn,i]*outcome[t-1,1:popn])[-1]
        if(length(temp) == 0){
          outcome[t,i] = outcome[(t-1), i]
        }else{
          outcome[t,i] =  sample(outcome[(t-1),which(A[i,] == 1)] ,1)  
        }
      }
        
    }
  }
  
  return(outcome[time_point + 1,])
}