library(igraph)

MoranI = function(A, Y){
  ## input
  # A : an adjacency matrix for the network structure
  # Y : observations (Y)
  
  m = length(Y)
  
  s0 = sum(A + t(A)) / 2
  s1 = sum((A + t(A))^2) / 2
  s2 = sum((apply(A,1,sum) + apply(A, 2, sum))^2)
  
  b2 = sum((Y - mean(Y))^2) / m
  b4 = sum((Y - mean(Y))^4) / m
  
  II = t(Y - mean(Y)) %*% A %*% (Y - mean(Y)) / (s0*b2)
  II.mean = -1/(m-1)
  II.meansq = s1*(m*b2^2-b4)/(s0^2*b2^2*(m-1)) +
    (s2-2*s1)*(2*b4-m*b2^2) / (s0^2*b2^2*(m-1)*(m-2)) +
    (s0^2 - s2 + s1) * (3*m*b2^2-6*b4) / (s0^2*b2^2*(m-1)*(m-2)*(m-3))
  II.var = II.meansq - II.mean^2
  II.std = (II - II.mean) / sqrt(II.var)
  return(c(II, II.mean, II.var, II.std))
  
}

## implement permutation test
make.permute.moran = function(A, Y, np){
  ## input
  # A : an adjacency matrix
  # outcomes : observations (Y)
  # np : the number of permutation
  ## output
  # moran : observed statistic of Moran's I
  # round(1-pnorm(moran),3) : p-value using asymptotic normality
  # pval : empirical p-value from np permutation.


  popn = length(Y)
  result = rep(NA, 5)
  pval = 1
  moran = MoranI(A, Y)[4]

  for(q in 2:np){
    set.seed(q) 	
    # adjacency matrix
    newA = matrix(NA, ncol = popn, nrow = popn)
    s = sample(c(1:popn), popn, replace =FALSE)
    
    # construct the new adjacency matrix
   	newA[s, s] = A
  
    pval = ifelse(MoranI(newA, Y)[4] >= moran , pval+1, pval)
  }
  
  pval = pval / np
  
  result = c(moran, round(1-pnorm(moran),3), pval)
  
  return(result) 
}