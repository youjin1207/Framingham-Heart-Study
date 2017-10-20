newphi = function(A, Y){
  ### input
  # A : an adjacency matrix for the network structure
  # Y : nominal observations
  ### output
  # Phi : (unstandardized) Phi
  
  A = as.matrix(A)
  n = length(Y)
  pvector = as.vector(table(Y)) / n
  
  B = matrix(0, nrow = nrow(A), ncol = ncol(A))
  for(i in 1:n){
      prob1 = pvector[which(Y[i] == names(table(Y)))]
      B[i,] = A[i,] / prob1
  }
  for(j in 1:n){
    prob2 = pvector[which(Y[j] == names(table(Y)))]
    B[,j] = B[,j] / prob2
  }
  
  Phi = 0
  for(i in 1:n){
    for(j in i:n){
      Phi = Phi + B[i,j]*(2*(Y[i] == Y[j]) - 1) 
    }
  }
  
  Phi = 2*Phi
  return(Phi)  
}


mean.newphi = function(s0, Y){
  ### input
  # s0 : the sum of every element in an adjacency matrix
  # Y : outcomes 
  # pvec : a sample proportions for k-categorical variable
  ### output
  # returns first moment of unstandardized phi

  n = length(Y)
  pvec = as.vector(table(Y)) / n 
  k = length(pvec)
  
  result = (s0/(n*(n-1)))*(n^2*k*(2-k) - n*sum(1/pvec))
  
  return(result)
}


m2.newphi = function(s0, s1, s2, Y){
  ### input
  # s0, s1, s2 : moments from adjacency matrix
  # Y : outcomes 
  # pvec : a sample proportions for k-categorical variable
  ### output
  # returns second moment of unstandardized phi

  n = length(Y)
  pvec = as.vector(table(Y)) / n 
  Q1 = sum(1/pvec)
  Q2 = sum(1/pvec^2)
  Q3 = sum(1/pvec^3)
  Q22 = sum((1/pvec)%*%t(1/pvec))
  k = length(pvec)
  
  E1 = (n^2*Q22 - n*Q3)/(n*(n-1))
  
  E2 = 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - 2*( 2*n^2*Q2 - n^2*k*Q2) + 2*n*Q3 - n^2*Q22
  E2 = E2/(n*(n-1)*(n-2))
  
  A1 = 4*n^4*k^2 - 4*n^4*k^3 + n^4*k^4 - (2*n^3*k*Q1 - n^3*k^2*Q1)
  A2 = 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - (2*n^2*Q2 - n^2*k*Q2)
  
  A = A1 - 2*A2
  
  B1 = 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - (2*n^2*Q2 - n^2*k*Q2)
  B2 = 2*n^2*Q2 - n^2*k*Q2 - n*Q3
  B3 = n^2*Q22 - n*Q3  
  
  B = B1 - B2 - B3
  
  C1 = 2*n^3*k*Q1 - n^3*k^2*Q1 - n^2*Q22
  C2 = 2*n^2*Q2 - n^2*k*Q2 - n*Q3
  
  C = C1 - 2*C2
  
  
  E3 = (A - 2*B - C) / (n*(n-1)*(n-2)*(n-3))
  
  result = s1*E1 + (s2 - 2*s1)*E2 + (s0^2 - s2 + s1)*E3
  
  return(result)
}


weight.matrix = function(A){
  ### input
  # A : an adjacency matrix, generally a weight matrix.
  ### output
  # returns necessary parts from adjacency matrix to derive standardized Phi,

  s0 = sum(A + t(A)) / 2
  s1 = sum((A + t(A))^2) / 2
  s2 = sum((apply(A,1,sum) + apply(A, 2, sum))^2)
  result = c(s0, s1, s2)

  return(result)
}

std.newphi = function(A, Y){
  ### input
  # A : an adjacency matrix
  # Y : observed outcomes
  ### output
  # prints out standardized phi

  rawphi = newphi(A, Y)
  weights = weight.matrix(A)
  s0 = weights[1]; s1 = weights[2]; s2 = weights[3]
  mean.rawphi = mean.newphi(s0, Y)
  var.rawphi = m2.newphi(s0, s1, s2, Y) - mean.rawphi^2
  std.phi = (rawphi - mean.rawphi) / sqrt(var.rawphi)
  return(std.phi)
}


## implement permutation test
make.permute.Phi = function(A, Y, np){
  ## input
  # A : an adjacency matrix
  # outcomes : observations (Y)
  # np : the number of permutation
  ## output
  # phi : observed statistic of phi
  # round(1-pnorm(phi),3) : p-value using asymptotic normality
  # pval : empirical p-value from np permutation.


  popn = length(Y)
  result = rep(NA, 5)
  pval = 1
  phi = std.newphi(A, Y)

  for(q in 2:np){
    set.seed(q)   
    # adjacency matrix
    newA = matrix(NA, ncol = popn, nrow = popn)
    s = sample(c(1:popn), popn, replace =FALSE)
    
    # construct the new adjacency matrix
    newA[s, s] = A
  
    pval = ifelse(std.newphi(newA,Y) >= phi , pval+1, pval)
  }
  
  pval = pval / np
  
  result = c(phi, round(1-pnorm(phi),3), pval)
  
  return(result) 
  
}