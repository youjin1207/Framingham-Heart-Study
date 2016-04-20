rm(list = ls())
library(igraph)
library(MASS)
library(mvrtn)
library(MultinomialCI)


## calculate geodesic distance in network space
## out of concern that ``shortest.paths`` does not work sometimes
all.distance <- function(G){
  # the class of G is igraph
  
  path.matrix <- shortest.paths(G, v=V(G), to = V(G))
  ad.matrix <- get.adjacency(G)
  na.set <- c()
  
  for(q in 1:length(V(G))){
    if (sum(is.na(path.matrix[q,]) ) > length(V(G))/2){
      na.set <- c(na.set, q)
    }
  }
  
  tmp.set <- na.set
  
  while(length(tmp.set) > 0){
    for(w in 1:length(na.set)){
      for(p in 1:length(V(G)) ){
        if (ad.matrix[na.set[w],p] == 1  & !(p %in% tmp.set)){
          path.matrix[na.set[w], -tmp.set] <- path.matrix[p, -tmp.set] + 1
          path.matrix[-tmp.set, na.set[w]] <- path.matrix[-tmp.set, p] + 1
          tmp.set <- tmp.set[tmp.set!=na.set[w]]
          break
        }
      }
    }
  }
  
  
  while(sum(is.na(path.matrix)) > 0){
    for(w in 1:length(na.set)){
      for(z in w:length(na.set)){
        if (is.na(path.matrix[na.set[w],na.set[z]])){
          for(p in 1:length(V(G))){
            if (ad.matrix[na.set[w], p] == 1){
              path.matrix[na.set[w], na.set[z]] <- path.matrix[p, na.set[z]] + 1
              path.matrix[na.set[z], na.set[w]] <- path.matrix[na.set[z], p] + 1
              break
            }
          }
        }
      }
    }
  }
  
  return(path.matrix)
}

## make a connection between disconnected components
connect_ngraph <- function(G){
  
  # the input G is an igraph
  
  connected_G <- G
  
  if (!is.connected(connected_G)){ # if there exists more than one component
    
    num <- no.clusters(connected_G) # number of clusters
    
    selected_node <- c() # the index vector selected from each cluster
    
    for (i in 1:num){
      
      sub <- induced.subgraph(connected_G, clusters(connected_G)$membership == i)
      
      selected_node[i] <- sample(V(sub)$name, 1) # select one sample from sub
      
      
    }
    
    
    for (j in 2:num){
      connected_G <- connected_G + 
        edges(c(V(connected_G)[V(connected_G)$name == selected_node[1]] ,
                V(connected_G)[V(connected_G)$name == selected_node[j]]))
      
    }
    
    
    
  }
  
  if (is.connected(connected_G)){
    return(connected_G)
  } else {
    return(NULL)
  }
  
}

## return confidence intervals(normal based) under independence assumption
ci_indep <- function(mat,mu){
  means=apply(mat,1,mean) # 1 indicates row
  vars=apply(mat,1,var)
  nsim=nrow(mat) # the number of iteration
  n=ncol(mat)
  lower=means-(qnorm(.975)*sqrt(vars/n))
  upper=means+(qnorm(.975)*sqrt(vars/n))
  ci=cbind(lower, upper)
  return(ci)
}

## coverage rate under independence assumption
cov_indep <- function(mat,mu){
  means=apply(mat,1,mean) # 1 indicates row
  vars=apply(mat,1,var)
  nsim=nrow(mat) # the number of iteration
  n=ncol(mat)
  lower=means-(qnorm(.975)*sqrt(vars/n))
  upper=means+(qnorm(.975)*sqrt(vars/n))
  cov=sum(lower<mu & upper>mu)/nsim
  return(cov)
}

## return pvalues
get_pvalue <- function(results, obs){
  n <- length(results) # the number of observed outcomes from each independent network
  # 'obs' is the results from the current observation
  pval <- sum(results >= obs) / n 
  return(pval) # retun the 'empirical' p value of the current observation
}

## generate latent space model - second version.
latent_space2 <- function(popn){ 
  # popn : the total number of population (node) in a single network
  
  sample <- rnorm(n = popn, 0, 1 ) # Generate popn's iid X's
  
  
  prob <- matrix(0, popn, popn) # initiate a probability (of sharing a tie) matrix
  adj <- matrix(0, popn, popn) # initiate an adjacency matrix
  
  for (i in 1:popn) {
    for (j in i:popn) {
      if (i == j) { 
        prob[i,j] = 0.0
      }
      else if (i != j) {
        dif <- abs(sample[i] - sample[j])
        
        if (dif < 0.05) dif = 0.05
        
        if (dif < 0.10) val = 1 / dif
        else val = -10 * dif
        
        prob[i,j] <- exp(val) / (1 + exp(val))
        
        
      }
      
      
      adj[i,j] <- rbinom(1, 1 , prob[i,j]) 
      adj[j,i] <- adj[i,j]
      
    }
    
    
  }
  
  return(adj)  # Return an adjancency matrix
  
}

## create a network from an adjacency matrix 
ngraph <- function(A){
  # A is an adjacency matrix of the outcomes
  # X is a data frame of all attributes

  n <- ncol(A) # n is a number of observed outcomes
  
  d<-c() 
  w<-c()
  
  for (i in 1:n){
    for (j in i:n){
      
      if (A[i,j] != 0) {
        d <- rbind(d,c(i,j))
        w <- cbind(w, A[i,j])
      } 
      
    }
  }
    
  d <- as.data.frame(d)
  # d is a data frame containing a symbolic edge list in the first two columns.
  w <- as.vector(w)
  # w is a vector of weight for each edge.
  
  # Create an igraph graph containing edge list and edge/vertex attributes.
  network_graph <- graph.data.frame(d, directed = "FALSE", vertices = c(1:n))
  
  # put weights on each edge.
  E(network_graph)$weight <-w
    
  return(network_graph)
}

## implement snowball sampling based on adjacency relationship
snowball_sampling <- function(G, samn){
  
  # G is a population igraph
  # samn is a sample size
  
  if (vcount(G) < samn){
    # exit if the population number is less than the sample size
    return("Population size is not enough for snowball sampling")
  }
  
  snow<-c()
  
  # the starter is the vertex with the largest degree
  starter <- sample(1:length(V(G)),1)
  current<-c()
  current[1] <- V(G)$name[starter]
  count <- 1
  
  snow[1] <- current[1]# vertex name with the largest degree
  
  while (count < samn ){
    
    nnode <- length(current)
    
    for (i in 1:nnode){
      
      ngh <- neighbors(G, current[i]) # vertex index
      
      snow <- c(snow, V(G)$name[ngh])
      
      snow <- unique(snow)
      
    }  
    
    tmp_sample <- snow[(count+1):length(snow)]
    
    if (samn < length(snow)){
      
      need <- samn - count # number of samples needed
      
      tmp_sample <- sample(tmp_sample ,need)
      
      snow[(count+1):samn] <- tmp_sample
      
      snow <- snow[-c((samn+1):length(snow))]
      
    }
    
    current <- tmp_sample  # a vector of indexes
    count <-length(snow)
    
    
  }
  
  if (count == samn){
    return(snow)
  } else {
    return("somthing goes wrong.")
  }
  
  
}

## generate peer influenced categorical outcomes
nominal_peer_influence <- function(G, time_point, mprob, multip){
  
  # G is a network (dependence structure)
  # time_point : a vector of time point you want to observe the outcomes
  # mprob : a maximum susceptibability probability
  
  popn <- vcount(G)
  # popn : the total number of population (node) in a single network
  
  G_matrix <- get.adjacency(G)
  
  
  max_t <- time_point[length(time_point)] + 1  # maximum number of time for observation
  max_t <- as.integer(max_t)
  
  
  outcome <- matrix(0, ncol = popn, nrow = max_t) # initialize the outcome data frame
  
  # Generate the initial popn's outcomes
  k <- length(multip)
  for(i in 1:popn){
    aa <- rmultinom(1, size = 1, multip)
    outcome[1,i] <- which(aa == 1) 
  }
  
  for (t in 2:max_t){
    for (i in 1:popn){
      
      p <- runif(1,0,mprob) # Generate a new susceptibility probability
      dummy <- rbinom(1,1,p)
      
      if(dummy == 0){
        outcome[t,i] <- outcome[(t-1),i]
      }else{
        temp <- table(G_matrix[1:popn,i]*outcome[t-1,1:popn])[-1]
        outcome[t,i] <- as.integer(names(temp)[temp == max(temp)])[1]
      }
      
    }
  }
  
  return(outcome[time_point + 1,])
}

## generate peer influenced continuous outcomes
peer_influence_dependence <- function(G, time_point, mprob){
  
  # G is a network (dependence structure)
  # time_point : a vector of time point you want to observe the outcomes
  # mprob : a maximum susceptibability probability
  
  popn <- vcount(G)
  # popn : the total number of population (node) in a single network
  
  G_matrix <- get.adjacency(G)
  
  
  max_t <- time_point[length(time_point)] + 1  # maximum number of time for observation
  max_t <- as.integer(max_t)
  
  
  outcome <- matrix(0, ncol = popn, nrow = max_t) # initialize the outcome data frame
  
  # Generate the initial popn's outcomes
  outcome[1,] <- rnorm(popn, 0, 1)
  
  for (t in 2:max_t){
    for (i in 1:popn){
      
      p <- runif(1,0,mprob) # Generate a new susceptibility probability
      
      
      deg <- degree(G, i) # degree of node i
      if (deg != 0){
        v <- sum(G_matrix[1:popn,i]*outcome[t-1,1:popn]) / deg # the average of all of subject i's neighbor's outcomes at time (t-1)
      } else
        v <- 0
      
      
      outcome[t,i] <- (1 - p)*outcome[t-1,i] + p*v
      
      
      
    }
  }
  
  return(outcome[time_point + 1,])
  
  
}

## generate latent variable dependent outcomess
new_nominal_latent <- function(popn, rho, multip){
  ## compared to nominal_latent_variable, it does not involve random ordering 
  # popn : the total number of population (node) in a single network
  # rho : the correlation between a latent variable X and a dummy variable Z
  # multip : multinomial probabiility
  Sigma <- matrix(c(1,rho, rho,1),2,2) # Covariance matrix
  sample <- mvrnorm(n = popn, mu = c(0,0), Sigma ) # Generate popn's samples (X,Y)
  
  prob <- matrix(0, popn, popn) # initiate a probability (of sharing a tie) matrix 
  adj <- matrix(0, popn, popn) # initiate an adjacency matrix
  
  for (i in 1:popn) {
    for (j in i:popn) {
      if (i == j) { 
        prob[i,j] = 0.0
      }
      else if (i != j) {
        dif <- abs(sample[i,1] - sample[j,1])
        
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
      adj[i,j] <- rbinom(1, 1 , prob[i,j]) 
      adj[j,i] <- adj[i,j]
    } 
  }
  
  # make sure that the sum of multip is one
  # make sure that the sum of multip is one
  multip <- multip / sum(multip)
  cump<-c()
  for(i in 1:length(multip)){
    cump[i] <- sum(multip[1:i])
  }
  
  setpoint <- qnorm(cump, 0, 1)
  
  nominal_value <- c()
  for(i in 1:popn){
    if(sample[i,2] <= setpoint[1]){
      nominal_value[i] <- 1
    }
    for(k in 2:length(setpoint)){
      if(sample[i,2] > setpoint[(k-1)] & sample[i,2] <= setpoint[k]){
        nominal_value[i] <- k
      } 
    }
  }
  
  
  G <- ngraph(adj)
  G <- connect_ngraph(G)
  
  V(G)$outcome <- nominal_value
  V(G)$latent <- sample[,1]
  
  return(G)  # Return an adjancency matrix
}

## calculate variance inflator (normal based)
vinflator <- function(outcomemat){
  # outcomemat is a matrix where each row for network and each column for a variable Y_{i}
  M <- nrow(outcomemat) # number of networks
  n <- ncol(outcomemat) # number of observations per one network
  varmean <- colMeans(outcomemat)
  s2 <- rep(NA, M) # sample variance of each network
  for(m in 1:M){
    s2[m] <- sum((outcomemat[m,1:n] - rowMeans(outcomemat)[m])^2) / (n-1)
  }
  covariance <- matrix(0, nrow = n , nrol = n) # covariance matrix
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      covariance[i,j] <- sum((outcomemat[1:M,i]-varmean[i])*(outcomemat[1:M,j]-varmean[j]))/(M-1)
      covariance[j,i] <- covariance[i,j]
    }
  }
  
  w <- (sum(covariance))/n
  
  result <- rep(NA, M)
  sigma2 <- rep(NA, M)
  for(m in 1:M){
    sigma2[m] <- s2[m] + w/(n-1)
    result[m] <- 1 + w/sigma2[m]
  }
  
  return(result)
}

## function to implement permutation test - version 1
make.permute1 <- function(Adj, outcomes, np){
  # Adj is a current adjacency matrix
  # outcomes is a current vector
  # np is the number of permutation
  
  popn <- length(outcomes)
  result <- rep(NA, 5)
  pval <- 1
  moran <- MoranI(Adj, outcomes)[4]
  moran.mean <- MoranI(Adj, outcomes)[2]
  moran.var <- MoranI(Adj, outcomes)[3]
  
  for(q in 2:np){
    
    # adjacency matrix
    A <- matrix(NA, ncol = popn, nrow = popn)
    s <- sample(c(1:popn), popn, replace =FALSE)
    
    # construct the new adjacency matrix
    for(j in 1:nrow(A)){
      for(k in 1:ncol(A)){
        A[s[j], s[k]] <- Adj[j,k]
      }
    }
    
    ## from adjacency matrix to igraph data
    net <- graph.adjacency(A, mode="undirected", weighted = NULL)
    A <- as.matrix(A)
    pval <- ifelse(MoranI(A, outcomes)[4] >= moran , pval+1, pval)
    
  }
  
  pval <- pval / np
  
  result<- c(moran, round(1-pnorm(moran),3) , pval, moran.mean, moran.var)
  
  ## return empirical p-value
  return(result)
  
}

## function to implement permutation test - version 2
make.permute2 <- function(Adj, outcomes, np){
  # Adj is a current adjacency matrix
  # outcomes is current outcome matrix (or vector)
  # np is the number of permutation
  
  popn <- ncol(outcomes)
  m <- nrow(outcomes)
  pval <- rep(1, m) 
  moran <- rep(0,m)
  moran.mean <- rep(0,m)
  moran.var <- rep(0,m)
  result <- matrix(NA, ncol = 5, nrow = m)
  
  for(k in 1:m){
    moran[k] <- MoranI(A_hat, outcomes[k,])[4]
    moran.mean[k] <- MoranI(A_hat, outcomes[k,])[2]
    moran.var[k] <- MoranI(A_hat, outcomes[k,])[3]
  }
  
  for(q in 1:np){
    
    # adjacency matrix
    A <- matrix(NA, ncol = popn, nrow = popn)
    s <- sample(c(1:popn), popn, replace =FALSE)
    
    # construct the new adjacency matrix
    for(j in 1:nrow(A)){
      for(k in 1:ncol(A)){
        A[s[j], s[k]] <- A_hat[j,k]
      }
    }
    
    ## from adjacency matrix to igraph data
    net <- graph.adjacency(A, mode="undirected", weighted = NULL)
    A <- as.matrix(A)
    
    for(k in 1:m){
      pval[k] <- ifelse(MoranI(A, outcomes[k,])[4] >= moran[k] , pval[k]+1, pval[k])
    }
  }
  
  pval <- pval / np
  
  for(k in 1:m){
    result[k,] <- c(moran[k], round(1-pnorm(moran[k]),3) , pval[k], moran.mean[k], moran.var[k])
  }
  
  ## return empirical p-value
  return(result)
  
}

## implement permutation test in categorical case - version 1
all.permute <- function(Adj, outcomes, np, pvec){
  
  Adj <- as.matrix(Adj)
  popn <- length(outcomes)
  pval.moran <- 1 
  pval.phi <- 1
  
  Phi <- newphi(Adj, outcomes)
  moran <- MoranI(Adj, outcomes)[4]
  
  for(q in 1:np){
    
    # adjacency matrix
    A <- matrix(NA, ncol = popn, nrow = popn)
    s <- sample(c(1:popn), popn, replace =FALSE)
    
    # construct the new adjacency matrix
    for(j in 1:nrow(A)){
      for(k in 1:ncol(A)){
        A[s[j], s[k]] <- Adj[j,k]
      }
    }
    
    ## from adjacency matrix to igraph data
    net <- graph.adjacency(A, mode="undirected", weighted = NULL)
    A <- as.matrix(A)
    
    pval.moran <- ifelse(MoranI(A, outcomes)[4] >= moran , pval.moran+1, pval.moran)
    pval.phi <- ifelse(newphi(A, outcomes) >= Phi , pval.phi+1, pval.phi)
    
  }
  
  pval.moran <- pval.moran / np
  pval.phi <- pval.phi / np
  
  result <- c(Phi, pval.phi, moran, pval.moran)
  
  return(result)
}

## implement permutation test in categorical case - version 2
all.permute2 <- function(Adj, outcomes, np){
  
  Adj <- as.matrix(Adj)
  popn <- ncol(outcomes)
  m <- nrow(outcomes)
  pval.moran <- rep(0, m)
  pval.phi <- rep(0,m)
  Phi <- rep(0,m)
  moran <- rep(0,m)
  
  for(k in 1:m){
    Phi[k] <- newphi(Adj, outcomes[k,])
    moran[k] <- MoranI(Adj, outcomes[k,])[4]
  }
  
  
  for(q in 1:np){
    
    # adjacency matrix
    A <- matrix(NA, ncol = popn, nrow = popn)
    s <- sample(c(1:popn), popn, replace =FALSE)
    
    # construct the new adjacency matrix
    for(j in 1:nrow(A)){
      for(k in 1:ncol(A)){
        A[s[j], s[k]] <- Adj[j,k]
      }
    }
    
    ## from adjacency matrix to igraph data
    net <- graph.adjacency(A, mode="undirected", weighted = NULL)
    A <- as.matrix(A)
    
    for(k in 1:m){
      pval.moran[k] <- ifelse(MoranI(A, outcomes[k,])[4] >= moran[k] , pval.moran[k]+1, pval.moran[k])
      pval.phi[k] <- ifelse(newphi(A, outcomes[k,]) >= Phi[k] , pval.phi[k]+1, pval.phi[k])
    }
  }
  
  
  pval.moran <- pval.moran / np
  pval.phi <- pval.phi / np
  
  result <- matrix(NA, nrow = m, ncol = 4)
  for(k in 1:m){
    result[k,] <- c(Phi[k], pval.phi[k], moran[k], pval.moran[k])
  }
  
  return(result)
}


## generate latent variable dependent outcomes
latent_variable_dependence2 <- function(popn , rho, dep){ 
  # n : the total number of population (node) in a single network
  # rho : the correlation between a latent variable X and an outcome Y
  
  Sigma <- matrix(c(1,rho, rho,1),2,2) # Covariance matrix
  sample <- mvrnorm(n = popn, mu = c(0,0), Sigma ) # Generate popn's samples (X,Y)
  
  # Create network graph with latent variable dependence structure
  
  prob <- matrix(0, popn, popn) # initiate a probability (of sharing a tie) matrix 
  adj <- matrix(0, popn, popn) # initiate an adjacency matrix
  
  for (i in 1:popn) {
    for (j in i:popn) {
      if (i == j) { 
        prob[i,j] = 0.0
      }
      else if (i != j) {
        dif <- abs(sample[i,1] - sample[j,1])
        
        if (dif < 0.05) dif = 0.05
        
        if (dif < 0.10) val = 1 / dif
        else val = -dep * dif
        
        prob[i,j] <- exp(val) / (1 + exp(val))
        
        
      }
      
      
      adj[i,j] <- rbinom(1, 1 , prob[i,j]) 
      adj[j,i] <- adj[i,j]
      
    }
    
    
  }
  
  if (sum(adj) != 0) {  # If at least one pair is connected
    G <- ngraph(adj)
    G <- connect_ngraph(G)
    
    V(G)$outcome <- sample[,2]
    
    return(G)  # Return an adjancency matrix
  } else{
    return("We cannot create a network")
  }
}

## Calculate Moran's I based on adjacency matrix and outcomes
MoranI <- function(A, Y){
  # A is an adjacency matrix for the network structure
  # Y is a set of outcomes
  
  m <- length(Y)
  
  s0 <- sum(A + t(A)) / 2
  s1 <- sum((A + t(A))^2) / 2
  s2 <- sum((apply(A,1,sum) + apply(A, 2, sum))^2)
  
  b2 <- sum((Y - mean(Y))^2) / m
  b4 <- sum((Y - mean(Y))^4) / m
  
  II <- t(Y - mean(Y)) %*% A %*% (Y - mean(Y)) / (s0*b2)
  II.mean <- -1/(m-1)
  II.meansq <- s1*(m*b2^2-b4)/(s0^2*b2^2*(m-1)) +
    (s2-2*s1)*(2*b4-m*b2^2) / (s0^2*b2^2*(m-1)*(m-2)) +
    (s0^2 - s2 + s1) * (3*m*b2^2-6*b4) / (s0^2*b2^2*(m-1)*(m-2)*(m-3))
  II.var <- II.meansq - II.mean^2
  II.std <- (II - II.mean) / sqrt(II.var)
  return(c(II, II.mean, II.var, II.std))
  
}

## return necessary parts from adjacency matrix to compute phi
weight.matrix <- function(A){
  # A is an adjacency matrix, generally a weight matrix
  s0 <- sum(A + t(A)) / 2
  s1 <- sum((A + t(A))^2) / 2
  s2 <- sum((apply(A,1,sum) + apply(A, 2, sum))^2)
  result <- c(s0, s1, s2)
  return(result)
}

## return unstandardized phi
newphi <- function(A, Y){
  # A is an adjacency matrix for the network structure
  # Y is a set of outcomes of nominal variable
  A <- as.matrix(A)
  n <- length(Y)
  pvector <- as.vector(table(Y)) / n
  
  Phi <- 0
  for(i in 1:n){
    for(j in 1:n){
      prob1 <- pvector[Y[i]]
      prob2 <- pvector[Y[j]]
      Phi <- Phi + A[i,j]*(2*sum(Y[i] == Y[j]) - 1)/ (prob1*prob2)
      
    }
  }
  
  return(Phi)  
}

## return first moment of unstandardized phi
mean.newphi <- function(s0, Y){
  # s0 is the sum of every element in an adjacency matrix
  # Y is an outcomes 
  # pvec is a sample proportions for k-categorical variable
  n <- length(Y)
  pvec <- as.vector(table(Y)) / n 
  k <- length(pvec)
  
  result <- (s0/(n*(n-1)))*(n^2*k*(2-k) - n*sum(1/pvec))
  
  return(result)
}

## return second moment of unstandardized phi
m2.newphi <- function(s0, s1, s2, Y){
  n <- length(Y)
  pvec <- as.vector(table(Y)) / n 
  Q1 <- sum(1/pvec)
  Q2 <- sum(1/pvec^2)
  Q3 <- sum(1/pvec^3)
  Q22 <- sum((1/pvec)%*%t(1/pvec))
  k <- length(pvec)
  
  E1 <- (n^2*Q22 - n*Q3)/(n*(n-1))
  
  E2 <- 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - 2*( 2*n^2*Q2 - n^2*k*Q2) + 2*n*Q3 - n^2*Q22
  E2 <- E2/(n*(n-1)*(n-2))
  
  A1 <- 4*n^4*k^2 - 4*n^4*k^3 + n^4*k^4 - (2*n^3*k*Q1 - n^3*k^2*Q1)
  A2 <- 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - (2*n^2*Q2 - n^2*k*Q2)
  
  A <- A1 - 2*A2
  
  B1 <- 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - (2*n^2*Q2 - n^2*k*Q2)
  B2 <- 2*n^2*Q2 - n^2*k*Q2 - n*Q3
  B3 <- n^2*Q22 - n*Q3  
  
  B <- B1 - B2 - B3
  
  C1 <- 2*n^3*k*Q1 - n^3*k^2*Q1 - n^2*Q22
  C2 <- 2*n^2*Q2 - n^2*k*Q2 - n*Q3
  
  C <- C1 - 2*C2
  
  
  E3 <- (A - 2*B - C) / (n*(n-1)*(n-2)*(n-3))
  
  result <- s1*E1 + (s2 - 2*s1)*E2 + (s0^2 - s2 + s1)*E3
  
  return(result)
}




save.image("data/functions.RData")