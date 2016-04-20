###########################################
library(igraph)
library(MASS)
library(mvrtn)
library(MultinomialCI)
library(xtable)
###########################################
load("data/functions.RData")
load("data/categorical_peer.RData")

popn <- 1000
np <- 500
# weight matrix
weight.matrix_peer <- matrix(NA, ncol = 3, nrow = np)
# outcomes
cate_peer_t0 <- matrix(NA, ncol = popn, nrow = np)
cate_peer_t1 <- matrix(NA, ncol = popn, nrow = np)
cate_peer_t2 <- matrix(NA, ncol = popn, nrow = np)
cate_peer_t3 <- matrix(NA, ncol = popn, nrow = np)

# test statistics and their pvalues
result_peer_t0 <- matrix(NA, ncol = 2, nrow = np)
result_peer_t1 <- matrix(NA, ncol = 2, nrow = np)
result_peer_t2 <- matrix(NA, ncol = 2, nrow = np)
result_peer_t3 <- matrix(NA, ncol = 2, nrow = np)


for(i in 1:25){
  # weight matrix
  weight.matrix_peer[(20*(i-1)+1):(20*i),] <- categorical_peer_list[[i]][[1]]
  
  # outcome matrix
  cate_peer_t0[(20*(i-1)+1):(20*i),] <- categorical_peer_list[[i]][[2]]
  cate_peer_t1[(20*(i-1)+1):(20*i),] <- categorical_peer_list[[i]][[3]]
  cate_peer_t2[(20*(i-1)+1):(20*i),] <- categorical_peer_list[[i]][[4]]
  cate_peer_t3[(20*(i-1)+1):(20*i),] <- categorical_peer_list[[i]][[5]]
  
  # test statistics and their pvalues
  result_peer_t0[(20*(i-1)+1):(20*i),] <- categorical_peer_list[[i]][[6]]
  result_peer_t1[(20*(i-1)+1):(20*i),] <- categorical_peer_list[[i]][[7]]
  result_peer_t2[(20*(i-1)+1):(20*i),] <- categorical_peer_list[[i]][[8]]
  result_peer_t3[(20*(i-1)+1):(20*i),] <- categorical_peer_list[[i]][[9]]
}

##############################################################################
## confidence interval
library(MultinomialCI)
truth <- c(0.1,0.2,0.3,0.25,0.15)
count <- rep(0, 4)
catet0_table <- matrix(NA, nrow = 500, ncol = 5)
catet1_table <- matrix(NA, nrow = 500, ncol = 5)
catet2_table <- matrix(NA, nrow = 500, ncol = 5)
catet3_table <- matrix(NA, nrow = 500, ncol = 5)



for(i in 1:np){
  catet0_table[i,] <- as.numeric(table(cate_peer_t0[i,]))
  catet1_table[i,] <- as.numeric(table(cate_peer_t1[i,]))
  catet2_table[i,] <- as.numeric(table(cate_peer_t2[i,]))
  catet3_table[i,] <- as.numeric(table(cate_peer_t3[i,]))
}



count <- rep(0,4)
for(i in 1:500){
  est0 <- multinomialCI(catet0_table[i,], 0.05)
  if(sum(est0[,1] <= truth)==5 & sum(truth <= est0[,2])==5){
    count[1] <- count[1] +1
  }
  est1 <- multinomialCI(catet1_table[i,], 0.05)
  if(sum(est1[,1] <= truth)==5 & sum(truth <= est1[,2])==5){
    count[2] <- count[2] +1
  }
  est2 <- multinomialCI(catet2_table[i,], 0.05)
  if(sum(est2[,1] <= truth)==5 & sum(truth <= est2[,2])==5){
    count[3] <- count[3] +1
  }
  est3 <- multinomialCI(catet3_table[i,], 0.05)
  if(sum(est3[,1] <= truth)==5 & sum(truth <= est3[,2])==5){
    count[4] <- count[4] +1
  }
}


#### standardize phi
stdphit0 <- rep(NA, np)
stdphit1 <- rep(NA, np)
stdphit2 <- rep(NA, np)
stdphit3 <- rep(NA, np)

## t = 0
for(i in 1:np){
  a <- m2.newphi(weight.matrix_peer[i,1], weight.matrix_peer[i,2], 
                 weight.matrix_peer[i,3],  cate_peer_t0[i,]) 
  b <- mean.newphi(weight.matrix_peer[i,1], cate_peer_t0[i,])
  stdphit0[i] <- result_peer_t0[i] - b
  stdphit0[i] <- stdphit0[i] / sqrt(a - b^2)
}


## t = 1
for(i in 1:np){
  a <- m2.newphi(weight.matrix_peer[i,1], weight.matrix_peer[i,2], 
                 weight.matrix_peer[i,3],  cate_peer_t1[i,]) 
  b <- mean.newphi(weight.matrix_peer[i,1], cate_peer_t1[i,])
  stdphit1[i] <- result_peer_t1[i] - b
  stdphit1[i] <- stdphit1[i] / sqrt(a - b^2)
}



## t = 2
for(i in 1:np){
  a <- m2.newphi(weight.matrix_peer[i,1], weight.matrix_peer[i,2], 
                 weight.matrix_peer[i,3],  cate_peer_t2[i,]) 
  b <- mean.newphi(weight.matrix_peer[i,1], cate_peer_t2[i,])
  stdphit2[i] <- result_peer_t2[i] - b
  stdphit2[i] <- stdphit2[i] / sqrt(a - b^2)
}



## t = 3
for(i in 1:np){
  a <- m2.newphi(weight.matrix_peer[i,1], weight.matrix_peer[i,2], 
                 weight.matrix_peer[i,3],  cate_peer_t3[i,]) 
  b <- mean.newphi(weight.matrix_peer[i,1], cate_peer_t3[i,])
  stdphit3[i] <- result_peer_t3[i] - b
  stdphit3[i] <- stdphit3[i] / sqrt(a - b^2)
}


## summary table
catepeer_summary <- matrix(NA, nrow = 4, ncol = 4)

rownames(catepeer_summary) <- c("t=0", "t=1", "t=2", "t=3")
colnames(catepeer_summary) <- c("95% CI coverage rate",  expression(phi), 
                                "pvalue(z)", "pvalue(permutation)")
catepeer_summary[,1] <- count / np

catepeer_summary[,2] <- c(mean(stdphit0), mean(stdphit1), mean(stdphit2), mean(stdphit3))

catepeer_summary[,3] <- c( sum( (1-pnorm(stdphit0)) <= 0.05 ) / 500, sum( (1-pnorm(stdphit1)) <= 0.05 ) / 500,
                           sum( (1-pnorm(stdphit2)) <= 0.05 ) / 500, sum( (1-pnorm(stdphit3)) <= 0.05 ) / 500)
catepeer_summary[,3] <- catepeer_summary[,3]*100

catepeer_summary[,4] <- c( sum( result_peer_t0[,2] <= 0.05 ) / 500, 
                           sum( result_peer_t1[,2] <= 0.05 ) / 500,
                           sum( result_peer_t2[,2] <= 0.05 ) / 500, 
                           sum( result_peer_t3[,2] <= 0.05 ) / 500)
catepeer_summary[,4] <- catepeer_summary[,4]*100


catepeer_summary <- as.data.frame(catepeer_summary)
print(xtable(catepeer_summary , digits = 2, row.names = TRUE))



############### pm = 0.30 ###########
load("data/categorical_peer2.RData")

popn <- 1000
np <- 500
# weight matrix
weight.matrix_peer <- matrix(NA, ncol = 3, nrow = np)
# outcomes
cate_peer_t0 <- matrix(NA, ncol = popn, nrow = np)
cate_peer_t1 <- matrix(NA, ncol = popn, nrow = np)
cate_peer_t2 <- matrix(NA, ncol = popn, nrow = np)
cate_peer_t3 <- matrix(NA, ncol = popn, nrow = np)

# test statistics and their pvalues
result_peer_t0 <- matrix(NA, ncol = 2, nrow = np)
result_peer_t1 <- matrix(NA, ncol = 2, nrow = np)
result_peer_t2 <- matrix(NA, ncol = 2, nrow = np)
result_peer_t3 <- matrix(NA, ncol = 2, nrow = np)


for(i in 1:25){
  # weight matrix
  weight.matrix_peer[(20*(i-1)+1):(20*i),] <- categorical_peer_list2[[i]][[1]]
  
  # outcome matrix
  cate_peer_t0[(20*(i-1)+1):(20*i),] <- categorical_peer_list2[[i]][[2]]
  cate_peer_t1[(20*(i-1)+1):(20*i),] <- categorical_peer_list2[[i]][[3]]
  cate_peer_t2[(20*(i-1)+1):(20*i),] <- categorical_peer_list2[[i]][[4]]
  cate_peer_t3[(20*(i-1)+1):(20*i),] <- categorical_peer_list2[[i]][[5]]
  
  # test statistics and their pvalues
  result_peer_t0[(20*(i-1)+1):(20*i),] <- categorical_peer_list2[[i]][[6]]
  result_peer_t1[(20*(i-1)+1):(20*i),] <- categorical_peer_list2[[i]][[7]]
  result_peer_t2[(20*(i-1)+1):(20*i),] <- categorical_peer_list2[[i]][[8]]
  result_peer_t3[(20*(i-1)+1):(20*i),] <- categorical_peer_list2[[i]][[9]]
}

##############################################################################

## confidence interval
library(MultinomialCI)
truth <- c(0.1,0.2,0.3,0.25,0.15)
count <- rep(0, 4)
catet0_table <- matrix(NA, nrow = 500, ncol = 5)
catet1_table <- matrix(NA, nrow = 500, ncol = 5)
catet2_table <- matrix(NA, nrow = 500, ncol = 5)
catet3_table <- matrix(NA, nrow = 500, ncol = 5)



for(i in 1:np){
  catet0_table[i,] <- as.numeric(table(cate_peer_t0[i,]))
  catet1_table[i,] <- as.numeric(table(cate_peer_t1[i,]))
  catet2_table[i,] <- as.numeric(table(cate_peer_t2[i,]))
  catet3_table[i,] <- as.numeric(table(cate_peer_t3[i,]))
}



count <- rep(0,4)
for(i in 1:500){

  est0 <- multinomialCI(catet0_table[i,], 0.05)
  if(sum(est0[,1] <= truth)==5 & sum(truth <= est0[,2])==5){
    count[1] <- count[1] +1
  }
  est1 <- multinomialCI(catet1_table[i,], 0.05)
  if(sum(est1[,1] <= truth)==5 & sum(truth <= est1[,2])==5){
    count[2] <- count[2] +1
  }
  est2 <- multinomialCI(catet2_table[i,], 0.05)
  if(sum(est2[,1] <= truth)==5 & sum(truth <= est2[,2])==5){
    count[3] <- count[3] +1
  }
  est3 <- multinomialCI(catet3_table[i,], 0.05)
  if(sum(est3[,1] <= truth)==5 & sum(truth <= est3[,2])==5){
    count[4] <- count[4] +1
  }
}


#### standardize phi
stdphit0 <- rep(NA, np)
stdphit1 <- rep(NA, np)
stdphit2 <- rep(NA, np)
stdphit3 <- rep(NA, np)

## t = 0
for(i in 1:np){
  a <- m2.newphi(weight.matrix_peer[i,1], weight.matrix_peer[i,2], 
                 weight.matrix_peer[i,3],  cate_peer_t0[i,]) 
  b <- mean.newphi(weight.matrix_peer[i,1], cate_peer_t0[i,])
  stdphit0[i] <- result_peer_t0[i] - b
  stdphit0[i] <- stdphit0[i] / sqrt(a - b^2)
}


## t = 1
for(i in 1:np){
  a <- m2.newphi(weight.matrix_peer[i,1], weight.matrix_peer[i,2], 
                 weight.matrix_peer[i,3],  cate_peer_t1[i,]) 
  b <- mean.newphi(weight.matrix_peer[i,1], cate_peer_t1[i,])
  stdphit1[i] <- result_peer_t1[i] - b
  stdphit1[i] <- stdphit1[i] / sqrt(a - b^2)
}



## t = 2
for(i in 1:np){
  a <- m2.newphi(weight.matrix_peer[i,1], weight.matrix_peer[i,2], 
                 weight.matrix_peer[i,3],  cate_peer_t2[i,]) 
  b <- mean.newphi(weight.matrix_peer[i,1], cate_peer_t2[i,])
  stdphit2[i] <- result_peer_t2[i] - b
  stdphit2[i] <- stdphit2[i] / sqrt(a - b^2)
}



## t = 3
for(i in 1:np){
  a <- m2.newphi(weight.matrix_peer[i,1], weight.matrix_peer[i,2], 
                 weight.matrix_peer[i,3],  cate_peer_t3[i,]) 
  b <- mean.newphi(weight.matrix_peer[i,1], cate_peer_t3[i,])
  stdphit3[i] <- result_peer_t3[i] - b
  stdphit3[i] <- stdphit3[i] / sqrt(a - b^2)
}




## summary table
catepeer_summary <- matrix(NA, nrow = 4, ncol = 4)

rownames(catepeer_summary) <- c("t=0", "t=1", "t=2", "t=3")
colnames(catepeer_summary) <- c("95% CI coverage rate",  expression(phi), 
                                "pvalue(z)", "pvalue(permutation)")

catepeer_summary[,1] <- count / np

catepeer_summary[,2] <- c(mean(stdphit0), mean(stdphit1), mean(stdphit2), mean(stdphit3))

catepeer_summary[,3] <- c( sum( (1-pnorm(stdphit0)) <= 0.05 ) / 500, sum( (1-pnorm(stdphit1)) <= 0.05 ) / 500,
                           sum( (1-pnorm(stdphit2)) <= 0.05 ) / 500, sum( (1-pnorm(stdphit3)) <= 0.05 ) / 500)
catepeer_summary[,3] <- catepeer_summary[,3]*100

catepeer_summary[,4] <- c( sum( result_peer_t0[,2] <= 0.05 ) / 500, 
                           sum( result_peer_t1[,2] <= 0.05 ) / 500,
                           sum( result_peer_t2[,2] <= 0.05 ) / 500, 
                           sum( result_peer_t3[,2] <= 0.05 ) / 500)
catepeer_summary[,4] <- catepeer_summary[,4]*100


catepeer_summary <- as.data.frame(catepeer_summary)
print(xtable(catepeer_summary , digits = 2, row.names = TRUE))