###########################################
library(igraph)
library(MASS)
library(mvrtn)
library(MultinomialCI)
library(xtable)
###########################################
load("data/catelatentr0.RData")
load("data/catelatentr15.RData")
load("data/catelatentr2.RData")
load("data/catelatentr3.RData")
########################################
popn <- 1000
np <- 500
# weight matrix
w.matrix0 <- matrix(NA, ncol = 3, nrow = np)
w.matrix1 <- matrix(NA, ncol = 3, nrow = np)
w.matrix2 <- matrix(NA, ncol = 3, nrow = np)
w.matrix3 <- matrix(NA, ncol = 3, nrow = np)


# outcomes
cate_latent_r0 <- matrix(NA, ncol = popn, nrow = np)
cate_latent_r1 <- matrix(NA, ncol = popn, nrow = np)
cate_latent_r2 <- matrix(NA, ncol = popn, nrow = np)
cate_latent_r3 <- matrix(NA, ncol = popn, nrow = np)

# test statistics and their pvalues
result_latent0 <- matrix(NA, ncol = 2, nrow = np )
result_latent1 <- matrix(NA, ncol = 2, nrow = np )
result_latent2 <- matrix(NA, ncol = 2, nrow = np )
result_latent3 <- matrix(NA, ncol = 2, nrow = np )


for(i in 1:25){
  # weight matrix
  w.matrix0[(20*(i-1)+1):(20*i),] <- catelatentr0_list[[i]][[1]]
  w.matrix1[(20*(i-1)+1):(20*i),] <- catelatentr15_list[[i]][[1]]
  w.matrix2[(20*(i-1)+1):(20*i),] <- catelatentr2_list[[i]][[1]]
  w.matrix3[(20*(i-1)+1):(20*i),] <- catelatentr3_list[[i]][[1]]
  
  # outcome matrix
  cate_latent_r0[(20*(i-1)+1):(20*i),] <- catelatentr0_list[[i]][[2]]
  cate_latent_r1[(20*(i-1)+1):(20*i),] <- catelatentr15_list[[i]][[2]]
  cate_latent_r2[(20*(i-1)+1):(20*i),] <- catelatentr2_list[[i]][[2]]
  cate_latent_r3[(20*(i-1)+1):(20*i),] <- catelatentr3_list[[i]][[2]]
  
  # test statistics and their pvalues
  result_latent0[(20*(i-1)+1):(20*i),] <- catelatentr0_list[[i]][[3]]
  result_latent1[(20*(i-1)+1):(20*i),] <- catelatentr15_list[[i]][[3]]
  result_latent2[(20*(i-1)+1):(20*i),] <- catelatentr2_list[[i]][[3]]
  result_latent3[(20*(i-1)+1):(20*i),] <- catelatentr3_list[[i]][[3]]
}


##############################################################################
## confidence interval

truth <- c(0.1,0.2,0.3,0.25,0.15)
count <- rep(0, 4)
cater0_table <- matrix(NA, nrow = 500, ncol = 5)
cater1_table <- matrix(NA, nrow = 500, ncol = 5)
cater2_table <- matrix(NA, nrow = 500, ncol = 5)
cater3_table <- matrix(NA, nrow = 500, ncol = 5)
np <- 500

for(i in 1:np){
  cater0_table[i,] <- as.numeric(table(cate_latent_r0[i,]))
  cater1_table[i,] <- as.numeric(table(cate_latent_r1[i,]))
  cater2_table[i,] <- as.numeric(table(cate_latent_r2[i,]))
  cater3_table[i,] <- as.numeric(table(cate_latent_r3[i,]))
}

count <- rep(0,4)
for(i in 1:500){
  est0 <- multinomialCI(cater0_table[i,], 0.05)
  if(sum(est0[,1] <= truth)==5 & sum(truth <= est0[,2])==5){
    count[1] <- count[1] +1
  }
  est1 <- multinomialCI(cater1_table[i,], 0.05)
  if(sum(est1[,1] <= truth)==5 & sum(truth <= est1[,2])==5){
    count[2] <- count[2] +1
  }
  
  est2 <- multinomialCI(cater2_table[i,], 0.05)
  if(sum(est2[,1] <= truth)==5 & sum(truth <= est2[,2])==5){
    count[3] <- count[3] +1
  }
  
  est3 <- multinomialCI(cater3_table[i,], 0.05)
  if(sum(est3[,1] <= truth)==5 & sum(truth <= est3[,2])==5){
    count[4] <- count[4] +1
  }
}


#### standardize phi
stdphir0 <- rep(NA, np)
stdphir1 <- rep(NA, np)
stdphir2 <- rep(NA, np)
stdphir3 <- rep(NA, np)

## r = 0.0
for(i in 1:np){
  a <- m2.newphi(w.matrix0[i,1], w.matrix0[i,2], 
                 w.matrix0[i,3],  cate_latent_r0[i,]) 
  b <- mean.newphi(w.matrix0[i,1], cate_latent_r0[i,])
  stdphir0[i] <- result_latent0[i,1] - b
  stdphir0[i] <- stdphir0[i] / sqrt(a - b^2)
}

## r = 0.1
for(i in 1:np){
  a <- m2.newphi(w.matrix1[i,1], w.matrix1[i,2], 
                 w.matrix1[i,3],  cate_latent_r1[i,]) 
  b <- mean.newphi(w.matrix1[i,1], cate_latent_r1[i,])
  stdphir1[i] <- result_latent1[i,1] - b
  stdphir1[i] <- stdphir1[i] / sqrt(a - b^2)
}


## r = 0.2
for(i in 1:np){
  a <- m2.newphi(w.matrix2[i,1], w.matrix2[i,2], 
                 w.matrix2[i,3],  cate_latent_r2[i,]) 
  b <- mean.newphi(w.matrix2[i,1], cate_latent_r2[i,])
  stdphir2[i] <- result_latent2[i,1] - b
  stdphir2[i] <- stdphir2[i] / sqrt(a - b^2)
}


## r = 0.3
for(i in 1:np){
  a <- m2.newphi(w.matrix3[i,1], w.matrix3[i,2], 
                 w.matrix3[i,3],  cate_latent_r3[i,]) 
  b <- mean.newphi(w.matrix3[i,1], cate_latent_r3[i,])
  stdphir3[i] <- result_latent3[i,1] - b
  stdphir3[i] <- stdphir3[i] / sqrt(a - b^2)
}


## summary table
catelatent_summary <- matrix(NA, nrow = 4, ncol = 4)

rownames(catelatent_summary) <- c("r=0.0", "r=0.15", "r=0.2", "r=0.3")
colnames(catelatent_summary) <- c("95% CI coverage rate",  expression(phi), 
                                  "pvalue(z)" , "pvalue(permutation)")

catelatent_summary[,1] <- count / np

catelatent_summary[,2] <- c(mean(stdphir0), mean(stdphir1), mean(stdphir2), mean(stdphir3))

catelatent_summary[,3] <- c( sum( (1-pnorm(stdphir0)) <= 0.05 ) / 500, sum( (1-pnorm(stdphir1)) <= 0.05 ) / 500,
                             sum( (1-pnorm(stdphir2)) <= 0.05 ) / 500, sum( (1-pnorm(stdphir3)) <= 0.05 ) / 500)
catelatent_summary[,3] <- catelatent_summary[,3]*100

catelatent_summary[,4] <- c( sum( result_latent0[,2] <= 0.05 ) / 500, sum( result_latent1[,2] <= 0.05 ) / 500,
                             sum( result_latent2[,2] <= 0.05 ) / 500, sum( result_latent3[,2] <= 0.05 ) / 500)
catelatent_summary[,4] <- catelatent_summary[,4]*100

catelatent_summary <- as.data.frame(catelatent_summary)
print(xtable(catelatent_summary , digits = 2, row.names = TRUE))

