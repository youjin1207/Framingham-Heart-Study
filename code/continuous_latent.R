###########################################
library(igraph)
library(MASS)
library(mvrtn)
library(MultinomialCI)
library(xtable)
###########################################
load("data/functions.RData")
load("data/conti_latent_r0_full.RData")
load("data/conti_latent_r07_full.RData")
load("data/conti_latent_r1_full.RData")
load("data/conti_latent_r15_full.RData")

###########################################
popn <- 1000
np <- 500
## For CI coverage
conti.r0.outcomes <- matrix(NA, ncol = popn, nrow = np)
conti.r1.outcomes <- matrix(NA, ncol = popn, nrow = np)
conti.r15.outcomes <- matrix(NA, ncol = popn, nrow = np)
conti.r07.outcomes <- matrix(NA, ncol = popn, nrow = np)
## Moran's I
moran.r0 <- matrix(NA, ncol = 5, nrow = np)
moran.r07 <- matrix(NA, ncol = 5, nrow = np)
moran.r1 <- matrix(NA, ncol = 5, nrow = np)
moran.r15 <- matrix(NA, ncol = 5, nrow = np)


for(i in 1:25){
  conti.r0.outcomes[(20*(i-1)+1):(20*i),] <- result_list_r0_full[[i]][[1]]
  moran.r0[(20*(i-1)+1):(20*i), ] <- result_list_r0_full[[i]][[2]]
  
  conti.r07.outcomes[(20*(i-1)+1):(20*i),] <- result_list_r07_full[[i]][[1]]
  moran.r07[(20*(i-1)+1):(20*i), ] <- result_list_r07_full[[i]][[2]]

  conti.r1.outcomes[(20*(i-1)+1):(20*i),] <- result_list_r1_full[[i]][[1]]
  moran.r1[(20*(i-1)+1):(20*i), ] <- result_list_r1_full[[i]][[2]]

  conti.r15.outcomes[(20*(i-1)+1):(20*i),] <- result_list_r15_full[[i]][[1]]
  moran.r15[(20*(i-1)+1):(20*i), ] <- result_list_r15_full[[i]][[2]]
}


# rho = 0.0
ci.r0 <- ci_indep(conti.r0.outcomes)
right.index <- which(ci.r0[,1]*ci.r0[,2] < 0)
wrong.r0 <- ci.r0[-right.index,]
right.r0 <- ci.r0[right.index,]
cci.r0 <- rbind(wrong.r0, right.r0)

# rho = 0.07
ci.r1 <- ci_indep(conti.r07.outcomes)
right.index <- which(ci.r1[,1]*ci.r1[,2] < 0)
wrong.r1 <- ci.r1[-right.index,]
right.r1 <- ci.r1[right.index,]
cci.r1 <- rbind(wrong.r1, right.r1)

# rho = 0.10
ci.r2 <- ci_indep(conti.r1.outcomes)
right.index <- which(ci.r2[,1]*ci.r2[,2] < 0)
wrong.r2 <- ci.r2[-right.index,]
right.r2 <- ci.r2[right.index,]
cci.r2 <- rbind(wrong.r2, right.r2)

# rho = 0.15
ci.r3 <- ci_indep(conti.r15.outcomes)
right.index <- which(ci.r3[,1]*ci.r3[,2] < 0)
wrong.r3 <- ci.r3[-right.index,]
right.r3 <- ci.r3[right.index,]
cci.r3 <- rbind(wrong.r3, right.r3)


pdf("figures/coverage_latent.pdf")
par(mfrow = c(1,4), mar =  c(10, 4, 4, 2))
# rho = 0.0
plot(x = c(cci.r0[1,1] , cci.r0[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.r0) , max(ci.r0)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "rho = 0.00",
     xlab = paste("Coverage :", cov_indep(conti.r0.outcomes, mu = 0) ) , 
     ylab = "Proportion of Simulations", col = "darkslateblue")
for(i in 2:nrow(cci.r0)){
  lines(x = c(cci.r0[i,1], cci.r0[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.r0), 1.0 - (i-1)/nrow(cci.r0)),
        lwd = 1 , type = "l", col = "darkslateblue")
}
abline(h = cov_indep(conti.r0.outcomes, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# rho = 0.07
plot(x = c(cci.r1[1,1] , cci.r1[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.r1) , max(ci.r1)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "rho = 0.07",
     xlab = paste("Coverage :", cov_indep(conti.r07.outcomes, mu = 0) ) , 
     ylab = "Proportion of Simulations", col = "lightpink")
for(i in 2:nrow(cci.r1)){
  lines(x = c(cci.r1[i,1], cci.r1[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.r1), 1.0 - (i-1)/nrow(cci.r1)),
        lwd = 1 , type = "l", col = "lightpink")
}
abline(h = cov_indep(conti.r07.outcomes, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# rho = 0.10
plot(x = c(cci.r2[1,1] , cci.r2[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.r2) , max(ci.r2)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "rho = 0.10",
     xlab = paste("Coverage :", cov_indep(conti.r1.outcomes, mu = 0) ) , 
     ylab = "Proportion of Simulations", col = "gold")
for(i in 2:nrow(cci.r2)){
  lines(x = c(cci.r2[i,1], cci.r2[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.r2), 1.0 - (i-1)/nrow(cci.r2)),
        lwd = 1 , type = "l", col = "gold")
}
abline(h = cov_indep(conti.r1.outcomes, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# rho = 0.15
plot(x = c(cci.r3[1,1] , cci.r3[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.r3) , max(ci.r3)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "rho = 0.15",
     xlab = paste("Coverage :", cov_indep(conti.r15.outcomes, mu = 0) ) , 
     ylab = "Proportion of Simulations", col = "skyblue")
for(i in 2:nrow(cci.r3)){
  lines(x = c(cci.r3[i,1], cci.r3[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.r3), 1.0 - (i-1)/nrow(cci.r3)),
        lwd = 1 , type = "l", col = "skyblue")
}
abline(h = cov_indep(conti.r15.outcomes, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)
dev.off()

########### results of permutation test #################
coverage_contipeer <- c()
coverage_contipeer <- c(cov_indep(conti.r0.outcomes, mu = 0), cov_indep(conti.r07.outcomes, mu = 0),
                        cov_indep(conti.r1.outcomes, mu = 0), cov_indep(conti.r15.outcomes, mu = 0))

avgmoran <- c()
avgmoran <- c(mean(moran.r0[,1]), mean(moran.r07[,1]),
              mean(moran.r1[,1]), mean(moran.r15[,1]))


## make a summary table
summary_latent_conti <- matrix(NA, nrow = 4 , ncol = 4)
rownames(summary_latent_conti) <- c("r=0.00", "r=0.07" ,"r=0.10", "r=0.15")
colnames(summary_latent_conti) <-  c("95% CI coverage",
                                     "avg.Moran's I","% of p-values(z)<=0.05", 
                                     "% of p-values(permuation) <=0.05")
summary_latent_conti[,1] <- c(cov_indep(conti.r0.outcomes, mu = 0), cov_indep(conti.r07.outcomes, mu = 0),
                              cov_indep(conti.r1.outcomes, mu = 0), cov_indep(conti.r15.outcomes, mu = 0))

summary_latent_conti[,2] <- avgmoran

summary_latent_conti[,3] <- c(sum(moran.r0[,2] <= 0.05), sum(moran.r07[,2] <= 0.05),
                              sum(moran.r1[,2] <= 0.05), sum(moran.r15[,2] <= 0.05))
summary_latent_conti[,3] <- summary_latent_conti[,3] / 500
summary_latent_conti[,3] <- summary_latent_conti[,3]*100
summary_latent_conti[,4] <- c(sum(moran.r0[,3] <= 0.05), sum(moran.r07[,3] <= 0.05),
                              sum(moran.r1[,3] <= 0.05), sum(moran.r15[,3] <= 0.05))
summary_latent_conti[,4] <- summary_latent_conti[,4] / 500
summary_latent_conti[,4] <- summary_latent_conti[,4]*100
summary_latent_conti <- as.data.frame(summary_latent_conti)
print(xtable(summary_latent_conti , digits = 2, row.names = TRUE))


