###########################################
library(igraph)
library(MASS)
library(mvrtn)
library(MultinomialCI)
library(xtable)
###########################################
load("data/functions.RData")
load("data/continuous_peer.RData")
###########################################
np <- 500
popn <- 1000

moran.t0 <- matrix(NA, ncol = 5, nrow = np)
moran.t1 <- matrix(NA, ncol = 5, nrow = np)
moran.t2 <- matrix(NA, ncol = 5, nrow = np)
moran.t3 <- matrix(NA, ncol = 5, nrow = np)
moran.t4 <- matrix(NA, ncol = 5, nrow = np)
moran.t5 <- matrix(NA, ncol = 5, nrow = np)

conti.t0.outcomes <- matrix(NA, ncol = popn, nrow = np)
conti.t1.outcomes <- matrix(NA, ncol = popn, nrow = np)
conti.t2.outcomes <- matrix(NA, ncol = popn, nrow = np)
conti.t3.outcomes <- matrix(NA, ncol = popn, nrow = np)
conti.t4.outcomes <- matrix(NA, ncol = popn, nrow = np)
conti.t5.outcomes <- matrix(NA, ncol = popn, nrow = np)

for(i in 1:20){
  conti.t0.outcomes[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[7]]
  conti.t1.outcomes[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[8]]
  conti.t2.outcomes[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[9]]
  conti.t3.outcomes[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[10]]
  conti.t4.outcomes[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[11]]
  conti.t5.outcomes[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[12]]
  moran.t0[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[1]]
  moran.t1[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[2]]
  moran.t2[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[3]]
  moran.t3[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[4]]
  moran.t4[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[5]]
  moran.t5[(25*(i-1) + 1) :(25*i) ,] <- continuous_peer_list[[i]][[6]]
}

####### coverage rate ########
# time = 0
ci.t0 <- ci_indep(conti.t0.outcomes)
right.index <- which(ci.t0[,1]*ci.t0[,2] < 0)
wrong.t0 <- ci.t0[-right.index,]
right.t0 <- ci.t0[right.index,]
cci.t0 <- rbind(wrong.t0, right.t0)

# time = 1
ci.t1 <- ci_indep(conti.t1.outcomes)
right.index <- which(ci.t1[,1]*ci.t1[,2] < 0)
wrong.t1 <- ci.t1[-right.index,]
right.t1 <- ci.t1[right.index,]
cci.t1 <- rbind(wrong.t1, right.t1)

# time = 2
ci.t2 <- ci_indep(conti.t2.outcomes)
right.index <- which(ci.t2[,1]*ci.t2[,2] < 0)
wrong.t2 <- ci.t2[-right.index,]
right.t2 <- ci.t2[right.index,]
cci.t2 <- rbind(wrong.t2, right.t2)

# time = 3
ci.t3 <- ci_indep(conti.t3.outcomes)
right.index <- which(ci.t3[,1]*ci.t3[,2] < 0)
wrong.t3 <- ci.t3[-right.index,]
right.t3 <- ci.t3[right.index,]
cci.t3 <- rbind(wrong.t3, right.t3)

# time = 4
ci.t4 <- ci_indep(conti.t4.outcomes)
right.index <- which(ci.t4[,1]*ci.t4[,2] < 0)
wrong.t4 <- ci.t4[-right.index,]
right.t4 <- ci.t4[right.index,]
cci.t4 <- rbind(wrong.t4, right.t4)

# time = 5
ci.t5 <- ci_indep(conti.t5.outcomes)
right.index <- which(ci.t5[,1]*ci.t5[,2] < 0)
wrong.t5 <- ci.t5[-right.index,]
right.t5 <- ci.t5[right.index,]
cci.t5 <- rbind(wrong.t5, right.t5)

pdf("figures/coverage_peer.pdf")
par(mfrow = c(1,4), mar =  c(10, 4, 4, 2))
# t = 0
plot(x = c(cci.t0[1,1] , cci.t0[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.t0) , max(ci.t0)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "time = 0",
     xlab = paste("Coverage :", cov_indep(conti.t0.outcomes, mu = 0) ) , ylab = "Proportion of Simulations", col = "darkslateblue")
for(i in 2:nrow(cci.t0)){
  lines(x = c(cci.t0[i,1], cci.t0[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t0), 1.0 - (i-1)/nrow(cci.t0)),
        lwd = 1 , type = "l", col = "darkslateblue")
}
abline(h = cov_indep(conti.t0.outcomes, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# t = 1
plot(x = c(cci.t1[1,1] , cci.t1[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.t1) , max(ci.t1)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "time = 1",
     xlab = paste("Coverage :", cov_indep(conti.t1.outcomes, mu = 0) ) , ylab = "Proportion of Simulations", col = "lightpink")
for(i in 2:nrow(cci.t1)){
  lines(x = c(cci.t1[i,1], cci.t1[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t1), 1.0 - (i-1)/nrow(cci.t1)),
        lwd = 1 , type = "l", col = "lightpink")
}
abline(h = cov_indep(conti.t1.outcomes, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# t = 2
plot(x = c(cci.t2[1,1] , cci.t2[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.t2) , max(ci.t2)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "time = 2",
     xlab = paste("Coverage :", cov_indep(conti.t2.outcomes, mu = 0) ) , ylab = "Proportion of Simulations", col = "grey")
for(i in 2:nrow(cci.t2)){
  lines(x = c(cci.t2[i,1], cci.t2[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t2), 1.0 - (i-1)/nrow(cci.t2)),
        lwd = 1 , type = "l", col = "gold")
}
abline(h = cov_indep(conti.t2.outcomes, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)


# t = 3
plot(x = c(cci.t3[1,1] , cci.t3[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.t3) , max(ci.t3)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "time = 3",
     xlab = paste("Coverage :", cov_indep(conti.t3.outcomes, mu = 0) ) , ylab = "Proportion of Simulations", col = "skyblue")
for(i in 2:nrow(cci.t3)){
  lines(x = c(cci.t3[i,1], cci.t3[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t3), 1.0 - (i-1)/nrow(cci.t3)),
        lwd = 1 , type = "l", col = "skyblue")
}
abline(h = cov_indep(conti.t3.outcomes, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

dev.off()





