library(xtable)
source("makeCI.R")

alpha = 0.05 # type-I error

load("../Data/PeerConti.RData")

outcome0 = matrix(0, 500, 1000); outcome1 = matrix(0, 500, 1000)
outcome2 = matrix(0, 500, 1000); outcome3 = matrix(0, 500, 1000)
moran0 = c(); moran1 = c(); moran2 = c(); moran3 = c()
Zpower = rep(0,4); Ppower = rep(0,4)
for(mm in 1:500){
  outcome0[mm,] = PeerConti[[mm]][[1]][1,]
  outcome1[mm,] = PeerConti[[mm]][[1]][2,]
  outcome2[mm,] = PeerConti[[mm]][[1]][3,]
  outcome3[mm,] = PeerConti[[mm]][[1]][4,]
  
  moran0[mm] = PeerConti[[mm]][[2]][1,1]
  moran1[mm] = PeerConti[[mm]][[2]][2,1]
  moran2[mm] = PeerConti[[mm]][[2]][3,1]
  moran3[mm] = PeerConti[[mm]][[2]][4,1]
  
  Zpower = Zpower + (PeerConti[[mm]][[2]][c(1:4), 2] <= alpha) / 500
  Ppower = Ppower + (PeerConti[[mm]][[2]][c(1:4), 3] <= alpha) / 500
}

####### coverage rate ########
# time = 0
ci.t0 = ci_indep(outcome0)
right.index = which(ci.t0[,1]*ci.t0[,2] < 0)
wrong.t0 = ci.t0[-right.index,]
right.t0 = ci.t0[right.index,]
cci.t0 = rbind(wrong.t0, right.t0)

# time = 1
ci.t1 = ci_indep(outcome1)
right.index = which(ci.t1[,1]*ci.t1[,2] < 0)
wrong.t1 = ci.t1[-right.index,]
right.t1 = ci.t1[right.index,]
cci.t1 = rbind(wrong.t1, right.t1)

# time = 2
ci.t2 = ci_indep(outcome2)
right.index = which(ci.t2[,1]*ci.t2[,2] < 0)
wrong.t2 = ci.t2[-right.index,]
right.t2 = ci.t2[right.index,]
cci.t2 = rbind(wrong.t2, right.t2)

# time = 3
ci.t3 = ci_indep(outcome3)
right.index = which(ci.t3[,1]*ci.t3[,2] < 0)
wrong.t3 = ci.t3[-right.index,]
right.t3 = ci.t3[right.index,]
cci.t3 = rbind(wrong.t3, right.t3)


pdf("../Figure/coverage_peer.pdf", width = 12, height = 8)
par(mfrow = c(1,4), oma = c(5, 5, 2, 2), cex.lab = 2, 
    cex.main = 3, cex.axis = 2, tcl = 0.5,
    mai = c(0.7, 0.3, 0.3, 0.3))
# t = 0
plot(x = c(cci.t0[1,1] , cci.t0[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.t1) , max(ci.t0)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "t = 0", 
     xlab = paste("Coverage :", cov_indep(outcome0, mu = 0) ) , ylab = "", col = "darkslateblue")
for(i in 2:nrow(cci.t0)){
  lines(x = c(cci.t0[i,1], cci.t0[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t0), 1.0 - (i-1)/nrow(cci.t0)),
        lwd = 1 , type = "l", col = "darkslateblue")
}
abline(h = cov_indep(outcome0, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# t = 1
plot(x = c(cci.t1[1,1] , cci.t1[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.t1) , max(ci.t0)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "t = 1",
     xlab = paste("Coverage :", cov_indep(outcome1, mu = 0) ), ylab = "", col = "lightpink")
for(i in 2:nrow(cci.t1)){
  lines(x = c(cci.t1[i,1], cci.t1[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t1), 1.0 - (i-1)/nrow(cci.t1)),
        lwd = 1 , type = "l", col = "lightpink")
}
abline(h = cov_indep(outcome1, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# t = 2
plot(x = c(cci.t2[1,1] , cci.t2[1,2]), y = c(1.0, 1.0), xlim = c(min(ci.t1) , max(ci.t0)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "t = 2",
     xlab = paste("Coverage :", cov_indep(outcome2, mu = 0) ), ylab="", col = "gold")
for(i in 2:nrow(cci.t2)){
  lines(x = c(cci.t2[i,1], cci.t2[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t2), 1.0 - (i-1)/nrow(cci.t2)),
        lwd = 1 , type = "l", col = "gold")
}
abline(h = cov_indep(outcome2, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# t = 3
plot(x = c(cci.t3[1,1] , cci.t3[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.t1) , max(ci.t0)),
     ylim = c(0,1.0), lwd = 1, type = "l", main = "t = 3",
     xlab = paste("Coverage :", cov_indep(outcome3, mu = 0) ), ylab="", col = "skyblue")
for(i in 2:nrow(cci.t3)){
  lines(x = c(cci.t3[i,1], cci.t3[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t3), 1.0 - (i-1)/nrow(cci.t3)),
        lwd = 1 , type = "l", col = "skyblue")
}
abline(h = cov_indep(outcome3, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)
mtext("95% CI coverage assuming independence", side = 1, line = 1.5, outer = TRUE, cex= 2.5)
mtext("Proportion of Simulations", side = 2, line = 1, outer = TRUE, cex= 2.5)
dev.off()

#### make table
summary_peer_conti <- matrix(NA, nrow = 4 , ncol = 3)
rownames(summary_peer_conti) <- c("t = 0", "t = 1" ,"t = 2", "t = 3")
colnames(summary_peer_conti) <-  c("95% CI coverage",
                                    "% of p-values(z)<=0.05", 
                                    "% of p-values(permuation) <=0.05")
summary_peer_conti[,1] <- c(cov_indep(outcome0, mu = 0), cov_indep(outcome1, mu = 0),
                              cov_indep(outcome2, mu = 0), cov_indep(outcome3, mu = 0))

summary_peer_conti[,2] <- Zpower*100
summary_peer_conti[,3] <- Ppower*100
summary_peer_conti = as.data.frame(summary_peer_conti)
print(xtable(summary_peer_conti , digits = 2, row.names = TRUE))


