library(xtable)
library(mvrtn)
library(MultinomialCI)
source("makeCI.R")

alpha = 0.05 # type-I error

load("../data/PeerConti.RData")

outcome0 = matrix(0, 500, 200); outcome1 = matrix(0, 500, 200)
outcome2 = matrix(0, 500, 200); outcome3 = matrix(0, 500, 200)
moran0 = c(); moran1 = c(); moran2 = c(); moran3 = c()
Zpower = rep(0,4); Ppower = rep(0,4)
for(mm in 1:500){
  outcome0[mm,] = PeerConti200[[mm]][[1]][1,]
  outcome1[mm,] = PeerConti200[[mm]][[1]][2,]
  outcome2[mm,] = PeerConti200[[mm]][[1]][3,]
  outcome3[mm,] = PeerConti200[[mm]][[1]][4,]
  
  moran0[mm] = PeerConti200[[mm]][[2]][1,1]
  moran1[mm] = PeerConti200[[mm]][[2]][2,1]
  moran2[mm] = PeerConti200[[mm]][[2]][3,1]
  moran3[mm] = PeerConti200[[mm]][[2]][4,1]
  
  Zpower = Zpower + (PeerConti200[[mm]][[2]][c(1:4), 2] <= alpha) / 500
  Ppower = Ppower + (PeerConti200[[mm]][[2]][c(1:4), 3] <= alpha) / 500
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

## figure
Ppower = round(Ppower*100)
coverages = c(cov_indep(outcome0, mu = 0), cov_indep(outcome1, mu = 0),
  cov_indep(outcome2, mu = 0), cov_indep(outcome3, mu = 0))
coverages = round(coverages*100)
pdf("../figures/coverage_peer200.pdf", width = 12, height = 8)
par(mfrow = c(1,4), oma = c(3, 4, 4, 2), cex.lab = 1.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5,
    mai = c(1.2, 0.5, 0.5, 0.5))
# t = 0
plot(x = c(cci.t0[1,1] , cci.t0[1,2]), y = c(1.0, 1.0), xlim = c(min(ci.t1) , max(ci.t0)),
     ylim = c(0,1.0), lwd = 1, type = "l",  
     xlab = paste("Coverage :", coverages[1], "%", "\n Reject independence :", Ppower[1], "%") , 
     ylab = "", col = "darkslateblue",
     mgp = c(7,1,0), xpd = FALSE,  xaxt = 'n')
axis(1, at=c(-0.3, 0, 0.3),  labels = c(-0.3, 0, 0.3), tck = 0.05)
for(i in 2:nrow(cci.t0)){
  lines(x = c(cci.t0[i,1], cci.t0[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t0), 1.0 - (i-1)/nrow(cci.t0)),
        lwd = 1 , type = "l", col = "darkslateblue")
}
abline(h = cov_indep(outcome0, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# t = 1
plot(x = c(cci.t1[1,1] , cci.t1[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.t1) , max(ci.t0)),
     ylim = c(0,1.0), lwd = 1, type = "l",
     xlab = paste("Coverage :", coverages[2], "%", "\n Reject independence :", Ppower[2], "%"), 
     ylab = "", col = "lightpink",
     mgp = c(7,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-0.3, 0, 0.3),  labels = c(-0.3, 0, 0.3), tck = 0.05)
for(i in 2:nrow(cci.t1)){
  lines(x = c(cci.t1[i,1], cci.t1[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t1), 1.0 - (i-1)/nrow(cci.t1)),
        lwd = 1 , type = "l", col = "lightpink")
}
abline(h = cov_indep(outcome1, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# t = 2
plot(x = c(cci.t2[1,1] , cci.t2[1,2]), y = c(1.0, 1.0), xlim = c(min(ci.t1) , max(ci.t0)),
     ylim = c(0,1.0), lwd = 1, type = "l",
     xlab = paste("Coverage :",coverages[3], "%", "\n Reject independence :", Ppower[3], "%"), 
     ylab="", col = "gold",
     mgp = c(7,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-0.3, 0, 0.3),  labels = c(-0.3, 0, 0.3), tck = 0.05)
for(i in 2:nrow(cci.t2)){
  lines(x = c(cci.t2[i,1], cci.t2[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t2), 1.0 - (i-1)/nrow(cci.t2)),
        lwd = 1 , type = "l", col = "gold")
}
abline(h = cov_indep(outcome2, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# t = 3
plot(x = c(cci.t3[1,1] , cci.t3[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.t1) , max(ci.t0)),
     ylim = c(0,1.0), lwd = 1, type = "l", 
     xlab = paste("Coverage :", coverages[4], "%", "\n Reject independence :", Ppower[4], "%"), 
     ylab="", col = "skyblue",
     mgp = c(7,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-0.3, 0, 0.3),  labels = c(-0.3, 0, 0.3), tck = 0.05)
for(i in 2:nrow(cci.t3)){
  lines(x = c(cci.t3[i,1], cci.t3[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.t3), 1.0 - (i-1)/nrow(cci.t3)),
        lwd = 1 , type = "l", col = "skyblue")
}
abline(h = cov_indep(outcome3, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)
mtext(expression(paste("95% confidence intervals for ", mu, " assuming independence", sep="")), side = 3, line = 0, outer = TRUE, cex= 2.5, xpd = TRUE)
mtext("Proportion of Simulations", side = 2, line = 1, outer = TRUE, cex= 2.5, adj = 0.7)
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


## bias & se table
average.bias = c(mean(apply(outcome0, 1, mean)), mean(apply(outcome1, 1, mean)),
                 mean(apply(outcome2, 1, mean)), mean(apply(outcome3, 1, mean)))
average.absbias = c(mean(abs(apply(outcome0, 1, mean))), mean(abs(apply(outcome1, 1, mean))),
                    mean(abs(apply(outcome2, 1, mean))), mean(abs(apply(outcome3, 1, mean))))
mean.se = c(mean(apply(outcome0, 1, sd)/sqrt(200)),
            mean(apply(outcome1, 1, sd)/sqrt(200)),
            mean(apply(outcome2, 1, sd)/sqrt(200)),
            mean(apply(outcome3, 1, sd)/sqrt(200)))
se.replicate = c(sd(apply(outcome0, 1, mean)), sd(apply(outcome1, 1, mean)),
                 sd(apply(outcome2, 1, mean)), sd(apply(outcome3, 1, mean)))

mat = rbind(average.bias, average.absbias, mean.se, se.replicate)
xtable(mat, digits = 3)

