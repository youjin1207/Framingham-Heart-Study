library(xtable)
source("makeCI.R")

alpha = 0.05 # type-I error

load("../data/LatentConti.RData")

outcome0 = matrix(0, 500, 200); outcome1 = matrix(0, 500, 200)
outcome2 = matrix(0, 500, 200); outcome3 = matrix(0, 500, 200)
moran0 = c(); moran1 = c(); moran2 = c(); moran3 = c()
Zpower = rep(0,4); Ppower = rep(0,4)
for(mm in 1:500){
  outcome0[mm,] = LatentConti[[mm]][[1]][1,]
  outcome1[mm,] = LatentConti[[mm]][[1]][2,]
  outcome2[mm,] = LatentConti[[mm]][[1]][3,]
  outcome3[mm,] = LatentConti[[mm]][[1]][4,]
  
  moran0[mm] = LatentConti[[mm]][[2]][1,1]
  moran1[mm] = LatentConti[[mm]][[2]][2,1]
  moran2[mm] = LatentConti[[mm]][[2]][3,1]
  moran3[mm] = LatentConti[[mm]][[2]][4,1]
  
  Zpower = Zpower + (LatentConti[[mm]][[2]][c(1:4), 2] <= alpha) / 500
  Ppower = Ppower + (LatentConti[[mm]][[2]][c(1:4), 3] <= alpha) / 500
}


# rho = 0.00
ci.r0 = ci_indep(outcome0)
right.index =  which(ci.r0[,1]*ci.r0[,2] < 0)
wrong.r0 = ci.r0[-right.index,]
right.r0 = ci.r0[right.index,]
cci.r0 = rbind(wrong.r0, right.r0)

# rho = 0.10
ci.r1 = ci_indep(outcome1)
right.index = which(ci.r1[,1]*ci.r1[,2] < 0)
wrong.r1 = ci.r1[-right.index,]
right.r1 = ci.r1[right.index,]
cci.r1 = rbind(wrong.r1, right.r1)

# rho = 0.15
ci.r2 = ci_indep(outcome2)
right.index = which(ci.r2[,1]*ci.r2[,2] < 0)
wrong.r2 = ci.r2[-right.index,]
right.r2 = ci.r2[right.index,]
cci.r2 = rbind(wrong.r2, right.r2)

# rho = 0.20
ci.r3 = ci_indep(outcome3)
right.index = which(ci.r3[,1]*ci.r3[,2] < 0)
wrong.r3 = ci.r3[-right.index,]
right.r3 = ci.r3[right.index,]
cci.r3 = rbind(wrong.r3, right.r3)


Ppower = round(Ppower*100)
coverages = c(cov_indep(outcome0, mu = 0), cov_indep(outcome1, mu = 0),
  cov_indep(outcome2, mu = 0), cov_indep(outcome3, mu = 0))
coverages = round(coverages*100)
pdf("../figures/coverage_latent200.pdf", width = 12, height = 8)
par(mfrow = c(1,4), oma = c(3, 4, 4, 2), cex.lab = 1.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5,
    mai = c(1.2, 0.5, 0.5, 0.5))
# t = 0
plot(x = c(cci.r0[1,1] , cci.r0[1,2]), y = c(1.0, 1.0), xlim = c(min(ci.r3) , max(ci.r3)),
     ylim = c(0,1.0), lwd = 1, type = "l",  
     xlab = paste("Coverage :", coverages[1], "%", "\n Reject independence :", Ppower[1], "%") , 
     ylab = "", col = "darkslateblue",
     mgp = c(7,1,0), xpd = FALSE,  xaxt = 'n')
axis(1, at=c(-0.3, 0, 0.3),  labels = c(-0.3, 0, 0.3), tck = 0.05)
for(i in 2:nrow(cci.r0)){
  lines(x = c(cci.r0[i,1], cci.r0[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.r0), 1.0 - (i-1)/nrow(cci.r0)),
        lwd = 1 , type = "l", col = "darkslateblue")
}
abline(h = cov_indep(outcome0, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)


plot(x = c(cci.r1[1,1] , cci.r1[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.r3) , max(ci.r3)),
     ylim = c(0,1.0), lwd = 1, type = "l",
     xlab = paste("Coverage :", coverages[2], "%", "\n Reject independence :", Ppower[2], "%"), 
     ylab = "", col = "lightpink",
     mgp = c(7,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-0.3, 0, 0.3),  labels = c(-0.3, 0, 0.3), tck = 0.05)
for(i in 2:nrow(cci.r1)){
  lines(x = c(cci.r1[i,1], cci.r1[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.r1), 1.0 - (i-1)/nrow(cci.r1)),
        lwd = 1 , type = "l", col = "lightpink")
}
abline(h = cov_indep(outcome1, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)


plot(x = c(cci.r2[1,1] , cci.r2[1,2]), y = c(1.0, 1.0), xlim = c(min(ci.r3) , max(ci.r3)),
     ylim = c(0,1.0), lwd = 1, type = "l",
     xlab = paste("Coverage :",coverages[3], "%", "\n Reject independence :", Ppower[3], "%"), 
     ylab="", col = "gold",
     mgp = c(7,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-0.3, 0, 0.3),  labels = c(-0.3, 0, 0.3), tck = 0.05)
for(i in 2:nrow(cci.r2)){
  lines(x = c(cci.r2[i,1], cci.r2[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.r2), 1.0 - (i-1)/nrow(cci.r2)),
        lwd = 1 , type = "l", col = "gold")
}
abline(h = cov_indep(outcome2, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)


plot(x = c(cci.r3[1,1] , cci.r3[1,2]), y = c(1.0, 1.0)   ,xlim = c(min(ci.r3) , max(ci.r3)),
     ylim = c(0,1.0), lwd = 1, type = "l", 
     xlab = paste("Coverage :", coverages[4], "%", "\n Reject independence :", Ppower[4], "%"), 
     ylab="", col = "skyblue",
     mgp = c(7,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-0.3, 0, 0.3),  labels = c(-0.3, 0, 0.3), tck = 0.05)
for(i in 2:nrow(cci.r3)){
  lines(x = c(cci.r3[i,1], cci.r3[i,2])
        , y = c(1.0 - (i-1)/nrow(cci.r3), 1.0 - (i-1)/nrow(cci.r3)),
        lwd = 1 , type = "l", col = "skyblue")
}
abline(h = cov_indep(outcome3, mu = 0), lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)
mtext(expression(paste("95% confidence intervals for ", mu, " assuming independence", sep="")), side = 3, line = 0, outer = TRUE, cex= 2.5, xpd = TRUE)
mtext("Proportion of Simulations", side = 2, line = 1, outer = TRUE, cex= 2.5, adj = 0.7)
dev.off()

#### make table
summary_latent_conti = matrix(NA, nrow = 4 , ncol = 3)
rownames(summary_latent_conti) = c("r = 0.00", "r = 0.20" ,"r = 0.30", "r = 0.40")
colnames(summary_latent_conti) =  c("95% CI coverage",
                                    "% of p-values(z)<=0.05", 
                                    "% of p-values(permuation) <=0.05")
summary_latent_conti[,1] = c(cov_indep(outcome0, mu = 0), cov_indep(outcome1, mu = 0),
                             cov_indep(outcome2, mu = 0), cov_indep(outcome3, mu = 0))

summary_latent_conti[,2] = Zpower*100
summary_latent_conti[,3] = Ppower*100
summary_latent_conti = as.data.frame(summary_latent_conti)
print(xtable(summary_latent_conti , digits = 2, row.names = TRUE))

