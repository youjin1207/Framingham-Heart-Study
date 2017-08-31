library(xtable)
library(MultinomialCI)

source("makeCI.R")

alpha = 0.05 # type-I error
np = 500 # number of replicates

load("../data/PeerCate.RData")
outcome0 = matrix(0, 500, 200); outcome1 = matrix(0, 500, 200)
outcome2 = matrix(0, 500, 200); outcome3 = matrix(0, 500, 200)
moran0 = c(); moran1 = c(); moran2 = c(); moran3 = c()
Zpower = rep(0,4); Ppower = rep(0,4)
for(mm in 1:500){
  outcome0[mm,] = PeerCate[[mm]][[1]][1,]
  outcome1[mm,] = PeerCate[[mm]][[1]][2,]
  outcome2[mm,] = PeerCate[[mm]][[1]][3,]
  outcome3[mm,] = PeerCate[[mm]][[1]][4,]
  
  phi0[mm] = PeerCate[[mm]][[2]][1,1]
  phi1[mm] = PeerCate[[mm]][[2]][2,1]
  phi2[mm] = PeerCate[[mm]][[2]][3,1]
  phi3[mm] = PeerCate[[mm]][[2]][4,1]
  
  Zpower = Zpower + (PeerCaqe[[mm]][[2]][c(1:4), 2] <= alpha) / 500
  Ppower = Ppower + (PeerCate[[mm]][[2]][c(1:4), 3] <= alpha) / 500
}


truth <- c(0.1, 0.2, 0.3, 0.25, 0.15)
count <- rep(0, 4)
catet0_table <- matrix(NA, nrow = 500, ncol = 5)
catet1_table <- matrix(NA, nrow = 500, ncol = 5)
catet2_table <- matrix(NA, nrow = 500, ncol = 5)
catet3_table <- matrix(NA, nrow = 500, ncol = 5)

for(i in 1:np){
  catet0_table[i,] <- as.numeric(table(outcome0[i,]))
  catet1_table[i,] <- as.numeric(table(outcome1[i,]))
  catet2_table[i,] <- as.numeric(table(outcome2[i,]))
  catet3_table[i,] <- as.numeric(table(outcome3[i,]))
}


count = rep(0,4)
for(i in 1:500){
  est0 = multinomialCI(catet0_table[i,], 0.05)
  if(sum(est0[,1] <= truth)==5 & sum(truth <= est0[,2])==5){
    count[1] = count[1] +1
  }
  est1 = multinomialCI(catet1_table[i,], 0.05)
  if(sum(est1[,1] <= truth)==5 & sum(truth <= est1[,2])==5){
    count[2] = count[2] +1
  }
  est2 = multinomialCI(catet2_table[i,], 0.05)
  if(sum(est2[,1] <= truth)==5 & sum(truth <= est2[,2])==5){
    count[3] <- count[3] +1
  }
  est3 = multinomialCI(catet3_table[i,], 0.05)
  if(sum(est3[,1] <= truth)==5 & sum(truth <= est3[,2])==5){
    count[4] <- count[4] +1
  }
}

### summary table
catepeer_summary <- matrix(NA, nrow = 4, ncol = 3)
rownames(catepeer_summary) <- c("t=0", "t=1", "t=2", "t=3")
colnames(catepeer_summary) <-  c("95% CI coverage",
                                 "% of p-values(z)<=0.05", 
                                "% of p-values(permuation) <=0.05")
catepeer_summary[,1] <- count / np
catepeer_summary[,2] <- Zpower*100
catepeer_summary[,3] <- Ppower*100

catepeer_summary <- as.data.frame(catepeer_summary)
print(xtable(catepeer_summary , digits = 2, row.names = TRUE))


pdf("../figures/peercate_hist.pdf", width = 15, height = 8)
par(mfrow = c(1,1),   mar = c(7,10,5,3),  cex.lab = 3, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(density(phi0), 
     xlab = expression(Phi),
     main = expression(paste("Distribution of ", Phi, " under direct transmission")), 
     col = rgb(1,0,0,0.5), ylab = "Density",
     xlim = c(min(phi0), max(phi3)), 
      mgp = c(5,1,0), lwd = 4, xaxt = "n")
polygon(density(phi0), col=rgb(1,0,0,0.5))
axis(side = 1, at = seq(-20, 20, 1), 
     labels = seq(-20, 20, 1), 
     tck = 0.05)
lines(density(phi1), col= rgb(0,0,1,0.5), lwd = 4)
polygon(density(phi1), col=rgb(0,0,1,0.5))

lines(density(phi2), col = rgb(0,0,0,0.5), lwd = 4)
polygon(density(phi2), col=rgb(0,0,0,0.5))

lines(density(phi3), col = rgb(1,0.5,0,0.5), lwd = 4)
polygon(density(phi3), col=rgb(1,0.5,0,0.5))
abline(v = 1.645, lwd = 2, col = "red")
legend("topright", c("t = 0", "t = 1", "t = 2", "t = 3"),
       col = c(rgb(1,0,0,0.5), rgb(0, 0, 1,0.5), rgb(0,0,0,0.5), rgb(1,0.5,0,0.5)),
       seg.len = 2, 
       lty = c(1,1), lwd = 3, bty = 'n', cex = 3, xpd = NA)
dev.off()



