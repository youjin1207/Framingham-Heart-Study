library(igraph)
library(MASS)
library(parallel)
library(doParallel)
library(netdep)
## download network data
#network_c2 = read.table("phs000153_c2.txt", sep = "\t", header = TRUE)
#pheno.data = read.csv("pheno_c2_ex1_5.csv",sep = ",", header = TRUE)
pheno.info = pheno.data
focus.data = cbind(pheno.info$shareid, pheno.info$E487)
focus.data = as.data.frame(na.omit(focus.data))
names(focus.data) = c("shareid", "pheno")
focus.data = as.data.frame(na.omit(focus.data))
focus.data$pheno = ifelse(focus.data$pheno > 1, 1, focus.data$pheno)
network_c2_1.5 = network_c2 
Adj = matrix(0, nrow = nrow(focus.data), ncol = nrow(focus.data))
ids = focus.data$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network_c2_1.5[which(network_c2_1.5$shareid == ids[i]), 34])
  friends2 = which(ids %in% network_c2_1.5[which(network_c2_1.5$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i ] = 1
}
sample.size = nrow(Adj) # the number of nodes

G = graph.adjacency(Adj)
CC = components(G, mode = c("weak", "strong"))$membership

test0_cor = list()
kappa = c(1,2,3)
# k=1,2,3
for(r in 1:500){
  set.seed(r)
  outcome = focus.data$pheno
  corr.mat = Adj*0.15*kappa[k]
  for(i in 1:(nrow(Adj)-1)){
    for(j in (i+1):(nrow(Adj))){
      if(CC[i] == CC[j] & Adj[i,j] == 0){ # if connected, but not adjacent
        corr.mat[i,j] = corr.mat[j,i] =  0.10
      }else if(Adj[i,j] == 0){
        corr.mat[i,j] = corr.mat[j,i] =  0
      }
    }
  }
  diag(corr.mat) = 1
  mean.vector = rep(0, nrow(Adj))
  # generate network-dependent X and Y
  covariate = mvrnorm(1, mean.vector, Sigma = corr.mat)
  outcome = mvrnorm(1, mean.vector, Sigma = corr.mat)
  
  all.result = lm(outcome ~ covariate)
  all.p = summary(all.result)$coefficients[2,4]
  lower = summary(all.result)$coefficients[2,1] - 1.96*summary(all.result)$coefficients[2,2]
  upper = summary(all.result)$coefficients[2,1] + 1.96*summary(all.result)$coefficients[2,2]
  est = summary(all.result)$coefficients[2,1]
  pval = as.numeric(all.p)  
  se = summary(all.result)$coefficients[2,2]
  moran.cov = make.permute.moran(Adj, covariate, 500)
  moran.resi = make.permute.moran(Adj, all.result$residuals, 500)
  moran.outcome = make.permute.moran(Adj, outcome, 500)

  test0_cor[[r]] = list(moran.cov = moran.cov, moran.outcome = moran.outcome, glm.result = c(lower, upper, pval, moran.resi, est, se))
} 


### read the results
test0.ci = test1.ci = test2.ci = test3.ci = test4.ci = matrix(0, 500, 2)
test0.pval = test1.pval = test2.pval = test3.pval = test4.pval = c()
test0.est = test1.est = test2.est = test3.est = test4.est = c()
test0.moran.cov = test1.moran.cov = test2.moran.cov = test3.moran.cov = test4.moran.cov  = c()
test0.moran.resi = test1.moran.resi = test2.moran.resi = test3.moran.resi = test4.moran.resi =  c()
test0.sd = test1.sd = test2.sd = test3.sd = test4.sd = c()
test0.se = test1.se = test2.se = test3.se = test4.se = c()

load("Data/test_cor0.RData")
for(i in 1:500){
  test0.ci[i,] = test_cor0[[i]]$glm.result[1:2]
  test0.pval[i] = test_cor0[[i]]$glm.result[3]
  test0.moran.cov[i] = test_cor0[[i]]$moran.cov[3]
  test0.moran.outcome[i] = test_cor0[[i]]$moran.outcome[3]
  test0.moran.resi[i] = test_cor0[[i]]$glm.result[6]
  test0.est[i] = test_cor0[[i]]$glm.result[7]
  test0.sd[i] = test_cor0[[i]]$glm.result[8]
}
in.index = which(test0.ci[,1]*test0.ci[,2] < 0)
test0.out = test0.ci[-in.index,]
test0.in = test0.ci[in.index,]
test0.cci = rbind(test0.out, test0.in)
test0.se = sd(test0.est)

load("Data/test_cor1.RData")
for(i in 1:500){
  test1.ci[i,] = test_cor1[[i]]$glm.result[1:2]
  test1.pval[i] = test_cor1[[i]]$glm.result[3]
  test1.moran.cov[i] = test_cor1[[i]]$moran.cov[3]
  test1.moran.resi[i] = test_cor1[[i]]$glm.result[6]
  test1.est[i] = test_cor1[[i]]$glm.result[7]
  test1.sd[i] = test_cor1[[i]]$glm.result[8]
}
in.index = which(test1.ci[,1]*test1.ci[,2] < 0)
test1.out = test1.ci[-in.index,]
test1.in = test1.ci[in.index,]
test1.cci = rbind(test1.out, test1.in)
test1.se = sd(test1.est)

mean((test1.est-1.96*test1.se)*(test1.est+1.96*test1.se) < 0)

load("Data/test_cor2.RData")
for(i in 1:500){
  test2.ci[i,] = test_cor2[[i]]$glm.result[1:2]
  test2.pval[i] = test_cor2[[i]]$glm.result[3]
  test2.moran.cov[i] = test_cor2[[i]]$moran.cov[3]
  test2.moran.outcome[i] = test_cor2[[i]]$moran.outcome[3]
  test2.moran.resi[i] = test_cor2[[i]]$glm.result[6]
  test2.est[i] = test_cor2[[i]]$glm.result[7]
  test2.sd[i] = test_cor2[[i]]$glm.result[8]
}
in.index = which(test2.ci[,1]*test2.ci[,2] < 0)
test2.out = test2.ci[-in.index,]
test2.in = test2.ci[in.index,]
test2.cci = rbind(test2.out, test2.in)
test2.se = sd(test2.est)

load("Data/test_cor3.RData")
for(i in 1:500){
  test3.ci[i,] = test_cor3[[i]]$glm.result[1:2]
  test3.pval[i] = test_cor3[[i]]$glm.result[3]
  test3.moran.cov[i] = test_cor3[[i]]$moran.cov[3]
  test3.moran.outcome[i] = test_cor3[[i]]$moran.outcome[3]
  test3.moran.resi[i] = test_cor3[[i]]$glm.result[6]
  test3.est[i] = test_cor3[[i]]$glm.result[7]
  test3.sd[i] = test_cor3[[i]]$glm.result[8]
}
in.index = which(test3.ci[,1]*test3.ci[,2] < 0)
test3.out = test3.ci[-in.index,]
test3.in = test3.ci[in.index,]
test3.cci = rbind(test3.out, test3.in)
test3.se = sd(test3.est)

### table
moran.power.cov = c(mean(test0.moran.cov <= 0.05), mean(test1.moran.cov <= 0.05), 
                    mean(test2.moran.cov <= 0.05), mean(test3.moran.cov <= 0.05)) 
moran.power.resi = c(mean(test0.moran.resi <= 0.05), mean(test1.moran.resi <= 0.05), 
                     mean(test2.moran.resi <= 0.05), mean(test3.moran.resi <= 0.05)) 
moran.power.outcome = c(mean(test0.moran.outcome <= 0.05), mean(test1.moran.outcome <= 0.05), 
                        mean(test2.moran.outcome <= 0.05), mean(test3.moran.outcome <= 0.05)) 
sd.mean = c(mean(test0.sd), mean(test1.sd), mean(test2.sd), mean(test3.sd))
sem = c(test0.se, test1.se, test2.se, test3.se)
beta.coverage = c(nrow(test0.in), nrow(test1.in),  nrow(test2.in), nrow(test3.in))/500
beta.coverage = beta.coverage*100
coverages = beta.coverage
moran.est = c(mean(test0.moran.cov), mean(test1.moran.cov), mean(test2.moran.cov), mean(test3.moran.cov))
moran.est = formatC(moran.est, 3, format = "f")

tab = rbind(c(mean(test0.est), mean(test1.est), mean(test2.est), mean(test3.est)),
            c(mean(abs(test0.est)), mean(abs(test1.est)), mean(abs(test2.est)), mean(abs(test3.est))),
            sd.mean, sem, beta.coverage, moran.power.outcome, moran.power.cov, moran.power.resi)
rownames(tab) = c("bias", "|bias|", "average SE", "SD(est)", "Coverage rate", "power of Y",  "power of X", "power of residuals")
print(xtable(tab[,c(4,1,2,3)], digits = 4))

### make a figure
pdf("Figure/confounding.pdf", width = 12, height = 8)
par(mfrow = c(1,4), oma = c(3, 4, 4, 2), cex.lab = 1.8, 
    cex.main = 1.8, cex.axis = 2, tcl = 0.5,
    mai = c(1.0, 0.5, 0.5, 0.5))
plot(x = c(test3.cci[1,1] , test3.cci[1,2]), y = c(1.0, 1.0)   ,xlim = c(-max(abs(test2.cci)) , max(abs(test2.cci)) ),
     ylim = c(0,1.0), lwd = 0.1, type = "l", 
     xlab = expression(paste("Coverage of ", beta, "= 0 : ", 94.8, "%")),
     ylab="", col = "deepskyblue4", main = expression(paste("Permuted A", "(", kappa, "=3)")),
     mgp = c(6,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-max(abs(test2.cci))/2, 0, max(abs(test2.cci))/2),  labels = c(round(-max(abs(test2.cci))/2, 1), 0, round(max(abs(test2.cci))/2, 1)  ), tck = 0.05)
for(i in 2:nrow(test3.cci)){
  lines(x = c(test3.cci[i,1], test3.cci[i,2])
        , y = c(1.0 - (i-1)/nrow(test3.cci), 1.0 - (i-1)/nrow(test3.cci)),
        lwd = 0.1 , type = "l", col = "deepskyblue4")
}
abline(h = coverages[4]/100, lty =2, col = "red", lwd = 2)
abline(v = 0, lty = 1, col = "black", lwd = 2)

plot(x = c(test0.cci[1,1] , test0.cci[1,2]), y = c(1.0, 1.0), xlim = c(-max(abs(test2.cci)), max(abs(test2.cci))),
     ylim = c(0,1.0), lwd = 0.1, type = "l",  
     xlab = expression(paste("Coverage of ", beta, "= 0 : ", 88.6, "%")),
     ylab = "", col = "yellow2", main = expression(paste(kappa, "=1")),
     mgp = c(6,1,0), xpd = FALSE,  xaxt = 'n')
axis(1, at=c(-max(abs(test2.cci))/2, 0, max(abs(test2.cci))/2),  labels = c(round(-max(abs(test2.cci))/2, 1), 0, round(max(abs(test2.cci))/2, 1)  ), tck = 0.05)
for(i in 2:nrow(test0.cci)){
  lines(x = c(test0.cci[i,1], test0.cci[i,2])
        , y = c(1.0 - (i-1)/nrow(test0.cci), 1.0 - (i-1)/nrow(test0.cci)),
        lwd = 0.1 , type = "l", col = "yellow2")
}
abline(h = coverages[1]/100, lty =2, col = "red", lwd = 2)
abline(v = 0, lty = 1, col = "black", lwd = 2)

plot(x = c(test1.cci[1,1] , test1.cci[1,2]), y = c(1.0, 1.0)   ,xlim = c(-max(abs(test2.cci)), max(abs(test2.cci))),
     ylim = c(0,1.0), lwd = 0.1, type = "l",
     xlab = expression(paste("Coverage of ", beta, "= 0 : ", 79.2, "%")),
     ylab = "", col = "purple", main = expression(paste(kappa, "=2")),
     mgp = c(6,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-max(abs(test2.cci))/2, 0, max(abs(test2.cci))/2),  labels = c(round(-max(abs(test2.cci))/2, 1), 0, round(max(abs(test2.cci))/2, 1)  ), tck = 0.05)
for(i in 2:nrow(test1.cci)){
  lines(x = c(test1.cci[i,1], test1.cci[i,2])
        , y = c(1.0 - (i-1)/nrow(test1.cci), 1.0 - (i-1)/nrow(test1.cci)),
        lwd = 0.1 , type = "l", col = "purple")
}
abline(h = coverages[2]/100, lty =2, col = "red", lwd = 2)
abline(v = 0, lty = 1, col = "black", lwd = 2)

plot(x = c(test2.cci[1,1] , test2.cci[1,2]), y = c(1.0, 1.0)   ,xlim = c(-max(abs(test2.cci)) , max(abs(test2.cci))),
     ylim = c(0,1.0), lwd = 0.1, type = "l", 
     xlab = expression(paste("Coverage of ", beta, "= 0 : ", "71.0", "%")),
     ylab="", col = "springgreen4", main = expression(paste(kappa, "=3")),
     mgp = c(6,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-max(abs(test2.cci))/2, 0, max(abs(test2.cci))/2),  labels = c(round(-max(abs(test2.cci))/2, 1), 0, round(max(abs(test2.cci))/2, 1)  ), tck = 0.05)
for(i in 2:nrow(test2.cci)){
  lines(x = c(test2.cci[i,1], test2.cci[i,2])
        , y = c(1.0 - (i-1)/nrow(test2.cci), 1.0 - (i-1)/nrow(test2.cci)),
        lwd = 0.1 , type = "l", col = "springgreen4")
}
abline(h = coverages[3]/100, lty =2, col = "red", lwd = 2)
abline(v = 0, lty = 1, col = "black", lwd = 2)


mtext(expression(paste("95% confidence intervals for ", beta, " assuming independence", sep="")), side = 3, line = 0, outer = TRUE, cex= 2.5, xpd = TRUE)
mtext("Proportion of Simulations", side = 2, line = 1, outer = TRUE, cex= 2.5, adj = 0.7)
dev.off()