library(xtable)
source("MoranI.R")

## read data
#network_c2 = read.table("data/phs000153_c2.txt", sep = "\t", header = TRUE)
#psy.data = read.table("data/phs000007.v29.pht000100.v6.p10.c2.psych1_3s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
network2 = network_c2[(network_c2$SPELLBEGIN == 1),]
#pheno_c2_ex1_5 = read.table("data/phs000007.v29.pht000034.v7.p10.c2.ex1_5s.HMB-IRB-NPU-MDS.txt",sep = "\t", header = TRUE)
#pheno.data = read.table("data/phs000007.v29.pht000034.v7.p10.c2.ex1_5s.HMB-IRB-NPU-MDS.txt",sep = "\t", header = TRUE)
pheno.info = pheno.data


bp.ci = employ.ci = illness.ci = eye.ci = matrix(0, 500, 2)
bp.pval = employ.pval = illness.pval = eye.pval = c()
bp.moran = employ.moran = illness.moran = eye.moran = c()
bp.est = employ.est = illness.est = eye.est = c()
## 1. BLOOD PRESSURE: SYSTOLIC - 1ST MD READING 
focus.data = cbind(pheno.info$shareid, pheno.info$E485)
focus.data = as.data.frame(na.omit(focus.data))
names(focus.data) = c("shareid", "pheno")
focus.data = as.data.frame(na.omit(focus.data))
network_c2_1.5 = network_c2[(network_c2$SPELLBEGIN <= 12*(95-71+1)) & (network_c2$SPELLEND > 12*(91-71)),  ]
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
make.permute.moran(Adj, focus.data$pheno, 500)
sample.size[1] = nrow(Adj)
num.edges[1] = sum(Adj) / 2
covariate = c();result = c(); moran.resi = c()
outcome = focus.data$pheno

for(r in 1:500){
	set.seed(r)
  	outcome = focus.data$pheno
  	outcome = (outcome - mean(outcome))/sd(outcome) 

  	corr.mat = Adj*0.2
  	corr.mat = ifelse(corr.mat < 0.2, 0.1, corr.mat)
  	diag(corr.mat) = 0.5
  	mean.vector = rep(1, nrow(Adj))
  	mean.vector[which(rowSums(Adj) == 0)] = -1
 
  	covariate = mvrnorm(1, mean.vector, Sigma = corr.mat)
  	all.result = lm(outcome ~ covariate)
  	all.p = summary(all.result)$coefficients[2,4]
  	bp.ci[r,1] = summary(all.result)$coefficients[2,1] - 1.96*summary(all.result)$coefficients[2,2]
  	bp.ci[r,2] = summary(all.result)$coefficients[2,1] + 1.96*summary(all.result)$coefficients[2,2]
  	bp.pval[r] = as.numeric(all.p)  
  	bp.moran[r] = make.permute.moran(Adj, summary(all.result)$residuals, 500)[3]
  	bp.est[r] = summary(all.result)$coefficients[2,1] 
} 


## 2. Employed
focus.data = cbind(pheno.info$shareid, pheno.info$E067)
focus.data = as.data.frame(na.omit(focus.data))
names(focus.data) = c("shareid", "pheno")
focus.data = as.data.frame(na.omit(focus.data))
focus.data$pheno = ifelse(focus.data$pheno == 0 , 0, 1)
network_c2_1.5 = network_c2[(network_c2$SPELLBEGIN <= 12*(95-71+1)) & (network_c2$SPELLEND > 12*(91-71)),  ]
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
outcome = focus.data$pheno
sample.size[2] = nrow(Adj)
num.edges[2] = sum(Adj) / 2

for(r in 1:500){
  set.seed(r)
  outcome = focus.data$pheno
  corr.mat = Adj*0.2
  corr.mat = ifelse(corr.mat < 0.2, 0.1, corr.mat)
  diag(corr.mat) = 0.5
  mean.vector = rep(1, nrow(Adj))
  mean.vector[which(rowSums(Adj) == 0)] = -1
  
  covariate = mvrnorm(1, mean.vector, Sigma = corr.mat)

  all.result = glm(outcome ~ covariate, family = binomial)
  all.p = summary(all.result)$coefficients[2,4]
  employ.ci[r,1] = summary(all.result)$coefficients[2,1] - 1.96*summary(all.result)$coefficients[2,2]
  employ.ci[r,2] = summary(all.result)$coefficients[2,1] + 1.96*summary(all.result)$coefficients[2,2]
  employ.pval[r] = as.numeric(all.p)   
  employ.moran[r] = make.permute.moran(Adj, all.result$residuals, 500)[3]
  employ.est[r] = summary(all.result)$coefficients[2,1]
} 


## 3. Visited Doctor
focus.data = cbind(pheno.info$shareid, pheno.info$E215)
focus.data = as.data.frame(na.omit(focus.data))
names(focus.data) = c("shareid", "pheno")
focus.data$pheno = ifelse(focus.data$pheno != 0, 1, focus.data$pheno)
focus.data = as.data.frame(na.omit(focus.data))
network_c2_1.5 = network_c2[(network_c2$SPELLBEGIN <= 12*(95-71+1)) & (network_c2$SPELLEND > 12*(91-71)),  ]
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
outcome = focus.data$pheno
sample.size[3] = nrow(Adj)
num.edges[3] = sum(Adj) / 2

for(r in 1:500){
  set.seed(r)
  outcome = focus.data$pheno
  covariate = c()
  corr.mat = Adj*0.2
  corr.mat = ifelse(corr.mat < 0.2, 0.1, corr.mat)
  diag(corr.mat) = 0.5
  mean.vector = rep(1, nrow(Adj))
  mean.vector[which(rowSums(Adj) == 0)] = -1
  covariate = mvrnorm(1, mean.vector, Sigma = corr.mat)
 
  all.result = glm(outcome ~ covariate, family = binomial)
  all.p = summary(all.result)$coefficients[2,4]
  illness.ci[r,1] = summary(all.result)$coefficients[2,1] - 1.96*summary(all.result)$coefficients[2,2]
  illness.ci[r,2] = summary(all.result)$coefficients[2,1] + 1.96*summary(all.result)$coefficients[2,2]
  illness.pval[r] = as.numeric(all.p)   
  illness.moran[r] = make.permute.moran(Adj, all.result$residuals, 500)[3]
  illness.est[r] = summary(all.result)$coefficients[2,1]
} 


## 4. Corneal Arcus   
focus.data = cbind(pheno.info$shareid, pheno.info$E487)
focus.data = as.data.frame(na.omit(focus.data))
names(focus.data) = c("shareid", "pheno")
focus.data = as.data.frame(na.omit(focus.data))
focus.data$pheno = ifelse(focus.data$pheno > 1, 1, focus.data$pheno)
network_c2_1.5 = network_c2[(network_c2$SPELLBEGIN <= 12*(95-71+1)) & (network_c2$SPELLEND > 12*(91-71)),  ]
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
sample.size[4] = nrow(Adj)
num.edges[4] = sum(Adj) / 2
outcome = focus.data$pheno

for(r in 1:500){
  set.seed(r)
  outcome = focus.data$pheno
  corr.mat = Adj*0.2
  corr.mat = ifelse(corr.mat < 0.2, 0.1, corr.mat)
  diag(corr.mat) = 0.5
  mean.vector = rep(1, nrow(Adj))
  mean.vector[which(rowSums(Adj) == 0)] = -1
  covariate = mvrnorm(1, mean.vector, Sigma = corr.mat)

  all.result = glm(outcome ~ covariate, family = binomial)
  all.p = summary(all.result)$coefficients[2,4]
  eye.ci[r,1] = summary(all.result)$coefficients[2,1] - 1.96*summary(all.result)$coefficients[2,2]
  eye.ci[r,2] = summary(all.result)$coefficients[2,1] + 1.96*summary(all.result)$coefficients[2,2]
  eye.pval[r] = as.numeric(all.p)  
  eye.moran[r] = make.permute.moran(Adj, all.result$residuals, 500)[3]

  eye.est[r] = summary(all.result)$coefficients[2,1]
} 


### summary table (in Supplementary Material)
network.summary = matrix(0, 2, 4)
network.summary[1,] = as.integer(sample.size)
network.summary[2,] = as.integer(num.edges)
rownames(network.summary) = c("Sample size (n)", "The number of edges")
colnames(network.summary) = c("Systolic blood pressure", "Employed",
                              "Visited to doctor", "Corneal arcus")
xtable(network.summary)
###


## make figure
###### Confidence Interval Plots #############
## figure
moran.power = c(mean(bp.moran <= 0.05), mean(employ.moran <= 0.05), 
                mean(illness.moran <= 0.05), mean(eye.moran <= 0.05)) 
beta.coverage = c(nrow(bp.in), nrow(employ.in),  nrow(illness.in), nrow(eye.in))/500
moran.power = moran.power*100
beta.coverage = beta.coverage*100
coverages = beta.coverage
moran.est = c(mean(bp.moran), mean(employ.moran), mean(illness.moran), mean(eye.moran))
moran.est = formatC(moran.est, 3, format = "f")
beta.coverage = formatC(beta.coverage, 1, format = "f")
pdf("../figures/fiveCIs.pdf", width = 12, height = 8)
par(mfrow = c(1,4), oma = c(3, 4, 4, 2), cex.lab = 1.8, 
    cex.main = 1.8, cex.axis = 2, tcl = 0.5,
    mai = c(1.2, 0.5, 0.5, 0.5))
# blood pressure
plot(x = c(bp.cci[1,1] , bp.cci[1,2]), y = c(1.0, 1.0), xlim = c(-max(abs(bp.cci)), max(abs(bp.cci))),
     ylim = c(0,1.0), lwd = 0.1, type = "l",  
     xlab = bquote(atop(paste("Coverage of ", beta, "= 0 : ", 85.2, "%"),  
              atop(bold("Average p-value for tests of"), bold("network dependence: 0.037")))),
     ylab = "", col = "palegreen4", main = "Systolic blood pressure",
     mgp = c(8,1,0), xpd = FALSE,  xaxt = 'n')
axis(1, at=c(-max(abs(bp.cci))/2, 0, max(abs(bp.cci))/2),  labels = c(round(-max(abs(bp.cci))/2, 1), 0, round(max(abs(bp.cci))/2, 1)  ), tck = 0.05)
for(i in 2:nrow(bp.cci)){
  lines(x = c(bp.cci[i,1], bp.cci[i,2])
        , y = c(1.0 - (i-1)/nrow(bp.cci), 1.0 - (i-1)/nrow(bp.cci)),
        lwd = 0.1 , type = "l", col = "palegreen4")
}
abline(h = coverages[1]/100, lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# employ
plot(x = c(employ.cci[1,1] , employ.cci[1,2]), y = c(1.0, 1.0)   ,xlim = c(-max(abs(employ.cci)), max(abs(employ.cci))),
     ylim = c(0,1.0), lwd = 0.1, type = "l",
     xlab =  bquote(atop(paste("Coverage of ", beta, "= 0 : ", 69.2, "%"),  
                         atop(bold("Average p-value for tests of"), 
                              bold("network dependence: 0.002")))), 
     ylab = "", col = "royalblue", main = "Employed",
     mgp = c(8,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-max(abs(employ.cci))/2, 0, max(abs(employ.cci))/2),  labels = c(round(-max(abs(employ.cci))/2, 1), 0, round(max(abs(employ.cci))/2, 1)  ), tck = 0.05)
for(i in 2:nrow(employ.cci)){
  lines(x = c(employ.cci[i,1], employ.cci[i,2])
        , y = c(1.0 - (i-1)/nrow(employ.cci), 1.0 - (i-1)/nrow(employ.cci)),
        lwd = 0.1 , type = "l", col = "royalblue")
}
abline(h = coverages[2]/100, lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# illness consumption
plot(x = c(illness.cci[1,1] , illness.cci[1,2]), y = c(1.0, 1.0)   ,xlim = c(-max(abs(illness.cci)) , max(abs(illness.cci))),
     ylim = c(0,1.0), lwd = 0.1, type = "l", 
     xlab =  bquote(atop(paste("Coverage of ", beta, "= 0 : ", 65.2, "%"),  
                         atop(bold("Average p-value for tests of"), 
                              bold("network dependence: 0.697")))), 
     ylab="", col = "tan1", main = "Visited doctor",
     mgp = c(8,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-max(abs(illness.cci))/2, 0, max(abs(illness.cci))/2),  labels = c(round(-max(abs(illness.cci))/2, 1), 0, round(max(abs(illness.cci))/2, 1)  ), tck = 0.05)
for(i in 2:nrow(illness.cci)){
  lines(x = c(illness.cci[i,1], illness.cci[i,2])
        , y = c(1.0 - (i-1)/nrow(illness.cci), 1.0 - (i-1)/nrow(illness.cci)),
        lwd = 0.1 , type = "l", col = "tan1")
}
abline(h = coverages[3]/100, lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)

# eye
plot(x = c(eye.cci[1,1] , eye.cci[1,2]), y = c(1.0, 1.0)   ,xlim = c(-max(abs(eye.cci)) , max(abs(eye.cci)) ),
     ylim = c(0,1.0), lwd = 0.1, type = "l", 
     xlab = bquote(atop(paste("Coverage of ", beta, "= 0 : ", 78.6, "%"),  
                         atop(bold("Average p-value for tests of"), 
                              bold("network dependence: 0.018")))), 
     ylab="", col = "mediumpurple4", main = "Corneal arcus",
     mgp = c(8,1,0), xpd = FALSE, xaxt = 'n')
axis(1, at=c(-max(abs(eye.cci))/2, 0, max(abs(eye.cci))/2),  labels = c(round(-max(abs(eye.cci))/2, 1), 0, round(max(abs(eye.cci))/2, 1)  ), tck = 0.05)
for(i in 2:nrow(eye.cci)){
  lines(x = c(eye.cci[i,1], eye.cci[i,2])
        , y = c(1.0 - (i-1)/nrow(eye.cci), 1.0 - (i-1)/nrow(eye.cci)),
        lwd = 0.1 , type = "l", col = "mediumpurple4")
}
abline(h = coverages[4]/100, lty =2, col = "black", lwd = 2)
abline(v = 0, lty = 1, col = "red", lwd = 2)
mtext(expression(paste("95% confidence intervals for ", beta, " assuming independence", sep="")), side = 3, line = 0, outer = TRUE, cex= 2.5, xpd = TRUE)
mtext("Proportion of Simulations", side = 2, line = 1, outer = TRUE, cex= 2.5, adj = 0.7)
dev.off()
