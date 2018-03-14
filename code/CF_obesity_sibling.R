library(igraph)
library(ColorPalette)
library(RColorBrewer)
GnPalette = colorRampPalette(brewer.pal(9,"Greens"))(20)

source("MoranI.R")
source("Phi.R")

# The time interval includes Original Cohort exams 12, 16, 19, 21, 23, 24, 26, 29 and Offspring Cohort exams 1 - 8. 
# Original Cohorts Exams : 12 (1971), 16, 19, 21, 23, 24, 26 (1999-2001)
# Offspring Cohorts Exams : 1 (1971-1975), 7 (1998 - 2001)

#network_c1 = read.table("data/phs000153.v9.pht000836.v8.p8.c1.vr_sntwk_2008_m_0641s.DS-SN-IRB-MDS.txt", sep = "\t", header = TRUE)
#network_c2 = read.table("data/phs000153_c2.txt", sep = "\t", header = TRUE)
network = rbind(network_c1, network_c2)

# collect pairs of brother or sister
brosis_c1 = network_c1[network_c1$ALTERTYPE == "BROTHER" | network_c1$ALTERTYPE == "SISTER", ]
brosis_c2 = network_c2[network_c2$ALTERTYPE == "BROTHER" | network_c2$ALTERTYPE == "SISTER",]
brosis = rbind(brosis_c1, brosis_c2)
sibling_c1 = cbind(brosis_c1$shareid, brosis_c1$sharealterid)
sibling_c2 = cbind(brosis_c2$shareid, brosis_c2$sharealterid)
sibling = rbind(sibling_c1, sibling_c2)
colnames(sibling) = c("shareid", "sharealterid")
sibling = as.data.frame(sibling)

## sibling is time-independent relationship so just relplicate seven times
tie.data = rbind(sibling, sibling,
                 sibling, sibling,
                 sibling, sibling,
                 sibling)
wave.ind = c(rep(1, nrow(sibling)), rep(2, nrow(sibling)),
             rep(3, nrow(sibling)), rep(4, nrow(sibling)),
             rep(5, nrow(sibling)), rep(6, nrow(sibling)),
             rep(7, nrow(sibling)))
tie.data = cbind(tie.data, wave.ind)
colnames(tie.data) = c("shareid", "sharealterid", "wave")
tie.data = as.data.frame(tie.data)

################# BMI status #########################
## c1 data
#wave1_c1 = read.table("data/phs000007.v29.pht000030.v7.p10.c1.ex1_1s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#wave2_c1 = read.table("data/phs000007.v29.pht000031.v7.p10.c1.ex1_2s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#wave3_c1 = read.table("data/phs000007.v29.pht000032.v6.p10.c1.ex1_3s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#wave4_c1 = read.table("data/phs000007.v29.pht000033.v8.p10.c1.ex1_4s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#wave5_c1 = read.table("data/phs000007.v29.pht000034.v7.p10.c1.ex1_5s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#wave6_c1 = read.table("data/phs000007.v29.pht000035.v8.p10.c1.ex1_6s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#wave7_c1 = read.table("data/phs000007.v29.pht000036.v8.p10.c1.ex1_7s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)


BMI1_c1 = (wave1_c1$A50)*703 / (wave1_c1$A51)^2
BMI2_c1 = (wave2_c1$B13) / (wave2_c1$B14 / 100)^2
BMI3_c1 = (wave3_c1$C416)*703 / (wave3_c1$C417)^2
BMI4_c1 = (wave4_c1$D401)*703 / (wave4_c1$D402)^2
BMI5_c1 = (wave5_c1$E024)*703 / (wave5_c1$E025)^2
BMI6_c1 = (wave6_c1$F007)*703 / (wave6_c1$F008)^2
BMI7_c1 = (wave7_c1$G440)*703 / (wave7_c1$G441)^2

BMI_c1 = rbind( cbind(wave1_c1$shareid, BMI1_c1, rep(1, nrow(wave1_c1))),
                cbind(wave2_c1$shareid, BMI2_c1, rep(2, nrow(wave2_c1))),
                cbind(wave3_c1$shareid, BMI3_c1, rep(3, nrow(wave3_c1))),
                cbind(wave4_c1$shareid, BMI4_c1, rep(4, nrow(wave4_c1))),
                cbind(wave5_c1$shareid, BMI5_c1, rep(5, nrow(wave5_c1))),
                cbind(wave6_c1$shareid, BMI6_c1, rep(6, nrow(wave6_c1))),
                cbind(wave7_c1$shareid, BMI7_c1, rep(7, nrow(wave7_c1))))
obesity.ind = ifelse(BMI_c1[,2] > 30, 1, 0)                
obesity_c1 = cbind(BMI_c1, obesity.ind) 
colnames(obesity_c1) = c("shareid", "BMI", "wave", "obesity")
obesity_c1 = as.data.frame(obesity_c1)
## c2 data
#wave1_c2 = read.table("data/phs000007.v29.pht000030.v7.p10.c2.ex1_1s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
#wave2_c2 = read.table("data/phs000007.v29.pht000031.v7.p10.c2.ex1_2s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
#wave3_c2 = read.table("data/phs000007.v29.pht000032.v6.p10.c2.ex1_3s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
#wave4_c2 = read.table("data/phs000007.v29.pht000033.v8.p10.c2.ex1_4s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
#wave5_c2 = read.table("data/phs000007.v29.pht000034.v7.p10.c2.ex1_5s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
#wave6_c2 = read.table("data/phs000007.v29.pht000035.v8.p10.c2.ex1_6s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
#wave7_c2 = read.table("data/phs000007.v29.pht000036.v8.p10.c2.ex1_7s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)

BMI1_c2 = (wave1_c2$A50)*703 / (wave1_c2$A51)^2
BMI2_c2 = (wave2_c2$B13) / (wave2_c2$B14 / 100)^2
BMI3_c2 = (wave3_c2$C416)*703 / (wave3_c2$C417)^2
BMI4_c2 = (wave4_c2$D401)*703 / (wave4_c2$D402)^2
BMI5_c2 = (wave5_c2$E024)*703 / (wave5_c2$E025)^2
BMI6_c2 = (wave6_c2$F007)*703 / (wave6_c2$F008)^2
BMI7_c2 = (wave7_c2$G440)*703 / (wave7_c2$G441)^2

BMI_c2 = rbind( cbind(wave1_c2$shareid, BMI1_c2, rep(1, nrow(wave1_c2))),
                cbind(wave2_c2$shareid, BMI2_c2, rep(2, nrow(wave2_c2))),
                cbind(wave3_c2$shareid, BMI3_c2, rep(3, nrow(wave3_c2))),
                cbind(wave4_c2$shareid, BMI4_c2, rep(4, nrow(wave4_c2))),
                cbind(wave5_c2$shareid, BMI5_c2, rep(5, nrow(wave5_c2))),
                cbind(wave6_c2$shareid, BMI6_c2, rep(6, nrow(wave6_c2))),
                cbind(wave7_c2$shareid, BMI7_c2, rep(7, nrow(wave7_c2))))
obesity.ind = ifelse(BMI_c2[,2] > 30, 1, 0)                
obesity_c2 = cbind(BMI_c2, obesity.ind) 
colnames(obesity_c2) = c("shareid", "BMI", "wave", "obesity")
obesity_c2 = as.data.frame(obesity_c2)


obesity.dat = rbind(obesity_c1, obesity_c2)
############# EDUCATION, AGE, SEX ###############
## year of education
edu_c1 = cbind(wave2_c1$shareid, wave2_c1$B43)
edu_c2 = cbind(wave2_c2$shareid, wave2_c2$B43)
edu.dat = rbind(edu_c1, edu_c2)
colnames(edu.dat) = c("shareid", "eduyear"); edu.dat = as.data.frame(edu.dat)
edu = rep(NA, nrow(tie.data))
for(i in 1:nrow(tie.data)){
  ind = which(edu.dat$shareid %in% tie.data$shareid[i])
  if(length(ind)!=0){
    edu[i] = edu.dat$eduyear[ind]
  }
}
## age and sex
#age.whole_c1 = read.table("data/age/phs000007.v29.pht003099.v4.p10.c1.vr_dates_2014_a_0912s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
age.off_c1 = age.whole_c1[age.whole_c1$idtype == 1,]
agesex_c1 = as.data.frame(cbind(age.off_c1$shareid, age.off_c1$sex, age.off_c1$age1,
                                age.off_c1$age2, age.off_c1$age3, age.off_c1$age4, 
                                age.off_c1$age5, age.off_c1$age6, age.off_c1$age7))
colnames(agesex_c1) = c("shareid", "sex", "age1", "age2", "age3", "age4", "age5", "age6", "age7")


## age and sex
#age.whole_c2 = read.table("data/age/phs000007.v29.pht003099.v4.p10.c2.vr_dates_2014_a_0912s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
age.off_c2 = age.whole_c2[age.whole_c2$idtype == 1,]
agesex_c2 = as.data.frame(cbind(age.off_c2$shareid, age.off_c2$sex, age.off_c2$age1,
                                age.off_c2$age2, age.off_c2$age3, age.off_c2$age4, 
                                age.off_c2$age5, age.off_c2$age6, age.off_c2$age7))
colnames(agesex_c2) = c("shareid", "sex", "age1", "age2", "age3", "age4", "age5", "age6", "age7")

## match education,age and sex
## age data at exam 1 are available for everyone so if ages are missing at other exams,
## we 'extrapolate' ages based on what we have at exam 1 for each subjects. 
agesex.dat = rbind(agesex_c1, agesex_c2)

for(i in 1:nrow(agesex.dat)){
  if(is.na(agesex.dat$age2[i])){
    agesex.dat$age2[i] = agesex.dat$age1[i] + 8
  }
  if(is.na(agesex.dat$age3[i])){
    agesex.dat$age3[i] = agesex.dat$age2[i] + 4
  }
  if(is.na(agesex.dat$age4[i])){
    agesex.dat$age4[i] = agesex.dat$age3[i] + 4
  }
  if(is.na(agesex.dat$age5[i])){
    agesex.dat$age5[i] = agesex.dat$age4[i] + 4
  }
  if(is.na(agesex.dat$age6[i])){
    agesex.dat$age6[i] = agesex.dat$age5[i] + 4
  }
  if(is.na(agesex.dat$age7[i])){
    agesex.dat$age7[i] = agesex.dat$age6[i] + 3
  }
}

age.current = rep(NA, nrow(tie.data)); sex = rep(NA, nrow(tie.data));
for(i in 1:nrow(tie.data)){
  ind = which(agesex.dat$shareid %in% tie.data$shareid[i])
  if(length(ind)!=0){
    age.current[i] = agesex.dat[ind,(tie.data$wave[i]+2)]
    sex[i] = agesex.dat$sex[ind]
  }
}

## make Y_{i,t}, Y_{j,t-1}, and Y_{j, t-2} at each wave t.
ego.obesity.current = rep(NA, nrow(tie.data)); ego.obesity.previous = rep(NA, nrow(tie.data));
alter.obesity.previous = rep(NA, nrow(tie.data)); alter.obesity.previous2 = rep(NA, nrow(tie.data))

for(i in 1:nrow(tie.data)){
  
  ego.current = which(obesity.dat$shareid %in% tie.data$shareid[i] & obesity.dat$wave == tie.data$wave[i])
  if(length(ego.current) == 1) ego.obesity.current[i] = obesity.dat$obesity[ego.current]
  
  ego.previous = which(obesity.dat$shareid %in% tie.data$shareid[i] & obesity.dat$wave == (tie.data$wave[i]-1))
  if(length(ego.previous) == 1) ego.obesity.previous[i] = obesity.dat$obesity[ego.previous]
  
  alter.previous = which(obesity.dat$shareid %in% tie.data$sharealterid[i] & obesity.dat$wave == (tie.data$wave[i]-1))
  if(length(alter.previous) == 1) alter.obesity.previous[i] = obesity.dat$obesity[alter.previous]
  
  alter.previous2 = which(obesity.dat$shareid %in% tie.data$sharealterid[i] & obesity.dat$wave == (tie.data$wave[i]-2))
  if(length(alter.previous2) == 1) alter.obesity.previous2[i] = obesity.dat$obesity[alter.previous2]
  
}


## merge every data for each pair of siblings.
all.data = cbind(tie.data, sex, age.current,edu,
                 ego.obesity.current, ego.obesity.previous,
                 alter.obesity.previous, alter.obesity.previous2)
wave1.data = all.data[all.data$wave == 1,]
wave2.data = all.data[all.data$wave == 2,]
wave3.data = all.data[all.data$wave == 3,]
wave4.data = all.data[all.data$wave == 4,]
wave5.data = all.data[all.data$wave == 5,]
wave6.data = all.data[all.data$wave == 6,]
wave7.data = all.data[all.data$wave == 7,]

#### wave 3 ####
ind = c()
for(i in 1:length(unique(wave3.data$shareid))){
  ind[i] = which(wave3.data$shareid == unique(wave3.data$shareid)[i])[1]
}
## choose one ego within one wave.
wave3.data.onepair = wave3.data[ind, ]
wave3.data.onepair = wave3.data.onepair[!is.na(wave3.data.onepair$ego.obesity.current), ]

ids = unique(wave3.data.onepair$shareid)
Adj = matrix(0, nrow = length(ids), ncol = length(ids))
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% wave3.data$sharealterid[which(wave3.data$shareid == ids[i])])
  friends2 = which(ids %in% wave3.data$shareid[which(wave3.data$sharealterid == ids[i])])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}
sibling.Adj.wave3 = Adj
sibling.wave3 = make.permute.moran(Adj, wave3.data.onepair$ego.obesity.current, 500)
n.wave3 = length(wave3.data.onepair$ego.obesity.current)

## fit logistic regression
logit.fit = glm(ego.obesity.current ~ age.current + sex + edu + 
                  ego.obesity.previous + alter.obesity.previous + alter.obesity.previous2,
                data = wave3.data.onepair, family = binomial)
no.na = names(logit.fit$residuals) # extract subjects who do not have any missing information.
wave3.data.onepair.logit = wave3.data.onepair[no.na,]

ids = unique(wave3.data.onepair.logit$shareid)
Adj = matrix(0, nrow = length(ids), ncol = length(ids))
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% wave3.data$sharealterid[which(wave3.data$shareid == ids[i])])
  friends2 = which(ids %in% wave3.data$shareid[which(wave3.data$sharealterid == ids[i])])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}
sibling.Adj.logit.wave3 = Adj
sibling.logit.wave3 = make.permute.moran(Adj, logit.fit$residuals, 500)
logit.n.wave3 = nrow(wave3.data.onepair.logit)


sibling.wave3.ego.obesity.previous = make.permute.moran(Adj, wave3.data.onepair$ego.obesity.previous, 500)
sibling.wave3.alter.obesity.previous = make.permute.moran(Adj, wave3.data.onepair$alter.obesity.previous, 500)
sibling.wave3.alter.obesity.previous2 = make.permute.moran(Adj, wave3.data.onepair$alter.obesity.previous2, 500)
sibling.wave3.alter.obesity.previous = make.permute.moran(Adj, wave3.data.onepair$alter.obesity.previous, 500)
sibling.wave3.age = make.permute.moran(Adj, wave3.data.onepair$age.current, 500)
sibling.wave3.sex = make.permute.moran(Adj, wave3.data.onepair$sex, 500)
sibling.wave3.edu = make.permute.moran(Adj, wave3.data.onepair$edu, 500)

### visualization ###
## applied to binary obesity status 
igraph.options(vertex.size = 3, edge.arrow.size = 0.1,
               vertex.label = NULL)
sibling.G.wave3 = graph.adjacency(sibling.Adj.wave3, "undirected")
V(sibling.G.wave3)$Y = wave3.data.onepair$ego.obesity.current

set.seed(123)
g = induced_subgraph(sibling.G.wave3, sample(c(1:length(V(sibling.G.wave3))), 1000))
E(g)$color = "grey"
ind = which(colSums(as.matrix(get.adjacency(g))) > 0)
subsibling.Adj.wave3 = as.matrix(get.adjacency(g))[ind, ind]
subsibling.G.wave3 = graph.adjacency(subsibling.Adj.wave3, "undirected")
V(subsibling.G.wave3)$Y = V(g)$Y[ind]
V(subsibling.G.wave3)$color = ifelse(V(subsibling.G.wave3)$Y == 1, "seagreen", "seashell")
E(subsibling.G.wave3)$color = "black"
pdf("../figures/subsibling_wave3_obesity.pdf")
igraph.options(vertex.size = 2.5, edge.arrow.size = 0.1,
               vertex.label = NULL)
par(mar=c(0,0,3,0), cex.main = 3)
set.seed(123)
plot(subsibling.G.wave3, layout = layout.fruchterman.reingold,
     vertex.label = "", main = "Sibling Relationships")
dev.off()

## applied to the residuals
igraph.options(vertex.size = 5, edge.arrow.size = 0.1,
               vertex.label = NULL)
sibling.logit.G.wave3  = graph.adjacency(sibling.Adj.logit.wave, "undirected")
V(sibling.logit.G.wave3)$Y = wave3.resi
for(i in 1:20){
  lower = quantile(V(sibling.logit.G.wave3)$Y, 0.05*(i-1))
  upper = quantile(V(sibling.logit.G.wave3)$Y, 0.05*i)
  V(sibling.logit.G.wave3)[V(sibling.logit.G.wave3)$Y >= lower & V(sibling.logit.G.wave3)$Y < upper]$color = GnPalette[i]
}
V(sibling.logit.G.wave3)[which.max(V(sibling.logit.G.wave3)$Y)]$color = GnPalette[20]

set.seed(123)
g = induced_subgraph(sibling.logit.G.wave3, sample(c(1:length(V(sibling.logit.G.wave3))), 500))
E(g)$color = "grey"
ind = which(colSums(as.matrix(get.adjacency(g))) > 0)
subsibling.Adj.logit.wave = as.matrix(get.adjacency(g))[ind, ind]
subsibling.logit.G.wave =  graph.adjacency(subsibling.Adj.logit.wave, "undirected")
V(subsibling.logit.G.wave)$Y = wave3.resi[ind]
for(i in 1:20){
  lower = quantile(V(subsibling.logit.G.wave)$Y, 0.05*(i-1))
  upper = quantile(V(subsibling.logit.G.wave)$Y, 0.05*i)
  V(subsibling.logit.G.wave)[V(subsibling.logit.G.wave)$Y >= lower & V(subsibling.logit.G.wave)$Y < upper]$color = GnPalette[i]
}
V(subsibling.logit.G.wave)[which.max(V(subsibling.logit.G.wave)$Y)]$color = GnPalette[20]

pdf("../figures/subsibling_logit_wave3.pdf")
E(subsibling.logit.G.wave)$color = "black"
igraph.options(vertex.size = 3, edge.arrow.width = 2,
               vertex.label = NULL, edge.arrow.col = "black")
par(mar=c(0,0,3,0), cex.main = 3)
set.seed(123)
plot(subsibling.logit.G.wave , layout = layout.fruchterman.reingold,
     vertex.label = "", main = "Sibling Relationships")
dev.off()

####wave 4 ####
ind = c()
for(i in 1:length(unique(wave4.data$shareid))){
  ind[i] = which(wave4.data$shareid == unique(wave4.data$shareid)[i])[1]
}

wave4.data.onepair = wave4.data[ind, ]
wave4.data.onepair = wave4.data.onepair[!is.na(wave4.data.onepair$ego.obesity.current), ]

ids = unique(wave4.data.onepair$shareid)
Adj = matrix(0, nrow = length(ids), ncol = length(ids))
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% wave4.data$sharealterid[which(wave4.data$shareid == ids[i])])
  friends2 = which(ids %in% wave4.data$shareid[which(wave4.data$sharealterid == ids[i])])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

sibling.wave4 = make.permute.moran(Adj, wave4.data.onepair$ego.obesity.current, 500)
n.wave4 = length(wave4.data.onepair$ego.obesity.current)

## logistic regression
logit.fit = glm(ego.obesity.current ~ age.current + sex + edu + 
                  ego.obesity.previous + alter.obesity.previous + alter.obesity.previous2,
                data = wave4.data.onepair, family = binomial)
no.na = names(logit.fit$residuals)
wave4.data.onepair.logit = wave4.data.onepair[no.na,]

ids = unique(wave4.data.onepair.logit$shareid)
Adj = matrix(0, nrow = length(ids), ncol = length(ids))
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% wave4.data$sharealterid[which(wave4.data$shareid == ids[i])])
  friends2 = which(ids %in% wave4.data$shareid[which(wave4.data$sharealterid == ids[i])])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

sibling.logit.wave4 = make.permute.moran(Adj, logit.fit$residuals, 500)
logit.n.wave4 = nrow(wave4.data.onepair.logit)

#### wave 5 ####
ind = c()
for(i in 1:length(unique(wave5.data$shareid))){
  ind[i] = which(wave5.data$shareid == unique(wave5.data$shareid)[i])[1]
}

wave5.data.onepair = wave5.data[ind, ]
wave5.data.onepair = wave5.data.onepair[!is.na(wave5.data.onepair$ego.obesity.current), ]

ids = unique(wave5.data.onepair$shareid)
Adj = matrix(0, nrow = length(ids), ncol = length(ids))
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% wave5.data$sharealterid[which(wave5.data$shareid == ids[i])])
  friends2 = which(ids %in% wave5.data$shareid[which(wave5.data$sharealterid == ids[i])])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

sibling.wave5 = make.permute.moran(Adj, wave5.data.onepair$ego.obesity.current, 500)
n.wave5 = length(wave5.data.onepair$ego.obesity.current)

## logistic regression
logit.fit = glm(ego.obesity.current ~ age.current + sex + edu + 
                  ego.obesity.previous + alter.obesity.previous + alter.obesity.previous2,
                data = wave5.data.onepair, family = binomial)
no.na = names(logit.fit$residuals)
wave5.data.onepair.logit = wave5.data.onepair[no.na,]

ids = unique(wave5.data.onepair.logit$shareid)
Adj = matrix(0, nrow = length(ids), ncol = length(ids))
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% wave5.data$sharealterid[which(wave5.data$shareid == ids[i])])
  friends2 = which(ids %in% wave5.data$shareid[which(wave5.data$sharealterid == ids[i])])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

sibling.logit.wave5 = make.permute.moran(Adj, logit.fit$residuals, 500)
logit.n.wave5 = nrow(wave5.data.onepair.logit)

#### wave 6 ####
ind = c()
for(i in 1:length(unique(wave6.data$shareid))){
  ind[i] = which(wave6.data$shareid == unique(wave6.data$shareid)[i])[1]
}

wave6.data.onepair = wave6.data[ind, ]
wave6.data.onepair = wave6.data.onepair[!is.na(wave6.data.onepair$ego.obesity.current), ]

ids = unique(wave6.data.onepair$shareid)
Adj = matrix(0, nrow = length(ids), ncol = length(ids))
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% wave6.data$sharealterid[which(wave6.data$shareid == ids[i])])
  friends2 = which(ids %in% wave6.data$shareid[which(wave6.data$sharealterid == ids[i])])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

sibling.wave6 = make.permute.moran(Adj, wave6.data.onepair$ego.obesity.current, 500)
n.wave6 = length(wave6.data.onepair$ego.obesity.current)

## logistic regression
logit.fit = glm(ego.obesity.current ~ age.current + sex + edu + 
                  ego.obesity.previous + alter.obesity.previous + alter.obesity.previous2,
                data = wave6.data.onepair, family = binomial)
no.na = names(logit.fit$residuals)
wave6.data.onepair.logit = wave6.data.onepair[no.na,]

ids = unique(wave6.data.onepair.logit$shareid)
Adj = matrix(0, nrow = length(ids), ncol = length(ids))
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% wave6.data$sharealterid[which(wave6.data$shareid == ids[i])])
  friends2 = which(ids %in% wave6.data$shareid[which(wave6.data$sharealterid == ids[i])])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

sibling.logit.wave6 = make.permute.moran(Adj, logit.fit$residuals, 500)
logit.n.wave6 = nrow(wave6.data.onepair.logit)

#### wave 7 ####
ind = c()
for(i in 1:length(unique(wave7.data$shareid))){
  ind[i] = which(wave7.data$shareid == unique(wave7.data$shareid)[i])[1]
}

wave7.data.onepair = wave7.data[ind, ]
wave7.data.onepair = wave7.data.onepair[!is.na(wave7.data.onepair$ego.obesity.current), ]

ids = unique(wave7.data.onepair$shareid)
Adj = matrix(0, nrow = length(ids), ncol = length(ids))
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% wave7.data$sharealterid[which(wave7.data$shareid == ids[i])])
  friends2 = which(ids %in% wave7.data$shareid[which(wave7.data$sharealterid == ids[i])])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

sibling.wave7 = make.permute.moran(Adj, wave7.data.onepair$ego.obesity.current, 500)
n.wave7 = length(wave7.data.onepair$ego.obesity.current)

## logistic regression
logit.fit = glm(ego.obesity.current ~ age.current + sex + edu + 
                  ego.obesity.previous + alter.obesity.previous + alter.obesity.previous2,
                data = wave7.data.onepair, family = binomial)
no.na = names(logit.fit$residuals)
wave7.data.onepair.logit = wave7.data.onepair[no.na,]

ids = unique(wave7.data.onepair.logit$shareid)
Adj = matrix(0, nrow = length(ids), ncol = length(ids))
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% wave7.data$sharealterid[which(wave7.data$shareid == ids[i])])
  friends2 = which(ids %in% wave7.data$shareid[which(wave7.data$sharealterid == ids[i])])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

sibling.logit.wave7 = make.permute.moran(Adj, logit.fit$residuals, 500)
logit.n.wave7 = nrow(wave7.data.onepair.logit)

######################################
sibling3 = cbind(sibling.wave3, n.wave3, sibling.logit.wave3, logit.n.wave3)
sibling4 = cbind(sibling.wave4, n.wave4, sibling.logit.wave4, logit.n.wave4)
sibling5 = cbind(sibling.wave5, n.wave5, sibling.logit.wave5, logit.n.wave5)
sibling6 = cbind(sibling.wave6, n.wave6, sibling.logit.wave6, logit.n.wave6)
sibling7 = cbind(sibling.wave7, n.wave7, sibling.logit.wave7, logit.n.wave7)


tab = rbind(sibling3, sibling4, sibling5, sibling6, sibling7)
colnames(tab) = c("W vs. Y", "n1", "W vs. Residuals", "n2")
print(xtable(tab))