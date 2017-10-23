library(igraph)
source("Phi.R")
source("MoranI.R")
############### Highest degrees (Offspring Cohort 3) ##############################
#network_c2 = read.table("data/phs000153_c2.txt", sep = "\t", header = TRUE)
psy.data = read.table("data/phs000007.v29.pht000100.v6.p10.c2.psych1_3s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
psy.data = psy.data[!is.na(psy.data$PY123),]
# approximate network ties existing during Offspring Cohort 3
network_c2_1.3 = network_c2[(network_c2$SPELLBEGIN <= 12*(87-71+1)) & (network_c2$SPELLEND > 12*(83-71)),  ]

Adj = matrix(0, nrow = nrow(psy.data), ncol = nrow(psy.data))
ids = psy.data$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network_c2_1.3[which(network_c2_1.3$shareid == ids[i]), 34])
  friends2 = which(ids %in% network_c2_1.3[which(network_c2_1.3$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i ] = 1
}

sum(Adj)/(2*nrow(Adj)) # network density

phi.degree = std.newphi(Adj, psy.data$PY123)
phi.degree.permute = make.permute.Phi(Adj, psy.data$PY123, 500)
moran.degree.permute = make.permute.moran(Adj, psy.data$PY123, 500)

## deleting category 6.none of the above
#psy.data = read.table("data/phs000007.v29.pht000100.v6.p10.c2.psych1_3s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
psy.data = psy.data[!is.na(psy.data$PY123),]
psy.data = psy.data[psy.data$PY123!=6,]
network_c2_1.3 = network_c2[(network_c2$SPELLBEGIN <= 12*(87-71+1)) & (network_c2$SPELLEND > 12*(83-71)),  ]

Adj = matrix(0, nrow = nrow(psy.data), ncol = nrow(psy.data))
ids = psy.data$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network_c2_1.3[which(network_c2_1.3$shareid == ids[i]), 34])
  friends2 = which(ids %in% network_c2_1.3[which(network_c2_1.3$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i ] = 1
}
sum(Adj)/(2*nrow(Adj))

phi.degree.complete = std.newphi(Adj, psy.data$PY123)
phi.degree.permute.complete = make.permute.Phi(Adj, psy.data$PY123, 500)
moran.degree.permute.complete = make.permute.moran(Adj, psy.data$PY123, 500)

############### Y : Kind of cold breakfast cereal (Offspring Cohort 5) ##############################
#food.data = read.table("data/phs000007.v29.pht000680.v5.p10.c2.ffreq1_5s.HMB-IRB-NPU-MDS.txt", sep ="\t", header = TRUE)
table(food.data$CER) # 
food.data = food.data[!is.na(food.data$CER),]
food.data = food.data[which(food.data$CER %in% c(0, 9, 14, 20, 31, 36, 39, 71, 76, 77, 79, 87, 93) ),]
table(food.data$CER)
# approximate network ties existing during Offspring Cohort Exam 5 : 1995 - 1998
network_c2_1.5 = network_c2[(network_c2$SPELLBEGIN <= 12*(95-71+1)) & (network_c2$SPELLEND > 12*(91-71)),  ]
Adj = matrix(0, nrow = nrow(food.data), ncol = nrow(food.data))
ids = food.data$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network_c2_1.5[which(network_c2_1.5$shareid == ids[i]), 34])
  friends2 = which(ids %in% network_c2_1.5[which(network_c2_1.5$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i ] = 1
}
sum(Adj)/(2*nrow(Adj)) # network denstiy

phi.cer = std.newphi(Adj, food.data$CER) 
phi.cer.permute = make.permute.Phi(Adj, food.data$CER, 500)
#moran.cer.permute = make.permute.moran(Adj, food.data$CER, 500)
