library(xtable)
source("MoranI.R")

#pheno_c2_ex0_11 = read.table("data/phs000007.v29.pht000013.v3.p10.c2.ex0_11s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
#network_c2 = read.table("data/phs000153.v9.pht000836.v8.p8.c2.vr_sntwk_2008_m_0641s.DS-SN-IRB-MDS.txt", sep = "\t", header = TRUE)

# The time interval includes Original Cohort exams 12, 16, 19, 21, 23, 24, 26, 29 and Offspring Cohort exams 1 - 8. 
network2 = network_c2[(network_c2$SPELLBEGIN == 1),]

## CHD outcome
#CHD.whole_c2 = read.table("data/CVD/phs000007.v29.pht003316.v6.p10.c2.vr_survcvd_2014_a_1023s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)

## only consider chd incidence during Exam 11 
before.study.chd = (CHD.whole_c2$chddate <= 0 & CHD.whole_c2$chd == 1)
during11.chd = (CHD.whole_c2$chddate >= 365*(69-48) & CHD.whole_c2$chddate <= 365*(73-48+1) & CHD.whole_c2$chd == 1)

CHD_c2 = cbind(CHD.whole_c2$shareid, during11.chd, before.study.chd)
colnames(CHD_c2)[1] = "shareid"
CHD_c2 = as.data.frame(CHD_c2)

## ex0_11
CHD_ex0_11 = CHD_c2[CHD_c2$shareid %in% pheno_c2_ex0_11$shareid, ]

##  age/sex
#age.whole_c2 = read.table("data/age/phs000007.v29.pht003099.v4.p10.c2.vr_dates_2014_a_0912s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
age.whole_c2 = as.data.frame(age.whole_c2)
age_c2 = as.data.frame(cbind(age.whole_c2$shareid, age.whole_c2$sex, age.whole_c2$age1))
colnames(age_c2) = c("shareid", "sex", "age1")
match_ex0_11 = age_c2[age_c2$shareid %in% pheno_c2_ex0_11$shareid, ]
agesex_ex0_11 = cbind(match_ex0_11[,1:2], match_ex0_11[,3] + 20)
colnames(agesex_ex0_11)[3] = c("age11")

Total = cbind(agesex_ex0_11, CHD_ex0_11[,2:3], pheno_c2_ex0_11$FD44)
names(Total) = c("shareid", "sex", "age11", "CHD", "CHD history", "HDL")

## male
Total.male = Total[Total$sex == 1,] 
Adj = matrix(0, nrow = nrow(Total.male), ncol = nrow(Total.male))
ids = Total.male$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network2[which(network2$shareid == ids[i]), 34])
  friends2 = which(ids %in% network2[which(network2$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}
male.chd = make.permute.moran(Adj, Total.male$CHD, 500)
male.chd.n = nrow(Adj)

Total.male.nona = na.omit(Total.male)
Adj = matrix(0, nrow = nrow(Total.male.nona), ncol = nrow(Total.male.nona))
ids = Total.male.nona$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network2[which(network2$shareid == ids[i]), 34])
  friends2 = which(ids %in% network2[which(network2$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}
result = glm(CHD ~ HDL, family=binomial(link='logit'), data = Total.male.nona)
male.resi = make.permute.moran(Adj, result$residual, 500)
male.hdl = make.permute.moran(Adj, Total.male.nona$HDL, 500)
male.resi.n = nrow(Adj)

## female
Total.female = Total[Total$sex == 2,] 
Adj = matrix(0, nrow = nrow(Total.female), ncol = nrow(Total.female))
ids = Total.female$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network2[which(network2$shareid == ids[i]), 34])
  friends2 = which(ids %in% network2[which(network2$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

female.chd = make.permute.moran(Adj, Total.female$CHD, 500)
female.chd.n = nrow(Adj)

Total.female.nona = na.omit(Total.female)
Adj = matrix(0, nrow = nrow(Total.female.nona), ncol = nrow(Total.female.nona))
ids = Total.female.nona$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network2[which(network2$shareid == ids[i]), 34])
  friends2 = which(ids %in% network2[which(network2$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}
result = glm(CHD ~ HDL, family=binomial(link='logit'), data = Total.female.nona)
female.hdl = make.permute.moran(Adj, Total.female.nona$HDL, 500)
female.resi = make.permute.moran(Adj, result$residual, 500)
female.resi.n = nrow(Adj)

## make a table
tab = matrix(0, nrow = 6, ncol = 4)
colnames(tab) = c("Sex", "n", "Moran's I", "P-value")
tab[1,] = c("Male", male.chd.n,
            formatC(male.chd[1], 2, format = "f"), 
            formatC(male.chd[3], 3, format = "f"))
tab[2,] = c("Female", female.chd.n,
            formatC(female.chd[1], 2, format = "f"), 
            formatC(female.chd[3], 3, format = "f"))
tab[3,] = c("Male", male.resi.n,
            formatC(male.hdl[1], 2, format = "f"), 
            formatC(male.hdl[3], 3, format = "f"))
tab[4,] = c("Female", female.resi.n,
            formatC(female.hdl[1], 2, format = "f"), 
            formatC(female.hdl[3], 3, format = "f"))
tab[5,] = c("Male", male.resi.n,
            formatC(male.resi[1], 2, format = "f"), 
            formatC(male.resi[3], 3, format = "f"))
tab[6,] = c("Female", female.resi.n,
            formatC(female.resi[1], 2, format = "f"), 
            formatC(female.resi[3], 3, format = "f"))
print(xtable(tab))