library(xtable)
source("MoranI.R")

## read data
#pheno_c2_ex0_17 = read.table("data/phs000007.v29.pht000019.v3.p10.c2.ex0_17s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
#network_c2 = read.table("data/phs000153_c2.txt", sep = "\t", header = TRUE)

# The time interval includes Original Cohort exams 12, 16, 19, 21, 23, 24, 26, 29 and Offspring Cohort exams 1 - 8. 
network = network_c2[ (network_c2$SPELLBEGIN <= 12*14) & (network_c2$SPELLEND > 12*10),  ]

## stroke outcome 
#stroke.whole_c2 = read.table("data/phs000007.v29.pht006023.v1.p10.c2.vr_survstk_2014_a_1031s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
stroke_ex0_17 = stroke.whole_c2[stroke.whole_c2$shareid %in% pheno_c2_ex0_17$shareid, ]
## before17 < 365*33 = 12045
## during17 ( >= 365*33  = 12045)  & (< 365*37 = 13505) 
before17.stroke = (stroke_ex0_17$strokedate < 12045) & (stroke_ex0_17$stroke == 1)
during17.stroke = (stroke_ex0_17$strokedate >= 12045) & (stroke_ex0_17$strokedate < 13505) & (stroke_ex0_17$stroke == 1)

stroke_c2 = cbind(stroke_ex0_17$shareid, before17.stroke, during17.stroke)
colnames(stroke_c2)[1] = "shareid"
stroke_c2 = as.data.frame(stroke_c2)

## the exclusion of subjects with rheumatic heart disease
RHD_ex0_17 = pheno_c2_ex0_17$FJ324 # 2:maybe, 0:no, 1:yes, .:unknown
# only include who are zero

## ECG : Atrial Fibrillation
AF_ex0_17 = pheno_c2_ex0_17$FJ291 #0:no, 1:yes, .:unknown

## Hypertension
hyper_ex0_17 = pheno_c2_ex0_17$FJ318 # 2:maybe, 0: no, 1:yes, .:unknown
hyper_ex0_17 = ifelse(hyper_ex0_17 != 0, 1, 0)

##  age/sex
#age.whole_c2 = read.table("data/age/phs000007.v29.pht003099.v4.p10.c2.vr_dates_2014_a_0912s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
age.whole_c2 = as.data.frame(age.whole_c2)
age_c2 = as.data.frame(cbind(age.whole_c2$shareid, age.whole_c2$sex, age.whole_c2$age1))
colnames(age_c2) = c("shareid", "sex", "age1")

## match age.sex, age should be 61-95
match_ex0_17 = age_c2[age_c2$shareid %in% pheno_c2_ex0_17$shareid, ]
agesex_ex0_17 = cbind(match_ex0_17[,1:2], match_ex0_17[,3] + 33)
colnames(agesex_ex0_17)[3] = c("age17")

## Coronary Heart Disease
#CHD.whole_c2 = read.table("data/CVD/phs000007.v29.pht003316.v6.p10.c2.vr_survcvd_2014_a_1023s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
# before stroke
CHD_ex0_17 = CHD.whole_c2[CHD.whole_c2$shareid %in% pheno_c2_ex0_17$shareid, ]
CHD_before.stroke = (CHD_ex0_17$chddate <= 13505) & (CHD_ex0_17$chd == 1)

## Coronary Heart failure
#CHF.whole_c2 = read.table("data/CVD/phs000007.v29.pht003316.v6.p10.c2.vr_survcvd_2014_a_1023s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
# before stroke
CHF_ex0_17 = CHD.whole_c2[CHF.whole_c2$shareid %in% pheno_c2_ex0_17$shareid, ]
CHF_before.stroke = (CHF_ex0_17$chfdate <= 13505) & (CHF_ex0_17$chf == 1)

Total = cbind(agesex_ex0_17, stroke_c2[,2:3], 
              CHD_before.stroke, CHF_before.stroke, 
              hyper_ex0_17, RHD_ex0_17 , AF_ex0_17)
colnames(Total) = c("shareid", "sex", "age17", "before17.stroke", "during17.stroke", 
                    "CHD_before.stroke", "CHF_before.stroke", 
                    "hyper", "RHD", "AF")
Total = na.omit(Total) ## N = 1662
## exclude who have had stroke
Total.freestroke = Total[Total$before17.stroke == 0, ]
## exclude who have rheumatic heart disease
Total.freeRHD = Total.freestroke[Total.freestroke$RHD == 0, ]


## male ##
Total.male = Total.freeRHD[Total.freeRHD$sex == 1,] 

male.result.AF = glm(during17.stroke ~ AF  , data = Total.male, family = binomial(link = "logit"))
male.result.whole = glm(during17.stroke ~ AF + hyper + CHD_before.stroke + CHF_before.stroke, 
                        data = Total.male, family = binomial(link = "logit"))


male.residual.AF = male.result.AF$residuals
male.residual.whole = male.result.whole$residuals

Adj = matrix(0, nrow = nrow(Total.male), ncol = nrow(Total.male))
ids = Total.male$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

all.male.stat = MoranI(Adj, Total.male$during17.stroke)
all.male.permute = make.permute.moran(Adj, Total.male$during17.stroke, 500)

all.male.AF.stat = MoranI(Adj, Total.male$AF)
all.male.AF.permute = make.permute.moran(Adj, Total.male$AF, 500)

all.male.whole.stat = MoranI(Adj, male.residual.whole)
all.male.whole.permute = make.permute.moran(Adj, male.residual.whole, 500)

## 60-69 years
Total.male.60 = Total.male[Total.male$age17 >= 60 & Total.male$age < 70,] 

male.60.result.AF = glm(during17.stroke ~ AF  , data = Total.male.60, family = binomial(link = "logit"))
male.60.result.whole = glm(during17.stroke ~ AF + hyper + CHD_before.stroke + CHF_before.stroke, 
                        data = Total.male.60, family = binomial(link = "logit"))

male.60.residual.AF = male.60.result.AF$residuals
male.60.residual.whole = male.60.result.whole$residuals

Adj = matrix(0, nrow = nrow(Total.male.60), ncol = nrow(Total.male.60))
ids = Total.male.60$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

all.male.60.stat = MoranI(Adj, Total.male.60$during17.stroke)
all.male.60.permute = make.permute.moran(Adj, Total.male.60$during17.stroke, 500)

all.male.60.AF.stat = MoranI(Adj, Total.male.60$AF)
all.male.60.AF.permute = make.permute.moran(Adj, Total.male.60$AF, 500)

all.male.60.whole.stat = MoranI(Adj, male.60.residual.whole)
all.male.60.whole.permute = make.permute.moran(Adj, male.60.residual.whole, 500)

## 70-79 years
Total.male.70 = Total.male[Total.male$age17 >= 70 & Total.male$age < 80,] 

male.70.result.AF = glm(during17.stroke ~ AF  , data = Total.male.70, family = binomial(link = "logit"))
male.70.result.whole = glm(during17.stroke ~ AF + hyper + CHD_before.stroke + CHF_before.stroke, 
                           data = Total.male.70, family = binomial(link = "logit"))

male.70.residual.AF = male.70.result.AF$residuals
male.70.residual.whole = male.70.result.whole$residuals

Adj = matrix(0, nrow = nrow(Total.male.70), ncol = nrow(Total.male.70))
ids = Total.male.70$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

all.male.70.stat = MoranI(Adj, Total.male.70$during17.stroke)
all.male.70.permute = make.permute.moran(Adj, Total.male.70$during17.stroke, 500)

all.male.70.AF.stat = MoranI(Adj, Total.male.70$AF)
all.male.70.AF.permute = make.permute.moran(Adj, Total.male.70$AF, 500)

all.male.70.whole.stat = MoranI(Adj, male.70.residual.whole)
all.male.70.whole.permute = make.permute.moran(Adj, male.70.residual.whole, 500)


## 80-89 years
Total.male.80 = Total.male[Total.male$age17 >= 80 & Total.male$age < 90,] 

male.80.result.AF = glm(during17.stroke ~ AF  , data = Total.male.80, family = binomial(link = "logit"))
male.80.result.whole = glm(during17.stroke ~ AF + hyper + CHD_before.stroke + CHF_before.stroke, 
                           data = Total.male.80, family = binomial(link = "logit"))

male.80.residual.AF = male.80.result.AF$residuals
male.80.residual.whole = male.80.result.whole$residuals

Adj = matrix(0, nrow = nrow(Total.male.80), ncol = nrow(Total.male.80))
ids = Total.male.80$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

all.male.80.stat = MoranI(Adj, Total.male.80$during17.stroke)
all.male.80.permute = make.permute.moran(Adj, Total.male.80$during17.stroke, 500)

all.male.80.AF.stat = MoranI(Adj, Total.male.80$AF)
all.male.80.AF.permute = make.permute.moran(Adj, Total.male.80$AF, 500)

all.male.80.whole.stat = MoranI(Adj, male.80.residual.whole)
all.male.80.whole.permute = make.permute.moran(Adj, male.80.residual.whole, 500)



## female ##
Total.female = Total.freeRHD[Total.freeRHD$sex == 2,] 

female.result.AF = glm(during17.stroke ~ AF  , data = Total.female, family = binomial(link = "logit"))
female.result.whole = glm(during17.stroke ~ AF + hyper + CHD_before.stroke + CHF_before.stroke, 
                        data = Total.female, family = binomial(link = "logit"))

female.residual.AF = female.result.AF$residuals
female.residual.whole = female.result.whole$residuals

Adj = matrix(0, nrow = nrow(Total.female), ncol = nrow(Total.female))
ids = Total.female$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

all.female.stat = MoranI(Adj, Total.female$during17.stroke)
all.female.permute = make.permute.moran(Adj, Total.female$during17.stroke, 500)

all.female.AF.stat = MoranI(Adj, Total.female$AF)
all.female.AF.permute = make.permute.moran(Adj, Total.female$AF, 500)

all.female.whole.stat = MoranI(Adj, female.residual.whole)
all.female.whole.permute = make.permute.moran(Adj, female.residual.whole, 500)

## 60-69 years
Total.female.60 = Total.female[Total.female$age17 >= 60 & Total.female$age < 70,] 

female.60.result.AF = glm(during17.stroke ~ AF  , data = Total.female.60, family = binomial(link = "logit"))
female.60.result.whole = glm(during17.stroke ~ AF + hyper + CHD_before.stroke + CHF_before.stroke, 
                           data = Total.female.60, family = binomial(link = "logit"))

female.60.residual.AF = female.60.result.AF$residuals
female.60.residual.whole = female.60.result.whole$residuals

Adj = matrix(0, nrow = nrow(Total.female.60), ncol = nrow(Total.female.60))
ids = Total.female.60$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

all.female.60.stat = MoranI(Adj, Total.female.60$during17.stroke)
all.female.60.permute = make.permute.moran(Adj, Total.female.60$during17.stroke, 500)

all.female.60.AF.stat = MoranI(Adj, Total.female.60$AF)
all.female.60.AF.permute = make.permute.moran(Adj, Total.female.60$AF, 500)

all.female.60.whole.stat = MoranI(Adj, female.60.residual.whole)
all.female.60.whole.permute = make.permute.moran(Adj, female.60.residual.whole, 500)

### 70-79 years
Total.female.70 = Total.female[Total.female$age17 >= 70 & Total.female$age < 80,] 

female.70.result.AF = glm(during17.stroke ~ AF  , data = Total.female.70, family = binomial(link = "logit"))
female.70.result.whole = glm(during17.stroke ~ AF + hyper + CHD_before.stroke + CHF_before.stroke, 
                           data = Total.female.70, family = binomial(link = "logit"))

female.70.residual.AF = female.70.result.AF$residuals
female.70.residual.whole = female.70.result.whole$residuals

Adj = matrix(0, nrow = nrow(Total.female.70), ncol = nrow(Total.female.70))
ids = Total.female.70$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

all.female.70.stat = MoranI(Adj, Total.female.70$during17.stroke)
all.female.70.permute = make.permute.moran(Adj, Total.female.70$during17.stroke, 500)

all.female.70.AF.stat = MoranI(Adj, Total.female.70$AF)
all.female.70.AF.permute = make.permute.moran(Adj, Total.female.70$AF, 500)

all.female.70.whole.stat = MoranI(Adj, female.70.residual.whole)
all.female.70.whole.permute = make.permute.moran(Adj, female.70.residual.whole, 500)


## 80-89 years
Total.female.80 = Total.female[Total.female$age17 >= 80 & Total.female$age < 90,] 

female.80.result.AF = glm(during17.stroke ~ AF  , data = Total.female.80, family = binomial(link = "logit"))
female.80.result.whole = glm(during17.stroke ~ AF + hyper + CHD_before.stroke + CHF_before.stroke, 
                           data = Total.female.80, family = binomial(link = "logit"))

female.80.residual.AF = female.80.result.AF$residuals
female.80.residual.whole = female.80.result.whole$residuals

Adj = matrix(0, nrow = nrow(Total.female.80), ncol = nrow(Total.female.80))
ids = Total.female.80$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

all.female.80.stat = MoranI(Adj, Total.female.80$during17.stroke)
all.female.80.permute = make.permute.moran(Adj, Total.female.80$during17.stroke, 500)

all.female.80.AF.stat = MoranI(Adj, Total.female.80$AF)
all.female.80.AF.permute = make.permute.moran(Adj, Total.female.80$AF, 500)

all.female.80.whole.stat = MoranI(Adj, female.80.residual.whole)
all.female.80.whole.permute = make.permute.moran(Adj, female.80.residual.whole, 500)


## male table 
male.table = matrix(0, nrow = 3, ncol = 7)
colnames(male.table) = c("Stroke (I)", "Stroke (p-value)", "AF (I)", "AF (p-value)",
                  "Residual(I)", "Stroke (p-value)", "n")
rownames(male.table) = c("60-69", "70-79", "80-89")
male.table[,1] = c(formatC(all.male.60.permute[1], 2, format = "f"), formatC(all.male.70.permute[1], 2, format = "f"),
                   formatC(all.male.80.permute[1], 2, format = "f"))
male.table[,2] = c(formatC(all.male.60.permute[3], 3, format = "f"), formatC(all.male.70.permute[3], 3, format = "f"),
                   formatC(all.male.80.permute[3], 3, format = "f")) 
male.table[,3] = c(formatC(all.male.60.AF.permute[1], 2, format = "f"), formatC(all.male.70.AF.permute[1], 2, format = "f"),
                   formatC(all.male.80.AF.permute[1], 2, format = "f"))
male.table[,4] = c(formatC(all.male.60.AF.permute[3], 3, format = "f"), formatC(all.male.70.AF.permute[3], 3, format = "f"),
                   formatC(all.male.80.AF.permute[3], 3, format = "f")) 
male.table[,5] = c(formatC(all.male.60.whole.permute[1], 2, format = "f"), formatC(all.male.70.whole.permute[1], 2, format = "f"),
                   formatC(all.male.80.whole.permute[1], 2, format = "f"))
male.table[,6] = c(formatC(all.male.60.whole.permute[3], 3, format = "f"), formatC(all.male.70.whole.permute[3], 3, format = "f"),
                   formatC(all.male.80.whole.permute[3], 3, format = "f")) 

male.table[,7] = c(nrow(Total.male.60), nrow(Total.male.70), nrow(Total.male.80))

print(xtable(male.table))

## female table 
female.table = matrix(0, nrow = 3, ncol = 7)
colnames(female.table) = c("Stroke (I)", "Stroke (p-value)", "AF (I)", "AF (p-value)",
                  "Residual(I)", "Stroke (p-value)", "n")
rownames(female.table) = c("60-69", "70-79", "80-89")
female.table[,1] = c(formatC(all.female.60.permute[1], 2, format = "f"), formatC(all.female.70.permute[1], 2, format = "f"),
                   formatC(all.female.80.permute[1], 2, format = "f"))
female.table[,2] = c(formatC(all.female.60.permute[3], 3, format = "f"), formatC(all.female.70.permute[3], 3, format = "f"),
                   formatC(all.female.80.permute[3], 3, format = "f")) 
female.table[,3] = c(formatC(all.female.60.AF.permute[1], 2, format = "f"), formatC(all.female.70.AF.permute[1], 2, format = "f"),
                   formatC(all.female.80.AF.permute[1], 2, format = "f"))
female.table[,4] = c(formatC(all.female.60.AF.permute[3], 3, format = "f"), formatC(all.female.70.AF.permute[3], 3, format = "f"),
                   formatC(all.female.80.AF.permute[3], 3, format = "f")) 
female.table[,5] = c(formatC(all.female.60.whole.permute[1], 2, format = "f"), formatC(all.female.70.whole.permute[1], 2, format = "f"),
                   formatC(all.female.80.whole.permute[1], 2, format = "f"))
female.table[,6] = c(formatC(all.female.60.whole.permute[3], 3, format = "f"), formatC(all.female.70.whole.permute[3], 3, format = "f"),
                   formatC(all.female.80.whole.permute[3], 3, format = "f")) 

female.table[,7] = c(nrow(Total.female.60), nrow(Total.female.70), nrow(Total.female.80))

print(xtable(female.table))