library(xtable)
source("MoranI.R")

## read data
#pheno_c1_ex0_15 = read.table("data/ex0_15/phs000007.v29.pht000017.v3.p10.c1.ex0_15s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#pheno_c1_ex0_16 = read.table("data/ex0_16/phs000007.v29.pht000018.v4.p10.c1.ex0_16s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#pheno_c1_ex1_1 = read.table("data/ex1_1/phs000007.v29.pht000030.v7.p10.c1.ex1_1s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#pheno_c1_ex1_2 = read.table("data/ex1_2/phs000007.v29.pht000031.v7.p10.c1.ex1_2s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)

#network_c1 = read.table("data/phs000153.v9.pht000836.v8.p8.c1.vr_sntwk_2008_m_0641s.DS-SN-IRB-MDS.txt", sep = "\t", header = TRUE)

network = network_c1[ (network_c1$SPELLBEGIN <= 12*12) & (network_c1$SPELLEND > 9*12),  ]

##CVD outcome
#CVD.whole_c1 = read.table("data/CVD/phs000007.v29.pht003316.v6.p10.c1.vr_survcvd_2014_a_1023s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)

before16.chd = (CVD.whole_c1$chddate < 11315) & (CVD.whole_c1$chd == 1)
during16.chd = (CVD.whole_c1$chddate >= 11315) & (CVD.whole_c1$chddate < 12775) & (CVD.whole_c1$chd == 1)

before16.chf = (CVD.whole_c1$chfdate < 11315) & (CVD.whole_c1$chf == 1)
during16.chf = (CVD.whole_c1$chfdate >= 11315) & (CVD.whole_c1$chfdate < 12775) & (CVD.whole_c1$chf == 1)

before16.cvd = (CVD.whole_c1$cvddate < 11315) & (CVD.whole_c1$cvd == 1)
during16.cvd = (CVD.whole_c1$cvddate >= 11315) & (CVD.whole_c1$cvddate < 12775) & (CVD.whole_c1$cvd == 1)

CVD_c1 = cbind(CVD.whole_c1$shareid, before16.chd, before16.chf, before16.cvd,
                     during16.chd, during16.chf, during16.cvd)
colnames(CVD_c1)[1] = "shareid"
CVD_c1 = as.data.frame(CVD_c1)
## before16 < 365*31 = 11315
## during16 (365*31  = 11315)  & (< 365*35 = 12775) 
## if sum(before16*) > 1 exclude the sample, sum(during16*) == 1, Y = 1.

match_ex0_16 = CVD_c1[CVD_c1$shareid %in% pheno_c1_ex0_16$shareid, ]
CVD_ex0_16 = as.data.frame(cbind(match_ex0_16, as.numeric(rowSums(match_ex0_16[2:4]) == 0), 
                   match_ex0_16[,7] ))
colnames(CVD_ex0_16)[8:9] = c("eligible", "event")

match_ex1_2 = CVD_c1[CVD_c1$shareid %in% pheno_c1_ex1_2$shareid, ]
CVD_ex1_2 = as.data.frame(cbind(match_ex1_2, as.numeric(rowSums(match_ex1_2[2:4]) == 0), 
                   match_ex1_2[,7] ))
colnames(CVD_ex1_2)[8:9] = c("eligible", "event")


## age/sex
#age.whole_c1 = read.table("data/age/phs000007.v29.pht003099.v4.p10.c1.vr_dates_2014_a_0912s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
age.whole_c1 = as.data.frame(age.whole_c1)
age_c1 = as.data.frame(cbind(age.whole_c1$shareid, age.whole_c1$sex, age.whole_c1$age1))
colnames(age_c1) = c("shareid", "sex", "age1")

# ex0_16
match_ex0_16 = age_c1[age_c1$shareid %in% pheno_c1_ex0_16$shareid, ]
agesex_ex0_16 = cbind(match_ex0_16[,1:2], match_ex0_16[,3] + 31,  match_ex0_16[,3] + 31 >= 40)
colnames(agesex_ex0_16)[3:4] = c("age16", ">=40")


# ex1_2
match_ex1_2 = age_c1[age_c1$shareid %in% pheno_c1_ex1_2$shareid, ]
agesex_ex1_2 = cbind(match_ex1_2[,1:2], match_ex1_2[,3] + 8,  match_ex1_2[,3] + 8 >= 40)
colnames(agesex_ex1_2)[3:4] = c("age16", ">=40")

## blood pressure (shareid, diastolic, systolic, pulse)

#pheno_c1_ex0_15 = read.table("data/ex0_15/phs000007.v29.pht000017.v3.p10.c1.ex0_15s.HMB-IRB-MDS.txt", 
#                             sep = "\t", header = TRUE)
#pheno_c1_ex0_16 =  read.table("data/ex0_16/phs000007.v29.pht000018.v4.p10.c1.ex0_16s.HMB-IRB-MDS.txt",
#                           sep = "\t", header = TRUE)

# ex0_16
diastolic_ex0_16 = cbind(pheno_c1_ex0_16$FI22, pheno_c1_ex0_16$FI24, pheno_c1_ex0_16$FI26)
diastolic = rowMeans(diastolic_ex0_16, na.rm = TRUE)
systolic_ex0_16 = cbind(pheno_c1_ex0_16$FI21, pheno_c1_ex0_16$FI23, pheno_c1_ex0_16$FI25)
systolic = rowMeans(systolic_ex0_16, na.rm = TRUE)
pulse = systolic - diastolic                   
blood_ex0_16 = as.data.frame(cbind(pheno_c1_ex0_16[,colnames(pheno_c1_ex0_16) == "shareid"],
                    diastolic, systolic, pulse))
colnames(blood_ex0_16)[1] = c("shareid")

# ex1_2
diastolic_ex1_2 = cbind(pheno_c1_ex1_2$B23, pheno_c1_ex1_2$B25, pheno_c1_ex1_2$B27)
diastolic = rowMeans(diastolic_ex1_2, na.rm = TRUE)
systolic_ex1_2 = cbind(pheno_c1_ex1_2$B22, pheno_c1_ex1_2$B24, pheno_c1_ex1_2$B26)
systolic = rowMeans(systolic_ex1_2, na.rm = TRUE)
pulse = systolic - diastolic                   
blood_ex1_2 = as.data.frame(cbind(pheno_c1_ex1_2[,colnames(pheno_c1_ex1_2) == "shareid"],
                     diastolic, systolic, pulse))
colnames(blood_ex1_2)[1] = c("shareid")

## anti-hyptertensive treatment, binary (shareid, hyper)
# ex0_16
treat.hyper_ex0_16 = pheno_c1_ex0_16$FI187
treat.hyper = ifelse(treat.hyper_ex0_16 == 2, 0, treat.hyper_ex0_16)
hyper_ex0_16 = as.data.frame(cbind(pheno_c1_ex0_16$shareid, treat.hyper))
colnames(hyper_ex0_16)[1] = "shareid"
# ex1_2
treat.hyper_ex1_2 = pheno_c1_ex1_2$B295
treat.hyper = ifelse(treat.hyper_ex1_2 == 2, 0, treat.hyper_ex1_2)
hyper_ex1_2 = as.data.frame(cbind(pheno_c1_ex1_2$shareid, treat.hyper))
colnames(hyper_ex1_2)[1] = "shareid"


## no.of cigarettes smoked per day (shareid, cigar)
# ex0_15
cigar_ex0_15 = cbind(pheno_c1_ex0_15$shareid, pheno_c1_ex0_15$FH102)
colnames(cigar_ex0_15) = c("shareid", "FH102")
cigar_ex0_15 = as.data.frame(cigar_ex0_15)
cigar = c()
for(i in 1:length(pheno_c1_ex0_16$shareid)){
  if(length(which(cigar_ex0_15$shareid == pheno_c1_ex0_16$shareid[i])) == 1){
    cigar[i] = c(cigar_ex0_15$FH102[which(cigar_ex0_15$shareid == pheno_c1_ex0_16$shareid[i])] )
  }else{
    cigar[i] = NA
  }
}
cigar = ifelse(cigar == 88, 0, cigar)
cigar_ex0_16 = as.data.frame(cbind(pheno_c1_ex0_16$shareid, cigar))
colnames(cigar_ex0_16)[1] = "shareid" # has four NA
# ex1_2
cigar_ex1_2 = as.data.frame(cbind(pheno_c1_ex1_2$shareid, pheno_c1_ex1_2$B377))
colnames(cigar_ex1_2) = c("shareid", "cigar")


## Total:HDL cholesterol (shareid, total, HDL, ratio)
# ex0_15
chole_ex0_15 = as.data.frame(cbind(pheno_c1_ex0_15$shareid, pheno_c1_ex0_15$FH332, pheno_c1_ex0_15$FH333))
colnames(chole_ex0_15) = c("shareid", "total", "HDL")
chole = matrix(0, nrow = length(pheno_c1_ex0_16$shareid), ncol = 4)
for(i in 1:length(pheno_c1_ex0_16$shareid)){
  if(length(which(chole_ex0_15$shareid == pheno_c1_ex0_16$shareid[i])) == 1){
    chole[i, 1:3] = as.numeric(chole_ex0_15[which(chole_ex0_15$shareid == pheno_c1_ex0_16$shareid[i]),1:3])
    chole[i, 4] =  chole[i, 2] / chole[i, 3]
  }else{
    chole[i, 1:4] = c(pheno_c1_ex0_16$shareid[i], NA, NA, NA)
  }
}
chole_ex0_16 = as.data.frame(chole)
colnames(chole_ex0_16) = c("shareid", "total", "HDL", "ratio")

# ex1_1
chole_ex1_1 = as.data.frame(cbind(pheno_c1_ex1_1$shareid, pheno_c1_ex1_1$A9, pheno_c1_ex1_1$A10))
colnames(chole_ex1_1) = c("shareid", "total", "HDL")
chole = matrix(0, nrow = length(pheno_c1_ex1_2$shareid), ncol = 4)

for(i in 1:length(pheno_c1_ex1_2$shareid)){
  if(length(which(chole_ex1_1$shareid == pheno_c1_ex1_2$shareid[i])) == 1){
    chole[i, 1:3] = as.numeric(chole_ex1_1[which(chole_ex1_1$shareid == pheno_c1_ex1_2$shareid[i]),1:3])
    chole[i, 4] =  chole[i, 2] / chole[i, 3]
  }else{
    chole[i, 1:4] = c(pheno_c1_ex1_2$shareid[i], NA, NA, NA)
  }
}
chole_ex1_2 = as.data.frame(chole)
colnames(chole_ex1_2) = c("shareid", "total", "HDL", "ratio")


## diabetes 
# (shareid, interim history of insulin, oral hypoglycenmic agents, at 200, at 140)

# ex0_16
treatment = cbind(pheno_c1_ex0_16$shareid, pheno_c1_ex0_16$FI53, pheno_c1_ex0_16$FI54)
treatment[,2] = ifelse(treatment[,2] !=4, 0, 1)
treatment[,3] = ifelse(treatment[,3] !=1, 0, 1)
#dbt.origin_c1 = read.table("data/Dia/phs000007.v29.pht000040.v4.p10.c1.dbt0_27s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
dbt.origin_c1 = as.data.frame(dbt.origin_c1)
dbt_ori_c1 = as.data.frame(cbind(dbt.origin_c1$shareid, dbt.origin_c1$dm200_hx16, dbt.origin_c1$dm140_curr16))
match_ex0_16 = dbt_ori_c1[dbt.origin_c1$shareid %in% pheno_c1_ex0_16$shareid, ]
dbt_ex0_16 = cbind(treatment, match_ex0_16[,2:3])
dbt_status = ifelse(rowSums(dbt_ex0_16[,2:4], na.rm = FALSE) == 0, 0, 1)
dbt_ex0_16 = as.data.frame(cbind(dbt_ex0_16, dbt_status))
colnames(dbt_ex0_16) = c("shareid", "treat.insulin","treat.oral","at.200", "at.140", "dbt.status")

# ex1_2
treatment = cbind(pheno_c1_ex1_2$shareid, pheno_c1_ex1_2$B68, pheno_c1_ex1_2$B65)
treatment[,2] = ifelse(treatment[,2] == 2 | treatment[,2] == 1 , 1, 0)
treatment[,3] = ifelse(treatment[,3] == 1, 1, 0)
#dbt.off_c1 = read.table("data/Dia/phs000007.v29.pht000041.v6.p10.c1.dbt1_7s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
dbt.off_c1 = as.data.frame(dbt.off_c1)
dbt.off_c1 = as.data.frame(cbind(dbt.off_c1$shareid, dbt.off_c1$DIAB2, dbt.off_c1$DMRX2))
names(dbt.off_c1) = c("shareid", "DIAB2", "DMRX2")
match_ex1_2 = dbt.off_c1[dbt.off_c1$shareid %in% pheno_c1_ex1_2$shareid, ]
dbt_ex1_2 = cbind(treatment, match_ex1_2[,2:3])
dbt_status = ifelse(rowSums(dbt_ex1_2[,2:4], na.rm = FALSE) == 0, 0, 1)
dbt_ex1_2 = as.data.frame(cbind(dbt_ex1_2, dbt_status))
colnames(dbt_ex1_2) = c("shareid", "treat.insulin", "treat.oral", "mellitus", "treatment", "dbt.status")


## BMI (shareid, weight(pound), height(inch), BMI(kg/m^2) )
BMI_ex0_16 = cbind(pheno_c1_ex0_16$shareid, pheno_c1_ex0_16$FI14, pheno_c1_ex0_16$FI15)
weight.kg = (pheno_c1_ex0_16$FI14) / 2.2
height.m  = (pheno_c1_ex0_16$FI15) / 39.37
BMI = weight.kg / (height.m)^2
BMI_ex0_16 = as.data.frame(cbind(BMI_ex0_16, BMI))
colnames(BMI_ex0_16) = c("shareid", "weight(pound)", "height(inch)", "BMI(kg/m^2)")
# ex1_2
weight.pound = (pheno_c1_ex1_2$B13)*2.2
height.inch  = (pheno_c1_ex1_2$B14/100)*39.37
BMI = pheno_c1_ex1_2$B13 / (pheno_c1_ex1_2$B14/100)^2
BMI_ex1_2 = as.data.frame(cbind(pheno_c1_ex1_2$shareid, weight.pound, height.inch, BMI))
colnames(BMI_ex1_2) = c("shareid", "weight(pound)", "height(inch)", "BMI(kg/m^2)")


## eletrocardiographic evidence of definite left ventricular hypertrophy
# ex0_16
LVH_ex0_16 = as.data.frame(cbind(pheno_c1_ex0_16$shareid, pheno_c1_ex0_16$FI182))
LVH_ex0_16[,2] = ifelse(LVH_ex0_16[,2] ==4, 1, 0)
colnames(LVH_ex0_16) = c("shareid", "ECG")
# ex1_2
LVH_ex1_2 = as.data.frame(cbind(pheno_c1_ex1_2$shareid, pheno_c1_ex1_2$B275))
LVH_ex1_2[,2] = ifelse(LVH_ex1_2[,2] == 1, 1, 0)
colnames(LVH_ex1_2) = c("shareid", "ECG")


## left ventricular mass (g/m) (sharid, LVID, VST, PWT, mass, adjusted)
# left ventricular mass(grams) = 1.04[(LVID + VST + PWT)^3 - (LVID)^3] - 13.6
# LVIDd and LVIDs – Left ventricular internal diameter end diastole and end systole
# at end-diastole
# ex0_16
mass_ex0_16 = as.data.frame(cbind(pheno_c1_ex0_16$shareid, pheno_c1_ex0_16$FI306,
                                  pheno_c1_ex0_16$FI290, pheno_c1_ex0_16$FI301))
colnames(mass_ex0_16) = c("shareid", "LVID", "VST", "PWT")
mass = 1.04*( (mass_ex0_16$LVID + mass_ex0_16$VST + mass_ex0_16$PWT )^3 - (mass_ex0_16$LVID)^3) - 13.6
mass = mass/1000
height.m = BMI_ex0_16$`height(inch)`/ 39.37
adjusted.mass = mass/height.m
mass_ex0_16 = as.data.frame(cbind(mass_ex0_16, mass, adjusted.mass))

# ex1_2
mass_ex1_2 = as.data.frame(cbind(pheno_c1_ex1_2$shareid, pheno_c1_ex1_2$B448,
                                  pheno_c1_ex1_2$B432, pheno_c1_ex1_2$B443))
colnames(mass_ex1_2) = c("shareid", "LVID", "VST", "PWT")
mass = 1.04*( (mass_ex1_2$LVID + mass_ex1_2$VST + mass_ex1_2$PWT )^3 - (mass_ex1_2$LVID)^3) - 13.6
mass = mass/1000
height.m = BMI_ex1_2$`height(inch)`/ 39.37
adjusted.mass = mass/height.m
mass_ex1_2 = as.data.frame(cbind(mass_ex1_2, mass, adjusted.mass))


## merge ex0_16 and ex1_2
Ex0_16 = as.data.frame(cbind(CVD_ex0_16$shareid, CVD_ex0_16$before16.cvd, CVD_ex0_16$during16.cvd,
                             CVD_ex0_16$eligible, CVD_ex0_16$event,
      agesex_ex0_16$age16, agesex_ex0_16$'>=40', agesex_ex0_16$sex, blood_ex0_16$diastolic,
      blood_ex0_16$pulse, hyper_ex0_16$treat.hyper, cigar_ex0_16$cigar,
      chole_ex0_16$ratio, dbt_ex0_16$dbt.status, BMI_ex0_16$'BMI(kg/m^2)',
      LVH_ex0_16$ECG, mass_ex0_16$'adjusted.mass'))
Ex1_2= as.data.frame(cbind(CVD_ex1_2$shareid, CVD_ex1_2$before16.cvd, CVD_ex1_2$during16.cvd,
              CVD_ex1_2$eligible, CVD_ex1_2$event,
              agesex_ex1_2$age16, agesex_ex1_2$'>=40', agesex_ex1_2$sex, blood_ex1_2$diastolic,
              blood_ex1_2$pulse, hyper_ex1_2$treat.hyper, cigar_ex1_2$cigar,
              chole_ex1_2$ratio, dbt_ex1_2$dbt.status, BMI_ex1_2$'BMI(kg/m^2)',
              LVH_ex1_2$ECG, mass_ex1_2$'adjusted.mass'))
Total = as.data.frame(rbind(Ex0_16, Ex1_2))
colnames(Total) = c("shareid", "before16.cvd", "during16.cvd",
                    "eligible", "event", "age", ">=40", "sex", 
                    "diastolic", "pulse", "treat.hyper", "cigar",
                    "ratio", "dbt.status", "BMI", "ECG", "adjusted.mass")

Total.old = Total[Total$'>=40' == 1, ] 
Total.old.eli = Total.old[Total.old$eligible == 1, ] 
Total.old.eli = na.omit(Total.old.eli) # 1181

result = glm(event ~ sex + age + diastolic + pulse + treat.hyper +
                    cigar + ratio + dbt.status + BMI + ECG + adjusted.mass, family=binomial(link='logit'), 
                  data = Total.old.eli)

residual = result$residuals

Adj = matrix(0, nrow = length(residual), ncol = length(residual))

ids = Total.old.eli$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

## male
Total.male = Total[Total$sex == 1,] 
Total.male.old = Total.male[Total.male$'>=40' == 1,  ] 
Total.male.old.eli  = Total.male.old[Total.male.old$before16.cvd == 0, ] 
Total.male.old.eli = na.omit(Total.male.old.eli)
male.result = glm(during16.cvd ~ age + diastolic + pulse + treat.hyper +
                      cigar + ratio + dbt.status + BMI + ECG + adjusted.mass, family=binomial(link='logit'), 
                    data = Total.male.old.eli)
male.residual = male.result$residuals
Adj = matrix(0, nrow = length(male.residual), ncol = length(male.residual))
ids = Total.male.old.eli$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i] = 1
}


male.cvd = make.permute.moran(Adj, Total.male.old.eli$during16.cvd, 500)
male.resi = make.permute.moran(Adj, male.residual, 500)
male.lvm = make.permute.moran(Adj, Total.male.old.eli$adjusted.mass, 500)


## famale
Total.female = Total[Total$sex == 2, ] 
Total.female.old = Total.female[Total.female$'>=40' == 1,  ] 
Total.female.old.eli = Total.female.old[Total.female.old$before16.cvd == 0, ] 
Total.female.old.eli = na.omit(Total.female.old.eli)
female.result = glm(during16.cvd ~ age + diastolic + pulse + treat.hyper +
      cigar + ratio + dbt.status + BMI + ECG + adjusted.mass, family=binomial(link='logit'), 
    data = Total.female.old.eli)
female.residual = female.result$residuals
Adj = matrix(0, nrow = length(female.residual), ncol = length(female.residual))
ids = Total.female.old.eli$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}

female.cvd = make.permute.moran(Adj, Total.female.old.eli$during16.cvd, 500)
female.resi = make.permute.moran(Adj, female.residual, 500)
female.lvm = make.permute.moran(Adj, Total.female.old.eli$adjusted.mass, 500)


###
tab = matrix(0, nrow = 6, ncol = 4)
colnames(tab) = c("Sex", "Y", "Moran's I", "P-value")
tab[,1] = c("Male", "Female", "Male", "Female", "Male", "Female")
tab[,2] = c(rep("Incidence of CVD",2), rep("LVM",2), rep("Residuals",2))
aa = c(male.cvd[1], female.cvd[1], male.lvm[1], female.lvm[1], male.resi[1], female.resi[1])
bb = c(male.cvd[3], female.cvd[3], male.lvm[1], female.lvm[1], male.resi[3], female.resi[3])
tab[,3] = formatC(aa, 2, format = "f")
tab[,4] = formatC(bb, 3, format = "f")
print(xtable(tab))