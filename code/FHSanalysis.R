library(igraph)
library(ColorPalette)
library(RColorBrewer)

source("MoranI.R")

### read data from FHS orginial data, which are not available without authorized access.
#pheno_c1_ex0_15 = read.table("data/ex0_15/phs000007.v29.pht000017.v3.p10.c1.ex0_15s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#pheno_c1_ex0_16 = read.table("data/ex0_16/phs000007.v29.pht000018.v4.p10.c1.ex0_16s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#pheno_c1_ex1_1 = read.table("data/ex1_1/phs000007.v29.pht000030.v7.p10.c1.ex1_1s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#pheno_c1_ex1_2 = read.table("data/ex1_2/phs000007.v29.pht000031.v7.p10.c1.ex1_2s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#network_c1 = read.table("data/phs000153.v9.pht000836.v8.p8.c1.vr_sntwk_2008_m_0641s.DS-SN-IRB-MDS.txt", sep = "\t", header = TRUE)
#CVD.whole_c1 = read.table("data/CVD/phs000007.v29.pht003316.v6.p10.c1.vr_survcvd_2014_a_1023s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#age.whole_c1 = read.table("data/age/phs000007.v29.pht003099.v4.p10.c1.vr_dates_2014_a_0912s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
# add time constraint to network data (formation and cessation)
network = network_c1[ (network_c1$SPELLBEGIN <= 12*12) & (network_c1$SPELLEND > 9*12),  ]



# indicator of having coronary heart disease disease before the exam and during the exam
## before16 < 365*31 = 11315
## during16 (365*31  = 11315)  & (< 365*35 = 12775) 
before16.chd = (CVD.whole_c1$chddate < 11315) & (CVD.whole_c1$chd == 1)
during16.chd = (CVD.whole_c1$chddate >= 11315) & (CVD.whole_c1$chddate < 12775) & (CVD.whole_c1$chd == 1)

# indicator of having coronary heart failure before the exam and during the exam
before16.chf = (CVD.whole_c1$chfdate < 11315) & (CVD.whole_c1$chf == 1)
during16.chf = (CVD.whole_c1$chfdate >= 11315) & (CVD.whole_c1$chfdate < 12775) & (CVD.whole_c1$chf == 1)

# indicator of having cardiovascular disease before the exam and during the exam
before16.cvd = (CVD.whole_c1$cvddate < 11315)  & (CVD.whole_c1$cvd == 1)
during16.cvd = (CVD.whole_c1$cvddate >= 11315) & (CVD.whole_c1$cvddate < 12775) & (CVD.whole_c1$cvd == 1)

# coronary heart disease, coronary heart failure, cardiovascular disease
CVD_c1 = cbind(CVD.whole_c1$shareid, before16.chd, before16.chf, before16.cvd,
               during16.chd, during16.chf, during16.cvd)
colnames(CVD_c1)[1] = "shareid"
CVD_c1 = as.data.frame(CVD_c1)


# match CVD data to the orginial cohort exam 16.
match_ex0_16 = CVD_c1[CVD_c1$shareid %in% pheno_c1_ex0_16$shareid, ]
CVD_ex0_16 = as.data.frame(cbind(match_ex0_16$shareid, match_ex0_16$before16.chd,
                           match_ex0_16$before16.chf, match_ex0_16$before16.cvd))
names(CVD_ex0_16) = c("shareid", "chd", "chf", "cvd")

# myocardial infarction
myo_ex0_16 = pheno_c1_ex0_16$FI181 
myo_ex0_16 = ifelse(myo_ex0_16 == 2, 0, myo_ex0_16) # change 'maybe' to 'no'

# pulmonary disease
pulm_ex0_16 = as.data.frame(cbind(pheno_c1_ex0_16$FI68, pheno_c1_ex0_16$FI70))
# if not 'yes' for now, make it 'no'
pulm_ex0_16[,1] = ifelse(pulm_ex0_16[,1] == 2, 0, pulm_ex0_16[,1])
pulm_ex0_16[,2] = ifelse(pulm_ex0_16[,2] != 1, 0, pulm_ex0_16[,2])
colnames(pulm_ex0_16) = c("pulm1", "pulm2")

# diagnosis of valve disease
valve_ex0_16 = cbind(pheno_c1_ex0_16$FI194, pheno_c1_ex0_16$FI195)
valve_ex0_16[,1] = ifelse(valve_ex0_16[,1] == 2, 0, valve_ex0_16[,1])
valve_ex0_16[,2] = ifelse(valve_ex0_16[,2] == 4, 1, 0)
valve_ex0_16 = as.data.frame(valve_ex0_16)
colnames(valve_ex0_16) = c("valve1", "valve2")

# under treatment of hypertension
treat.hyper_ex0_16 = pheno_c1_ex0_16$FI187
treat.hyper = ifelse(treat.hyper_ex0_16 == 2, 0, treat.hyper_ex0_16)
hyper_ex0_16 = as.data.frame(treat.hyper)


# three different measures for diastolic blood pressure
diastolic_ex0_16 = cbind(pheno_c1_ex0_16$FI22, pheno_c1_ex0_16$FI24, pheno_c1_ex0_16$FI26)
# take the average of three
diastolic = rowMeans(diastolic_ex0_16, na.rm = TRUE)

# three different measures for systolic blood pressure
systolic_ex0_16 = cbind(pheno_c1_ex0_16$FI21, pheno_c1_ex0_16$FI23, pheno_c1_ex0_16$FI25)
# take the average of three
systolic = rowMeans(systolic_ex0_16, na.rm = TRUE)
blood_ex0_16 = as.data.frame(cbind(diastolic, systolic))

# indicator for having a diatolic blood pressure higher than 99 mm HG
highdia = ifelse(blood_ex0_16$diastolic >= 90, 1, 0)
# indicator for having a systolic blood pressure higher than 140 mm HG
highsys = ifelse(blood_ex0_16$systolic >= 140, 1, 0)
blood_ex0_16 = as.data.frame(cbind(blood_ex0_16, highdia, highsys))

## combine the indicators for exclusion.
eligibility_ex0_16 = as.data.frame(cbind(CVD_ex0_16, myo_ex0_16, pulm_ex0_16, valve_ex0_16,
                                         hyper_ex0_16,  blood_ex0_16$highdia, blood_ex0_16$highsys ))


# do the above precedure for offspring cohort exam 2
match_ex1_2 = CVD_c1[CVD_c1$shareid %in% pheno_c1_ex1_2$shareid, ]
CVD_ex1_2 =  as.data.frame(cbind(match_ex1_2$shareid, match_ex1_2$before16.chd,
                           match_ex1_2$before16.chf, match_ex1_2$before16.cvd))
names(CVD_ex1_2) = c("shareid", "chd", "chf", "cvd")

# myocardial infarction
myo_ex1_2 = pheno_c1_ex1_2$B274 
myo_ex1_2 = ifelse(myo_ex1_2 == 2, 0, myo_ex1_2)

# pulmonary disease
pulm_ex1_2 = cbind(pheno_c1_ex1_2$B324, pheno_c1_ex1_2$B135, pheno_c1_ex1_2$B137)
pulm_ex1_2[,1] = ifelse(pulm_ex1_2[,1] == 2, 0, pulm_ex1_2[,1])
pulm_ex1_2[,2] = ifelse(pulm_ex1_2[,2] == 2, 0, pulm_ex1_2[,2])
pulm_ex1_2[,3] = ifelse(pulm_ex1_2[,3] == 2 | pulm_ex1_2[,3] == 3 , 0, pulm_ex1_2[,3])
colnames(pulm_ex1_2) = c("pulm1", "pulm2", "pulm3")

# diagnosis of valve disease
valve_ex1_2 = as.data.frame(cbind(pheno_c1_ex1_2$B302, pheno_c1_ex1_2$B303))
valve_ex1_2[,1] = ifelse(valve_ex1_2[,1] == 2, 0, valve_ex1_2[,1])
valve_ex1_2[,2] = ifelse(valve_ex1_2[,2] == 2, 0, valve_ex1_2[,2])
colnames(valve_ex1_2) = c("valve1", "valve2")

# under treatment of hypertension
treat.hyper_ex1_2 = pheno_c1_ex1_2$B295
treat.hyper = ifelse(treat.hyper_ex1_2 == 2, 0, treat.hyper_ex1_2)
hyper_ex1_2 = as.data.frame(treat.hyper)

# three different measures for diastolic blood pressure
diastolic_ex1_2 = cbind(pheno_c1_ex1_2$B23, pheno_c1_ex1_2$B25, pheno_c1_ex1_2$B27)
diastolic = rowMeans(diastolic_ex1_2, na.rm = TRUE)

# three different measures for systolic blood pressure
systolic_ex1_2 = cbind(pheno_c1_ex1_2$B22, pheno_c1_ex1_2$B24, pheno_c1_ex1_2$B26)
systolic = rowMeans(systolic_ex1_2, na.rm = TRUE)

blood_ex1_2 = as.data.frame(cbind(diastolic, systolic))
highdia = ifelse(blood_ex1_2$diastolic >= 90, 1, 0)
highsys = ifelse(blood_ex1_2$systolic >= 140, 1, 0)
blood_ex1_2 = as.data.frame(cbind(blood_ex1_2, highdia, highsys))

eligibility_ex1_2 = as.data.frame(cbind(CVD_ex1_2, myo_ex1_2, pulm_ex1_2, valve_ex1_2,
                                         hyper_ex1_2,  blood_ex1_2$highdia, blood_ex1_2$highsys ))


### age and sex
age_c1 = as.data.frame(cbind(age.whole_c1$shareid, age.whole_c1$sex, age.whole_c1$age1))
colnames(age_c1) = c("shareid", "sex", "age1")

# match to original cohort exam 16
match_ex0_16 = age_c1[age_c1$shareid %in% pheno_c1_ex0_16$shareid, ]
# add 31 years because nearly 31 years passed from original cohort exam 1
agesex_ex0_16 = cbind(match_ex0_16[,1:2], match_ex0_16[,3] + 31)
colnames(agesex_ex0_16)[3] = c("age16")

# match to offspring cohort exam 2
match_ex1_2 = age_c1[age_c1$shareid %in% pheno_c1_ex1_2$shareid, ]
# add 8 years because nearly 8 years passed from offspring cohort exam 1
agesex_ex1_2 = cbind(match_ex1_2[,1:2], match_ex1_2[,3] + 8)
colnames(agesex_ex1_2)[3] = c("age16")

### blood pressure for each of cohort
BP_ex0_16 = cbind(blood_ex0_16$diastolic, blood_ex0_16$systolic)
BP_ex1_2 = cbind(blood_ex1_2$diastolic, blood_ex1_2$systolic)

### BMI for each of cohort
weight.kg = (pheno_c1_ex0_16$FI14) / 2.2
height.m  = (pheno_c1_ex0_16$FI15) / 39.37
BMI = weight.kg / (height.m)^2
BMI_ex0_16 = as.data.frame(cbind(weight.kg, height.m*100, BMI))
colnames(BMI_ex0_16) = c( "weight", "height", "BMI")

weight.pound = (pheno_c1_ex1_2$B13)*2.2
height.inch  = (pheno_c1_ex1_2$B14/100)*39.37
BMI = pheno_c1_ex1_2$B13 / (pheno_c1_ex1_2$B14/100)^2
BMI_ex1_2 = as.data.frame(cbind(pheno_c1_ex1_2$B13, pheno_c1_ex1_2$B14, BMI))
colnames(BMI_ex1_2) = c("weight", "height", "BMI")


### Left ventricular mass (LVM)
# left ventricular mass (g/m) (sharid, LVID, VST, PWT, mass, adjusted)
# left ventricular mass(grams) = 1.04[(LVID + VST + PWT)^3 - (LVID)^3] - 13.6


## original cohort exam 16 
# FI306 : ECHO: LVID-D-STD (MM)
# FI297 : ECHO: IV-SEPT-THICK-STD 
# FI301 : ECHO: LV-POST-WALL-THIC-STD (MM)
mass_ex0_16 = as.data.frame(cbind(pheno_c1_ex0_16$shareid, pheno_c1_ex0_16$FI306,
                                  pheno_c1_ex0_16$FI297, pheno_c1_ex0_16$FI301))
colnames(mass_ex0_16) = c("shareid", "LVID", "VST", "PWT")
mass = 1.04*( (mass_ex0_16$LVID + mass_ex0_16$VST + mass_ex0_16$PWT )^3 - (mass_ex0_16$LVID)^3)/1000
mass = mass - 13.6
height.m = BMI_ex0_16$height / 100
adjusted.mass = mass/height.m
mass_ex0_16 = as.data.frame(cbind(mass_ex0_16, mass, adjusted.mass))


# offspring cohort exam 2
## B448 : ECHO: LVID-D-STD (MM)
## B439 : ECHO: IV-SEPT-THICK-STD (MM)
## B443 : ECHO: LV-POST-WALL-THIC-STD (MM)
mass_ex1_2 = as.data.frame(cbind(pheno_c1_ex1_2$shareid, pheno_c1_ex1_2$B448,
                                 pheno_c1_ex1_2$B439, pheno_c1_ex1_2$B443))
colnames(mass_ex1_2) = c("shareid", "LVID", "VST", "PWT")
mass = 1.04*( (mass_ex1_2$LVID + mass_ex1_2$VST + mass_ex1_2$PWT )^3 - (mass_ex1_2$LVID)^3)/1000
mass = mass - 13.6
height.m = BMI_ex1_2$height / 100
adjusted.mass = mass/height.m
mass_ex1_2 = as.data.frame(cbind(mass_ex1_2, mass, adjusted.mass))


### merge two cohorts
# if at least one exclusin criteria is met, exclusion indicator is one.
eli_ex0_16 = rowSums(eligibility_ex0_16[,-1]) > 0
eli_ex1_2 = rowSums(eligibility_ex1_2[,-1]) > 0

Ex0_16 = as.data.frame(cbind(CVD_ex0_16$shareid, eli_ex0_16,
                             agesex_ex0_16$sex, agesex_ex0_16$age16, 
                             blood_ex0_16$diastolic, blood_ex0_16$systolic,
                             BMI_ex0_16$weight, BMI_ex0_16$height, BMI_ex0_16$BMI, 
                             mass_ex0_16$mass, mass_ex0_16$adjusted.mass))

Ex1_2 = as.data.frame(cbind(CVD_ex1_2$shareid, eli_ex1_2,
                            agesex_ex1_2$sex, agesex_ex1_2$age16, 
                            blood_ex1_2$diastolic, blood_ex1_2$systolic,
                            BMI_ex1_2$weight, BMI_ex1_2$height, BMI_ex1_2$BMI, 
                            mass_ex1_2$mass, mass_ex1_2$adjusted.mass))

colnames(Ex0_16) = c("shareid", "exclude", "sex", "age",  "diastolic", "systolic",
                    "weight", "height", "BMI", "mass", "adjusted.mass")
colnames(Ex1_2) = c("shareid", "exclude", "sex", "age",  "diastolic", "systolic",
                    "weight", "height", "BMI", "mass", "adjusted.mass")

Total = as.data.frame(rbind(Ex0_16, Ex1_2))
colnames(Total) = c("shareid", "exclude", "sex", "age", "diastolic", "systolic",
                    "weight", "height", "BMI", "mass", "adjusted.mass")

Total.eli = Total[Total$exclude == 0, ] # subset of eligible subjects (so-called healty subjects)
Total.eli.obs = na.omit(Total.eli)
Total.male.eli =  Total.eli.obs[Total.eli.obs$sex == 1, ]
Total.female.eli =  Total.eli.obs[Total.eli.obs$sex == 2, ]



Total.male.eli$BMI = ifelse(Total.male.eli$BMI < 23, 1, Total.male.eli$BMI)
Total.male.eli$BMI = ifelse(Total.male.eli$BMI >= 23 & Total.male.eli$BMI < 26, 2, Total.male.eli$BMI)
Total.male.eli$BMI = ifelse(Total.male.eli$BMI >= 26 & Total.male.eli$BMI < 30, 3, Total.male.eli$BMI)
Total.male.eli$BMI = ifelse(Total.male.eli$BMI >= 30, 4, Total.male.eli$BMI)


male.result.BMI = lm(adjusted.mass ~ as.factor(BMI), data = Total.male.eli)
male.result.ageBMI = lm(adjusted.mass ~ age + as.factor(BMI), data = Total.male.eli)
male.result.ageBMIBlood = lm(adjusted.mass ~ age + as.factor(BMI) + systolic, data = Total.male.eli)


male.residual.BMI = male.result.BMI$residuals
male.residual.ageBMI = male.result.ageBMI$residuals
male.residual.ageBMIBlood = male.result.ageBMIBlood$residuals

Adj = matrix(0, nrow = length(male.residual.BMI), ncol = length(male.residual.BMI))
ids = Total.male.eli$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i] = 1
}

#isolate = which(rowSums(Adj) == 0)
Male.mass = MoranI(Adj, Total.male.eli$adjusted.mass)
Male.BMI = MoranI(Adj, Total.male.eli$BMI)
Male.model1 = MoranI(Adj, male.residual.BMI)
Male.model2 = MoranI(Adj, male.residual.ageBMI)
Male.model3 = MoranI(Adj, male.residual.ageBMIBlood)

Male.mass.permute = make.permute.moran(Adj, Total.male.eli$adjusted.mass, 500)
Male.BMI.permute = make.permute.moran(Adj, Total.male.eli$BMI, 500)
Male.age.permute = make.permute.moran(Adj, Total.male.eli$age, 500)
Male.systolic.permute = make.permute.moran(Adj, Total.male.eli$systolic, 500)

Male.model1.permute = make.permute.moran(Adj, male.residual.BMI, 500)
Male.model2.permute = make.permute.moran(Adj, male.residual.ageBMI, 500)
Male.model3.permute = make.permute.moran(Adj, male.residual.ageBMIBlood, 500)

######## famale
Total.female = Total[Total$sex == 2, ] 
Total.female.eli = Total.female[Total.female$exclude == 0,  ] 
Total.female.eli = na.omit(Total.female.eli)

Total.female.eli$BMI = ifelse(Total.female.eli$BMI < 23, 1, Total.female.eli$BMI)
Total.female.eli$BMI = ifelse(Total.female.eli$BMI >= 23 & Total.female.eli$BMI < 26, 2, Total.female.eli$BMI)
Total.female.eli$BMI = ifelse(Total.female.eli$BMI >= 26 & Total.female.eli$BMI < 30, 3, Total.female.eli$BMI)
Total.female.eli$BMI = ifelse(Total.female.eli$BMI >= 30, 4, Total.female.eli$BMI)

female.result.BMI = lm(adjusted.mass ~ as.factor(BMI), data = Total.female.eli)
female.result.ageBMI = lm(adjusted.mass ~ age + as.factor(BMI), data = Total.female.eli)
female.result.ageBMIBlood = lm(adjusted.mass ~ age + as.factor(BMI) + systolic, data = Total.female.eli)

female.residual.BMI = female.result.BMI$residuals
female.residual.ageBMI = female.result.ageBMI$residuals
female.residual.ageBMIBlood = female.result.ageBMIBlood$residuals

Adj = matrix(0, nrow = nrow(Total.female.eli), ncol = nrow(Total.female.eli))
ids = Total.female.eli$shareid

for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network[which(network$shareid == ids[i]), 34])
  friends2 = which(ids %in% network[which(network$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i] = 1
}

Female.Adj = Adj
Female.mass = MoranI(Adj, Total.female.eli$adjusted.mass)
Female.BMI = MoranI(Adj, Total.female.eli$BMI)
Female.model1 = MoranI(Adj, female.residual.BMI)
Female.model2 = MoranI(Adj, female.residual.ageBMI)
Female.model3 = MoranI(Adj, female.residual.ageBMIBlood)

Female.mass.permute = make.permute.moran(Adj, Total.female.eli$adjusted.mass, 500)
Female.BMI.permute = make.permute.moran(Adj, Total.female.eli$BMI, 500)
Female.age.permute = make.permute.moran(Adj, Total.female.eli$age, 500)
Female.systolic.permute = make.permute.moran(Adj, Total.female.eli$systolic, 500)

Female.model1.permute = make.permute.moran(Adj, female.residual.BMI, 500)
Female.model2.permute = make.permute.moran(Adj, female.residual.ageBMI, 500)
Female.model3.permute = make.permute.moran(Adj, female.residual.ageBMIBlood, 500)



#############################################
summary.table = matrix(0, nrow =  7, ncol = 2)
colnames(summary.table) = c("Male (n = 684)", "Female (n = 1003)")
rownames(summary.table) = c("Age (year)", "Weight (kg)",
                            "Height (cm)", "BMI (kg/m^2)", 
                            "Systolic BP (mm Hg)", "LVM (g)",
                            "Adjusted LVM (g/m)") 
summary.table[1,] = c( paste(round(mean(Total.male.eli$age),2), " (", round(sd(Total.male.eli$age),2), ")", sep=""), 
                       paste(round(mean(Total.female.eli$age),2), " (", round(sd(Total.female.eli$age),2), ")", sep=""))
summary.table[2,] = c( paste(round(mean(Total.male.eli$weight),2), " (", round(sd(Total.male.eli$weight),2), ")", sep=""), 
                      paste(round(mean(Total.female.eli$weight),2), " (", round(sd(Total.female.eli$weight),2), ")", sep=""))
summary.table[3,] = c( paste(round(mean(Total.male.eli$height),2), " (", round(sd(Total.male.eli$height),2), ")", sep=""), 
                       paste(round(mean(Total.female.eli$height),2), " (", round(sd(Total.female.eli$height),2), ")", sep=""))
summary.table[4,] = c( paste(round(mean(Total.male.eli$BMI),2), " (", round(sd(Total.male.eli$BMI),2), ")", sep=""), 
                       paste(round(mean(Total.female.eli$BMI),2), " (", round(sd(Total.female.eli$BMI),2), ")", sep=""))
summary.table[5,] = c(paste(round(mean(Total.male.eli$systolic),2), " (", round(sd(Total.male.eli$systolic),2), ")", sep=""), 
                       paste(round(mean(Total.female.eli$systolic),2), " (", round(sd(Total.female.eli$systolic),2), ")", sep=""))
summary.table[6,] = c(paste(round(mean(Total.male.eli$mass),2), " (", round(sd(Total.male.eli$mass),2), ")", sep=""), 
                      paste(round(mean(Total.female.eli$mass),2), " (", round(sd(Total.female.eli$mass),2), ")", sep=""))
summary.table[7,] = c(paste(round(mean(Total.male.eli$adjusted.mass),2), " (", round(sd(Total.male.eli$adjusted.mass),2), ")", sep=""), 
                      paste(round(mean(Total.female.eli$adjusted.mass),2), " (", round(sd(Total.female.eli$adjusted.mass),2), ")", sep=""))

xtable(print(summary.table))

############
Male.table = matrix(0, nrow = 5, ncol = 2)
colnames(Male.table) = c("I_{std} (p-value)", "p-value (permutation)")
rownames(Male.table) = c("Adjusted LVM (g/m)",  "BMI (kg/m^2)",
                         "Residual of Adjusted LVM ~ BMI",
                         "Residual of Adjusted LVM ~ BMI + age",
                         "Residual of Adjusted LVM ~ BMI + age + systolic BP" )
Male.table[1,] = c( paste(round(Male.mass.permute[1],2), " (" , round(Male.mass.permute[2],3) ,")", sep=""),
                    round(Male.mass.permute[3], 3))
Male.table[2,] = c( paste(round(Male.BMI.permute[1],2), " (" , round(Male.BMI.permute[2],3) ,")", sep=""),
                    round(Male.BMI.permute[3], 3))
Male.table[3,] = c( paste(round(Male.model1.permute[1],2), " (" , round(Male.model1.permute[2],3) ,")", sep=""),
                    round(Male.model1.permute[3], 3))
Male.table[4,] = c( paste(round(Male.model2.permute[1],2), " (" , round(Male.model2.permute[2],3) ,")", sep=""),
                    round(Male.model2.permute[3], 3))
Male.table[5,] = c( paste(round(Male.model3.permute[1],2), " (" , round(Male.model3.permute[2],3) ,")", sep=""),
                    round(Male.model3.permute[3], 3))

xtable(print(Male.table))


############
Female.table = matrix(0, nrow = 5, ncol = 2)
colnames(Female.table) = c("I_{std} (p-value)", "p-value (permutation)")
rownames(Female.table) = c("Adjusted LVM (g/m)",  "BMI (kg/m^2)",
                         "Residual of Adjusted LVM ~ BMI",
                         "Residual of Adjusted LVM ~ BMI + age",
                         "Residual of Adjusted LVM ~ BMI + age + systolic BP" )
Female.table[1,] = c( paste(round(Female.mass.permute[1],2), " (" , round(Female.mass.permute[2],3) ,")", sep=""),
                    round(Female.mass.permute[3], 3))
Female.table[2,] = c( paste(round(Female.BMI.permute[1],2), " (" , round(Female.BMI.permute[2],3) ,")", sep=""),
                    round(Female.BMI.permute[3], 3))
Female.table[3,] = c( paste(round(Female.model1.permute[1],2), " (" , round(Female.model1.permute[2],3) ,")", sep=""),
                    round(Female.model1.permute[3], 3))
Female.table[4,] = c( paste(round(Female.model2.permute[1],2), " (" , round(Female.model2.permute[2],3) ,")", sep=""),
                    round(Female.model2.permute[3], 3))
Female.table[5,] = c( paste(round(Female.model3.permute[1],2), " (" , round(Female.model3.permute[2],3) ,")", sep=""),
                    round(Female.model3.permute[3], 3))
xtable(print(Female.table))

###############################################################
summary(male.result.ageBMIBlood)
table(Total.male.eli$BMI)
mean.age = 40; mean.sys = 118
male.estimates = summary(male.result.ageBMIBlood)$coefficients[,1]
male.BMI1 = male.estimates[1] + mean.age*male.estimates[2] + 118*male.estimates[6]
male.BMI2 = male.estimates[1] + mean.age*male.estimates[2] + male.estimates[3] + 118*male.estimates[6]
male.BMI3 = male.estimates[1] + mean.age*male.estimates[2] + male.estimates[4] + 118*male.estimates[6]
male.BMI4 = male.estimates[1] + mean.age*male.estimates[2] + male.estimates[5] + 118*male.estimates[6]


female.estimates = summary(female.result.ageBMIBlood)$coefficients[,1]
female.BMI1 = female.estimates[1] + mean.age*female.estimates[2] + 118*female.estimates[6]
female.BMI2 = female.estimates[1] + mean.age*female.estimates[2] + female.estimates[3] + 118*female.estimates[6]
female.BMI3 = female.estimates[1] + mean.age*female.estimates[2] + female.estimates[4] + 118*female.estimates[6]
female.BMI4 = female.estimates[1] + mean.age*female.estimates[2] + female.estimates[5] + 118*female.estimates[6]


count.mat = matrix(0, nrow = 2, ncol = 4)
colnames(count.mat) = c("<23", "23-25.99", "26-29.99", ">=30")
rownames(count.mat) = c("Male", "Female")
count.mat[1,] = c(male.BMI1, male.BMI2, male.BMI3, male.BMI4)
count.mat[2,] = c(female.BMI1, female.BMI2, female.BMI3, female.BMI4)

female.mat = formatC(summary(female.result.ageBMIBlood)$coefficients, digits = 2, format = "f")
male.mat = formatC(summary(male.result.ageBMIBlood)$coefficients, digits = 2, format = "f")
print(xtable(rbind(female.mat, male.mat)))
################################################################
RdPalette = colorRampPalette(brewer.pal(9,"Reds"))(20)
BlPalette = colorRampPalette(brewer.pal(9,"Blues"))(20)


### Male
Y = male.residual.ageBMIBlood
A = Male.Adj
G = graph.adjacency(A, "undirected")
V(G)$Y = Y
Male.sub = induced.subgraph(G, components(G)$membership == 1)
for(i in 1:20){
  lower = quantile(V(Male.sub)$Y, 0.05*(i-1))
  upper = quantile(V(Male.sub)$Y, 0.05*i)
  V(Male.sub)[V(Male.sub)$Y >= lower & V(Male.sub)$Y < upper]$color = BlPalette[i]
}
V(Male.sub)[which.max(V(Male.sub)$Y)]$color = BlPalette[20]



igraph.options(vertex.size = 3, edge.arrow.size = 0.1,
               vertex.label = NULL)
pdf("../figures/subMale_residual.pdf")
par(mar=c(0,0,3,0), cex.main = 3)
set.seed(123)
plot(Male.sub, layout = layout.fruchterman.reingold,
     vertex.label = "", main = "Male Network")
dev.off()



### Female
Y = female.residual.ageBMIBlood
A = Female.Adj
G = graph.adjacency(A, "undirected")
V(G)$Y = Y
Female.sub = induced.subgraph(G, components(G)$membership == 3)
for(i in 1:20){
  lower = quantile(V(Female.sub)$Y, 0.05*(i-1))
  upper = quantile(V(Female.sub)$Y, 0.05*i)
  V(Female.sub)[V(Female.sub)$Y >= lower & V(Female.sub)$Y < upper]$color = RdPalette[i]
}
V(Female.sub)[which.max(V(Female.sub)$Y)]$color = RdPalette[20]


igraph.options(vertex.size = 3, edge.arrow.size = 0.1,
               vertex.label = NULL)
pdf("../figures/subFemale_residual.pdf")
par(mar=c(0,0,3,0), cex.main = 3)
set.seed(123)
plot(Female.sub, layout = layout.fruchterman.reingold,
     vertex.label = "", main = "Female Network")
dev.off()

