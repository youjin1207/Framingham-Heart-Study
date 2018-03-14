library(survival)
library(igraph)

source("MoranI.R")


## Read phenotype & network data 
#pheno_c2_ex0_18 = read.table("data/phs000007.v29.pht000020.v3.p10.c2.ex0_18s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
#network_c2 = read.table("data/phs000153.v9.pht000836.v8.p8.c2.vr_sntwk_2008_m_0641s.DS-SN-IRB-MDS.txt", sep = "\t", header = TRUE)

## The time interval includes Original Cohort exams 12, 16, 19, 21, 23, 24, 26, 29 and Offspring Cohort exams 1 - 8. 
## approximate network ties existing during 1983 - 1985 (Original Cohort exam 18)
network2 = network_c2[(network_c2$SPELLBEGIN <= 12*(85-71+1)) & (network_c2$SPELLEND > 12*(84-71)), ]

## 0. death data
#surv.whole_c2 = read.table("data/phs000007.v29.pht003335.v5.p10.c2.vr_soe4srv_2014_a_1027s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
#hrv_c2 = read.table("data/phs000007.v29.pht000197.v6.p10.c2.hrv_1987s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
## table(hrv_c2$hrv_exam) : either 3 or 18
hrv_c2_ex0_18 = hrv_c2[hrv_c2$hrv_exam==18,]

# Date of Event (number of days since Exam 1)
# 365*35, 365*29
# exclude subjects who had already died before 18th exam : 3793 dead before 18th exam
surv.ex0_18_c2 = surv.whole_c2[!((surv.whole_c2$EVENT >= 21) & (surv.whole_c2$EVENT <= 29) & (surv.whole_c2$DATE < 365*35) ), ] 
death_ex0_18 = surv.ex0_18_c2[surv.ex0_18_c2$shareid %in% hrv_c2_ex0_18$shareid, ]

death.id = unique(death_ex0_18$shareid) # for each subject
death.time = c(); death.indi = c()
for(i in 1:length(death.id)){
  tmp.table = death_ex0_18[death_ex0_18$shareid == death.id[i],] 
  death.indi[i] = (sum(tmp.table[nrow(tmp.table),]$EVENT %in% c(21:29)) > 0) & 
  (tmp.table[nrow(tmp.table),]$DATE >= 365*35) & (tmp.table[nrow(tmp.table),]$DATE < 365*39)
  death.time[i] = ifelse(death.indi[i] == TRUE, tmp.table[nrow(tmp.table),]$DATE, 365*39)
}

# if dealth.ind == 0, death.time is censored time.
death.data = cbind(death.id, death.time, death.indi)
death.data = as.data.frame(death.data)

## HRV Data ##
hrv.data = hrv_c2_ex0_18[hrv_c2_ex0_18$shareid %in% death.id, ]
hrv.data = cbind(hrv.data$shareid, hrv.data$hrv_sdnn, hrv.data$hrv_pnn50,
                 hrv.data$hrv_rmssd, hrv.data$hrv_vlf_power, hrv.data$hrv_lf_power,
                 hrv.data$hrv_hf_power, hrv.data$hrv_tot_power)
colnames(hrv.data) = c("shareid", "SDNN", "pNN50", "rMSSD", "VLF", "LF", "HF", "TP")
hrv.data = as.data.frame(hrv.data)

## make a comparison in SD of HRVs between ours and those from the original paper.
SD = c(apply(log(hrv.data), 2, sd)[2:8], sd(log(hrv.data[,6]/hrv.data[,7])))
SD.table = matrix(0, nrow = 2, ncol = 8)
colnames(SD.table) = c("lnSDNN", "lnpNN50", "lnr-MSSD", "lnVLF", "lnLF", "lnHF", "lnTP", "lnLF/HF")
rownames(SD.table) = c("Original paper", "Our data")
SD.table[1,] = c(0.33, 1.32, 0.44, 0.76, 0.82, 0.85, 0.73, 0.57)
SD.table[2,] = formatC(SD, digits = 2, format = "f")
print(xtable(SD.table))

## unadjusted association with all-cause mortality (n = 516)
surv.object = Surv(death.data$death.time, death.data$death.indi)
fit.SDNN = coxph(surv.object ~ log(SDNN), data = hrv.data)
fit.pNN50 = coxph(surv.object ~ log(pNN50), data = hrv.data)
fit.rMSSD = coxph(surv.object ~ log(rMSSD), data = hrv.data)
fit.VLF = coxph(surv.object ~ log(VLF), data = hrv.data)
fit.LF = coxph(surv.object ~ log(LF), data = hrv.data)
fit.HF = coxph(surv.object ~ log(HF), data = hrv.data)
fit.TP = coxph(surv.object ~ log(TP), data = hrv.data)
fit.ratio = coxph(surv.object ~ log(LF/HF), data = hrv.data)


Adj = matrix(0, nrow = nrow(hrv.data), ncol = nrow(hrv.data))
ids = hrv.data$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network2[which(network2$shareid == ids[i]), 34])
  friends2 = which(ids %in% network2[which(network2$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}
moran.SDNN = make.permute.moran(Adj, fit.SDNN$residuals, 500)
moran.pNN50 = make.permute.moran(Adj, fit.pNN50$residuals, 500)
moran.rMSSD = make.permute.moran(Adj, fit.rMSSD$residuals, 500)
moran.VLF = make.permute.moran(Adj, fit.VLF$residuals, 500)
moran.LF = make.permute.moran(Adj, fit.LF$residuals, 500)
moran.HF = make.permute.moran(Adj, fit.HF$residuals, 500)
moran.TP = make.permute.moran(Adj, fit.TP$residuals, 500)
moran.ratio = make.permute.moran(Adj, fit.ratio$residuals, 500)

## collect subject-specific data
#  age/sex
#age.whole_c2 = read.table("data/age/phs000007.v29.pht003099.v4.p10.c2.vr_dates_2014_a_0912s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
age.whole_c2 = as.data.frame(age.whole_c2)
age_c2 = as.data.frame(cbind(age.whole_c2$shareid, age.whole_c2$sex, age.whole_c2$age1))
colnames(age_c2) = c("shareid", "sex", "age1")

# match age.sex
## age should be 49-82
match_ex0_18 = age_c2[age_c2$shareid %in% hrv.data$shareid, ]
agesex_ex0_18 = cbind(match_ex0_18[,1:2], match_ex0_18[,3] + 35) # age at exam 18 -- extrapolated from age measured at exam 1
colnames(agesex_ex0_18)[3] = c("age18")
agesex.data = cbind(hrv.data, agesex_ex0_18)

## Age- and sex-adjusted association with all-cause mortality (n=516)
surv.object = Surv(death.data$death.time, death.data$death.indi)
fit.agesex.SDNN = coxph(surv.object ~ log(SDNN) + as.factor(sex) + age18, data = agesex.data)
fit.agesex.pNN50 = coxph(surv.object ~ log(pNN50) + as.factor(sex) + age18, data = agesex.data)
fit.agesex.rMSSD = coxph(surv.object ~ log(rMSSD) + as.factor(sex) + age18, data = agesex.data)
fit.agesex.VLF = coxph(surv.object ~ log(VLF) + as.factor(sex) + age18, data = agesex.data)
fit.agesex.LF = coxph(surv.object ~ log(LF) + as.factor(sex) + age18, data = agesex.data)
fit.agesex.HF = coxph(surv.object ~ log(HF) + as.factor(sex) + age18, data = agesex.data)
fit.agesex.TP = coxph(surv.object ~ log(TP) + as.factor(sex) + age18, data = agesex.data)
fit.agesex.ratio = coxph(surv.object ~ log(LF/HF) + as.factor(sex) + age18, data = agesex.data)

Adj = matrix(0, nrow = nrow(agesex.data), ncol = nrow(agesex.data))
ids = agesex.data$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network2[which(network2$shareid == ids[i]), 34])
  friends2 = which(ids %in% network2[which(network2$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}
moran.agesex.SDNN = make.permute.moran(Adj, fit.agesex.SDNN$residuals, 500)
moran.agesex.pNN50 = make.permute.moran(Adj, fit.agesex.pNN50$residuals, 500)
moran.agesex.rMSSD = make.permute.moran(Adj, fit.agesex.rMSSD$residuals, 500)
moran.agesex.VLF = make.permute.moran(Adj, fit.agesex.VLF$residuals, 500)
moran.agesex.LF = make.permute.moran(Adj, fit.agesex.LF$residuals, 500)
moran.agesex.HF = make.permute.moran(Adj, fit.agesex.HF$residuals, 500)
moran.agesex.TP = make.permute.moran(Adj, fit.agesex.TP$residuals, 500)
moran.agesex.ratio = make.permute.moran(Adj, fit.agesex.ratio$residuals, 500)

## collect other clinical risk factors
pheno.ex0_18 = pheno_c2_ex0_18[pheno_c2_ex0_18$shareid %in% agesex.data$shareid, ]

# history of myocardial infarction
myo = ifelse(pheno.ex0_18$FK392 == 0, 0, 1) # 0 : no
# diuretic use
diuretic = ifelse(pheno.ex0_18$FK86 == 1, 1, 0) # 1: yes, now!
# presence or ventricular premature beats
beats = ifelse(pheno.ex0_18$FK384 == 0, 0, 1) # 0 : no
# history of congestive heart failutre
#CHF.whole_c2 = read.table("data/CVD/phs000007.v29.pht003316.v6.p10.c2.vr_survcvd_2014_a_1023s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
CHF.ex0_18 = CHF.whole_c2[CHF.whole_c2$shareid %in% agesex.data$shareid, ]
# before18 < 365*35 = 12775
CHF.before18 = CHF.ex0_18$chfdate < 12775 & CHF.ex0_18$chfdate > 0 & CHF.ex0_18$chf == 1

risk.factor = cbind(myo, CHF.before18, beats, diuretic)
colnames(risk.factor) = c("myo", "CHF", "beats", "diuretic")
risk.factor = as.data.frame(risk.factor)

all.data = cbind(agesex.data, death.data, risk.factor)
all.data = na.omit(all.data) # four subjects with missing data excluded.

## Age, sex, and clinical risk factors adjusted for association with all-cause mortality
surv.object = Surv(all.data$death.time, all.data$death.indi)
fit.all.SDNN = coxph(surv.object ~ log(SDNN) + as.factor(sex) + age18 + as.factor(myo) + as.factor(CHF) + as.factor(beats) + as.factor(diuretic), data = all.data)
fit.all.pNN50 = coxph(surv.object ~ log(pNN50) + as.factor(sex) + age18 + as.factor(myo) + as.factor(CHF) + as.factor(beats) + as.factor(diuretic), data = all.data)
fit.all.rMSSD = coxph(surv.object ~ log(rMSSD) + as.factor(sex) + age18 + as.factor(myo) + as.factor(CHF) + as.factor(beats) + as.factor(diuretic), data = all.data)
fit.all.VLF = coxph(surv.object ~ log(VLF) + as.factor(sex) + age18 + as.factor(myo) + as.factor(CHF) + as.factor(beats) + as.factor(diuretic), data = all.data)
fit.all.LF = coxph(surv.object ~ log(LF) + as.factor(sex) + age18 + as.factor(myo) + as.factor(CHF) + as.factor(beats) + as.factor(diuretic), data = all.data)
fit.all.HF = coxph(surv.object ~ log(HF) + as.factor(sex) + age18 + as.factor(myo) + as.factor(CHF) + as.factor(beats) + as.factor(diuretic), data = all.data)
fit.all.TP = coxph(surv.object ~ log(TP) + as.factor(sex) + age18 + as.factor(myo) + as.factor(CHF) + as.factor(beats) + as.factor(diuretic), data = all.data)
fit.all.ratio = coxph(surv.object ~ log(LF/HF) + as.factor(sex) + age18 + as.factor(myo) + as.factor(CHF) + as.factor(beats) + as.factor(diuretic), data = all.data)


Adj = matrix(0, nrow = nrow(all.data), ncol = nrow(all.data))
ids = all.data$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network2[which(network2$shareid == ids[i]), 34])
  friends2 = which(ids %in% network2[which(network2$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 # make the network undirected
  Adj[friends2, i ] = 1
}


moran.permute.SDNN = make.permute.moran(Adj, log(all.data$SDNN), 500)
moran.permute.pNN50 = make.permute.moran(Adj, log(all.data$pNN50), 500)
moran.permute.rMSSD = make.permute.moran(Adj, log(all.data$rMSSD), 500)
moran.permute.VLF = make.permute.moran(Adj, log(all.data$VLF), 500)
moran.permute.LF = make.permute.moran(Adj, log(all.data$LF), 500)
moran.permute.HF = make.permute.moran(Adj, log(all.data$HF), 500)
moran.permute.TP = make.permute.moran(Adj, log(all.data$TP), 500)
moran.permute.ratio = make.permute.moran(Adj, log(all.data$LF/  all.data$HF), 500)

moran.permute.sex = make.permute.moran(Adj, all.data$sex, 500)
moran.permute.age = make.permute.moran(Adj, all.data$age18, 500)
moran.permute.myo = make.permute.moran(Adj, all.data$myo, 500)
moran.permute.beats = make.permute.moran(Adj, all.data$beats, 500)
moran.permute.CHF = make.permute.moran(Adj, all.data$CHF, 500)
moran.permute.diuretic = make.permute.moran(Adj, all.data$diuretic, 500)

moran.all.SDNN = make.permute.moran(Adj, fit.all.SDNN$residuals, 500)
moran.all.pNN50 = make.permute.moran(Adj, fit.all.pNN50$residuals, 500)
moran.all.rMSSD = make.permute.moran(Adj, fit.all.rMSSD$residuals, 500)
moran.all.VLF = make.permute.moran(Adj, fit.all.VLF$residuals, 500)
moran.all.LF = make.permute.moran(Adj, fit.all.LF$residuals, 500)
moran.all.HF = make.permute.moran(Adj, fit.all.HF$residuals, 500)
moran.all.TP = make.permute.moran(Adj, fit.all.TP$residuals, 500)
moran.all.ratio = make.permute.moran(Adj, fit.all.ratio$residuals, 500)


## Estimated hazard ratio for all mortality per one standard deviation decrement 
martingale.result = matrix(0, nrow = 27, ncol = 3)
colnames(martingale.result) =  c("Hazard ratio",  "95% CI", "p-value")
rownames(martingale.result) = c("Unadjusted association with all-cause mortality",
                                "ln SDNN", "ln pNN50", "ln r-MSSD", "ln VLF", "ln LF", "ln HF", "ln TP", "ln LF/HF",
                                "Age- and sex-adjusted association with all-cause mortality",
                                "ln SDNN", "ln pNN50", "ln r-MSSD", "ln VLF", "ln LF", "ln HF", "ln TP", "ln LF/HF",
                                "Age, sex, and clinical risk factors adjusted for association with all-cause morality",
                                "ln SDNN", "ln pNN50", "ln r-MSSD", "ln VLF", "ln LF", "ln HF", "ln TP", "ln LF/HF")
martingale.result[2:9,1] = formatC(c(exp(-summary(fit.SDNN)$coefficients[1]*SD[1]), exp(-summary(fit.pNN50)$coefficients[1]*SD[2]),
                             exp(-summary(fit.rMSSD)$coefficients[1]*SD[3]), exp(-summary(fit.VLF)$coefficients[1]*SD[4]),
                             exp(-summary(fit.LF)$coefficients[1]*SD[5]), exp(-summary(fit.HF)$coefficients[1]*SD[6]),
                             exp(-summary(fit.TP)$coefficients[1]*SD[7]), exp(-summary(fit.ratio)$coefficients[1]*SD[8])), format = 'f', digits = 2)

martingale.result[11:18,1] = formatC(c(exp(-summary(fit.agesex.SDNN)$coefficients[1,1]*SD[1]), exp(-summary(fit.agesex.pNN50)$coefficients[1,1]*SD[2]),
                                     exp(-summary(fit.agesex.rMSSD)$coefficients[1,1]*SD[3]), exp(-summary(fit.agesex.VLF)$coefficients[1,1]*SD[4]),
                                     exp(-summary(fit.agesex.LF)$coefficients[1,1]*SD[5]), exp(-summary(fit.agesex.HF)$coefficients[1,1]*SD[6]),
                                     exp(-summary(fit.agesex.TP)$coefficients[1,1]*SD[7]), exp(-summary(fit.agesex.ratio)$coefficients[1,1]*SD[8])), format = 'f', digits = 2)

martingale.result[20:27,1] = formatC(c(exp(-summary(fit.all.SDNN)$coefficients[1,1]*SD[1]), exp(-summary(fit.all.pNN50)$coefficients[1,1]*SD[2]),
                                       exp(-summary(fit.all.rMSSD)$coefficients[1,1]*SD[3]), exp(-summary(fit.all.VLF)$coefficients[1,1]*SD[4]),
                                       exp(-summary(fit.all.LF)$coefficients[1,1]*SD[5]), exp(-summary(fit.all.HF)$coefficients[1,1]*SD[6]),
                                       exp(-summary(fit.all.TP)$coefficients[1,1]*SD[7]), exp(-summary(fit.all.ratio)$coefficients[1,1]*SD[8])), format = 'f', digits = 2)
up =  formatC(c(exp(-SD[1]*(summary(fit.SDNN)$coefficients[1] - 1.96*summary(fit.SDNN)$coefficients[3])),
         exp(-SD[2]*(summary(fit.pNN50)$coefficients[1] - 1.96*summary(fit.pNN50)$coefficients[3])),
         exp(-SD[3]*(summary(fit.rMSSD)$coefficients[1] - 1.96*summary(fit.rMSSD)$coefficients[3])),
         exp(-SD[4]*(summary(fit.VLF)$coefficients[1] - 1.96*summary(fit.VLF)$coefficients[3])),
         exp(-SD[5]*(summary(fit.LF)$coefficients[1] - 1.96*summary(fit.LF)$coefficients[3])),
         exp(-SD[6]*(summary(fit.HF)$coefficients[1] - 1.96*summary(fit.HF)$coefficients[3])),
         exp(-SD[7]*(summary(fit.TP)$coefficients[1] - 1.96*summary(fit.TP)$coefficients[3])),
         exp(-SD[8]*(summary(fit.ratio)$coefficients[1] - 1.96*summary(fit.ratio)$coefficients[3]))), format = 'f', digits = 2)
low =  formatC(c(exp(-SD[1]*(summary(fit.SDNN)$coefficients[1] + 1.96*summary(fit.SDNN)$coefficients[3])),
        exp(-SD[2]*(summary(fit.pNN50)$coefficients[1] + 1.96*summary(fit.pNN50)$coefficients[3])),
        exp(-SD[3]*(summary(fit.rMSSD)$coefficients[1] + 1.96*summary(fit.rMSSD)$coefficients[3])),
        exp(-SD[4]*(summary(fit.VLF)$coefficients[1] + 1.96*summary(fit.VLF)$coefficients[3])),
        exp(-SD[5]*(summary(fit.LF)$coefficients[1] + 1.96*summary(fit.LF)$coefficients[3])),
        exp(-SD[6]*(summary(fit.HF)$coefficients[1] + 1.96*summary(fit.HF)$coefficients[3])),
        exp(-SD[7]*(summary(fit.TP)$coefficients[1] + 1.96*summary(fit.TP)$coefficients[3])),
        exp(-SD[8]*(summary(fit.ratio)$coefficients[1] + 1.96*summary(fit.ratio)$coefficients[3]))), format = "f", digits = 2)
martingale.result[2:9,2] = c( paste("(", low[1],", ", up[1], ")", sep=""),
                              paste("(", low[2],", ", up[2], ")", sep=""),
                              paste("(", low[3],", ", up[3], ")", sep=""),
                              paste("(", low[4],", ", up[4], ")", sep=""),
                              paste("(", low[5],", ", up[5], ")", sep=""),
                              paste("(", low[6],", ", up[6], ")", sep=""),
                              paste("(", low[7],", ", up[7], ")", sep=""),
                              paste("(", low[8],", ", up[8], ")", sep=""))
up =  formatC(c(exp(-SD[1]*(summary(fit.agesex.SDNN)$coefficients[1,1] - 1.96*summary(fit.agesex.SDNN)$coefficients[1,3])),
                exp(-SD[2]*(summary(fit.agesex.pNN50)$coefficients[1,1] - 1.96*summary(fit.agesex.pNN50)$coefficients[1,3])),
                exp(-SD[3]*(summary(fit.agesex.rMSSD)$coefficients[1,1] - 1.96*summary(fit.agesex.rMSSD)$coefficients[1,3])),
                exp(-SD[4]*(summary(fit.agesex.VLF)$coefficients[1,1] - 1.96*summary(fit.agesex.VLF)$coefficients[1,3])),
                exp(-SD[5]*(summary(fit.agesex.LF)$coefficients[1,1] - 1.96*summary(fit.agesex.LF)$coefficients[1,3])),
                exp(-SD[6]*(summary(fit.agesex.HF)$coefficients[1,1] - 1.96*summary(fit.agesex.HF)$coefficients[1,3])),
                exp(-SD[7]*(summary(fit.agesex.TP)$coefficients[1,1] - 1.96*summary(fit.agesex.TP)$coefficients[1,3])),
                exp(-SD[8]*(summary(fit.agesex.ratio)$coefficients[1,1] - 1.96*summary(fit.agesex.ratio)$coefficients[1,3]))), format = 'f', digits = 2)
low =  formatC(c(exp(-SD[1]*(summary(fit.agesex.SDNN)$coefficients[1,1] + 1.96*summary(fit.agesex.SDNN)$coefficients[1,3])),
                 exp(-SD[2]*(summary(fit.agesex.pNN50)$coefficients[1,1] + 1.96*summary(fit.agesex.pNN50)$coefficients[1,3])),
                 exp(-SD[3]*(summary(fit.agesex.rMSSD)$coefficients[1,1] + 1.96*summary(fit.agesex.rMSSD)$coefficients[1,3])),
                 exp(-SD[4]*(summary(fit.agesex.VLF)$coefficients[1,1] + 1.96*summary(fit.agesex.VLF)$coefficients[1,3])),
                 exp(-SD[5]*(summary(fit.agesex.LF)$coefficients[1,1] + 1.96*summary(fit.agesex.LF)$coefficients[1,3])),
                 exp(-SD[6]*(summary(fit.agesex.HF)$coefficients[1,1] + 1.96*summary(fit.agesex.HF)$coefficients[1,3])),
                 exp(-SD[7]*(summary(fit.agesex.TP)$coefficients[1,1] + 1.96*summary(fit.agesex.TP)$coefficients[1,3])),
                 exp(-SD[8]*(summary(fit.agesex.ratio)$coefficients[1,1] + 1.96*summary(fit.agesex.ratio)$coefficients[1,3]))), format = "f", digits = 2)
martingale.result[11:18,2] = c( paste("(", low[1],", ", up[1], ")", sep=""),
                              paste("(", low[2],", ", up[2], ")", sep=""),
                              paste("(", low[3],", ", up[3], ")", sep=""),
                              paste("(", low[4],", ", up[4], ")", sep=""),
                              paste("(", low[5],", ", up[5], ")", sep=""),
                              paste("(", low[6],", ", up[6], ")", sep=""),
                              paste("(", low[7],", ", up[7], ")", sep=""),
                              paste("(", low[8],", ", up[8], ")", sep=""))
up =  formatC(c(exp(-SD[1]*(summary(fit.all.SDNN)$coefficients[1,1] - 1.96*summary(fit.all.SDNN)$coefficients[1,3])),
                exp(-SD[2]*(summary(fit.all.pNN50)$coefficients[1,1] - 1.96*summary(fit.all.pNN50)$coefficients[1,3])),
                exp(-SD[3]*(summary(fit.all.rMSSD)$coefficients[1,1] - 1.96*summary(fit.all.rMSSD)$coefficients[1,3])),
                exp(-SD[4]*(summary(fit.all.VLF)$coefficients[1,1] - 1.96*summary(fit.all.VLF)$coefficients[1,3])),
                exp(-SD[5]*(summary(fit.all.LF)$coefficients[1,1] - 1.96*summary(fit.all.LF)$coefficients[1,3])),
                exp(-SD[6]*(summary(fit.all.HF)$coefficients[1,1] - 1.96*summary(fit.all.HF)$coefficients[1,3])),
                exp(-SD[7]*(summary(fit.all.TP)$coefficients[1,1] - 1.96*summary(fit.all.TP)$coefficients[1,3])),
                exp(-SD[8]*(summary(fit.all.ratio)$coefficients[1,1] - 1.96*summary(fit.all.ratio)$coefficients[1,3]))), format = 'f', digits = 2)
low =  formatC(c(exp(-SD[1]*(summary(fit.all.SDNN)$coefficients[1,1] + 1.96*summary(fit.all.SDNN)$coefficients[1,3])),
                 exp(-SD[2]*(summary(fit.all.pNN50)$coefficients[1,1] + 1.96*summary(fit.all.pNN50)$coefficients[1,3])),
                 exp(-SD[3]*(summary(fit.all.rMSSD)$coefficients[1,1] + 1.96*summary(fit.all.rMSSD)$coefficients[1,3])),
                 exp(-SD[4]*(summary(fit.all.VLF)$coefficients[1,1] + 1.96*summary(fit.all.VLF)$coefficients[1,3])),
                 exp(-SD[5]*(summary(fit.all.LF)$coefficients[1,1] + 1.96*summary(fit.all.LF)$coefficients[1,3])),
                 exp(-SD[6]*(summary(fit.all.HF)$coefficients[1,1] + 1.96*summary(fit.all.HF)$coefficients[1,3])),
                 exp(-SD[7]*(summary(fit.all.TP)$coefficients[1,1] + 1.96*summary(fit.all.TP)$coefficients[1,3])),
                 exp(-SD[8]*(summary(fit.all.ratio)$coefficients[1,1] + 1.96*summary(fit.all.ratio)$coefficients[1,3]))), format = "f", digits = 2)
martingale.result[20:27,2] = c( paste("(", low[1],", ", up[1], ")", sep=""),
                                paste("(", low[2],", ", up[2], ")", sep=""),
                                paste("(", low[3],", ", up[3], ")", sep=""),
                                paste("(", low[4],", ", up[4], ")", sep=""),
                                paste("(", low[5],", ", up[5], ")", sep=""),
                                paste("(", low[6],", ", up[6], ")", sep=""),
                                paste("(", low[7],", ", up[7], ")", sep=""),
                                paste("(", low[8],", ", up[8], ")", sep=""))

## pvalue
martingale.result[2:9,3] = formatC(c(summary(fit.SDNN)$coefficients[5], summary(fit.pNN50)$coefficients[5],
                                     summary(fit.rMSSD)$coefficients[5], summary(fit.VLF)$coefficients[5],
                                     summary(fit.LF)$coefficients[5], summary(fit.HF)$coefficients[5],
                                     summary(fit.TP)$coefficients[5], summary(fit.ratio)$coefficients[5]), format = 'f', digits = 4)
martingale.result[11:18,3] = formatC(c(summary(fit.agesex.SDNN)$coefficients[1,5], summary(fit.agesex.pNN50)$coefficients[1,5],
                                     summary(fit.agesex.rMSSD)$coefficients[1,5], summary(fit.agesex.VLF)$coefficients[1,5],
                                     summary(fit.agesex.LF)$coefficients[1,5], summary(fit.agesex.HF)$coefficients[1,5],
                                     summary(fit.agesex.TP)$coefficients[1,5], summary(fit.agesex.ratio)$coefficients[1,5]), format = 'f', digits = 4)
martingale.result[20:27,3] = formatC(c(summary(fit.all.SDNN)$coefficients[1,5], summary(fit.all.pNN50)$coefficients[1,5],
                                       summary(fit.all.rMSSD)$coefficients[1,5], summary(fit.all.VLF)$coefficients[1,5],
                                       summary(fit.all.LF)$coefficients[1,5], summary(fit.all.HF)$coefficients[1,5],
                                       summary(fit.all.TP)$coefficients[1,5], summary(fit.all.ratio)$coefficients[1,5]), format = 'f', digits = 4)


## Moran's I and its p-values
moran.table = matrix(0, nrow = 9, ncol = 8)
unadjusted = formatC(rbind(moran.SDNN, moran.pNN50, moran.rMSSD, moran.VLF, moran.LF, moran.HF, moran.TP, moran.ratio), digits = 3, format = 'f')
adjusted.agesex = formatC(rbind(moran.agesex.SDNN, moran.agesex.pNN50, moran.agesex.rMSSD, moran.agesex.VLF, moran.agesex.LF, moran.agesex.HF, moran.agesex.TP, moran.agesex.ratio), digits = 3, format = 'f')
adjusted.all = formatC(rbind(moran.all.SDNN, moran.all.pNN50, moran.all.rMSSD, moran.all.VLF, moran.all.LF, moran.all.HF, moran.all.TP, moran.all.ratio), digits = 3, format = 'f')

moran.table[1, ] = unadjusted[,1]
moran.table[2, ] = unadjusted[,2]
moran.table[3, ] = unadjusted[,3]
moran.table[4, ] = adjusted.agesex[,1]
moran.table[5, ] = adjusted.agesex[,2]
moran.table[6, ] = adjusted.agesex[,3]
moran.table[7, ] = adjusted.all[,1]
moran.table[8, ] = adjusted.all[,2]
moran.table[9, ] = adjusted.all[,3]

print(moran.table)

print(xtable(moran.table))