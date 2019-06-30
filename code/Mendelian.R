library(AER)
library(netdep)
### read CARe data ###
## https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000282.v19.p11
#CareSampleInfo = read.table("phg000076.v6.FHS_CARe.sample-info.MULTI/phg000076.v6_release_manifest.txt", sep = "\t", header = TRUE)

## set 1
CareSampleInfo.c1 = CareSampleInfo[CareSampleInfo$Tar_Name == "phg000076.v6.FHS_CARe.genotype-calls-indfmt.c1.HMB-IRB-MDS.set1.tar", ]
Care.set1.c1.rs9939609 = c() 
for(i in 1:nrow(CareSampleInfo.c1)){
  filename = paste0("/phg000076.v6.FHS_CARe.genotype-calls-indfmt.c1.set1/",
                    as.character(CareSampleInfo.c1$Sample_ID[i]), ".DBGAP.TXT")
  tmp = read.table(filename, sep = " ", header = FALSE)
  colnames(tmp) = c("SNP_id", "allele1", "allele2", "GenCall_score", "norm_intensity1", "norm_intensity2")
  Care.set1.c1.rs9939609 =  rbind(Care.set1.c1.rs9939609, c(as.character(CareSampleInfo.c1$Subject_ID[i]), 
                                                            paste0(tmp$allele1[which(tmp$SNP_id == "rs9939609")], tmp$allele2[which(tmp$SNP_id == "rs9939609")])))
}
colnames(Care.set1.c1.rs9939609) = c("Subject_ID", "rs9939609")
Care.set1.c1.rs9939609 = as.data.frame(Care.set1.c1.rs9939609)

## set 2
CareSampleInfo.c1 = CareSampleInfo[CareSampleInfo$Tar_Name == "phg000076.v6.FHS_CARe.genotype-calls-indfmt.c1.HMB-IRB-MDS.set2.tar", ]
Care.set2.c1.rs9939609 = c() 
for(i in 1:nrow(CareSampleInfo.c1)){
  filename = paste0("/phg000076.v6.FHS_CARe.genotype-calls-indfmt.c1.set2/",
                    as.character(CareSampleInfo.c1$Sample_ID[i]), ".DBGAP.TXT")
  tmp = read.table(filename, sep = " ", header = FALSE)
  colnames(tmp) = c("SNP_id", "allele1", "allele2", "GenCall_score", "norm_intensity1", "norm_intensity2")
  Care.set2.c1.rs9939609 =  rbind(Care.set2.c1.rs9939609, c(as.character(CareSampleInfo.c1$Subject_ID[i]), 
                                                            paste0(tmp$allele1[which(tmp$SNP_id == "rs9939609")], tmp$allele2[which(tmp$SNP_id == "rs9939609")])))
}

colnames(Care.set2.c1.rs9939609) = c("Subject_ID", "rs9939609")
Care.set2.c1.rs9939609 = as.data.frame(Care.set2.c1.rs9939609)

Care.c1.rs9939609 = rbind(Care.set1.c1.rs9939609, Care.set2.c1.rs9939609 )
save(Care.c1.rs9939609, file = "Data/Care.set1.c1.rs9939609.RData")

## read network data
#network_c1 = read.table("/phs000153.v9.pht000836.v8.p8.c1.vr_sntwk_2008_m_0641s.DS-SN-IRB-MDS.txt", sep = "\t", header = TRUE)
spouse_network =  network_c1[network_c1$ALTERTYPE == "SPOUSE",]
spouse_network  = spouse_network[spouse_network$idtype == 1 & spouse_network$alteridtype == 1,]
spouse_network = spouse_network[spouse_network$SPELLBEGIN == 1,]

## birth era
#agedat = read.table("phs000007.v29.pht003099.v4.p10.c2.vr_dates_2014_a_0912s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
agedat$birthera = ifelse(agedat$age1 >= 31, 1, ifelse(agedat$age1 <= 25, 2, 3))
agedat$yrbirth = 1973 - agedat$age1 

## download phenotype data for each cohort
#dat1 = read.table("/phs000007.v30.pht000030.v8.p11.c1.ex1_1s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#dat2 = read.table("/phs000007.v30.pht000031.v8.p11.c1.ex1_2s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#dat3 = read.table("/phs000007.v30.pht000032.v7.p11.c1.ex1_3s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#dat4 = read.table("/phs000007.v30.pht000033.v9.p11.c1.ex1_4s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#dat5 = read.table("/phs000007.v30.pht000034.v8.p11.c1.ex1_5s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#dat6 = read.table("/phs000007.v30.pht000035.v9.p11.c1.ex1_6s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)
#dat7 = read.table("/phs000007.v30.pht000036.v9.p11.c1.ex1_7s.HMB-IRB-MDS.txt", sep = "\t", header = TRUE)

## extract bmi and smoke (time-dependent) from each of phenotype data
## dat1
dat1$bmi = (dat1$A50/2.205)/(dat1$A51/39.37)^2
dat1$smoke = (dat1$A96 == 1 & dat1$A101 == 0)
## dat2
dat2$bmi = dat2$B13 / (dat2$B16/39.37)^2
dat2$smoke = (dat2$B87 == 1)
## dat3
dat3$bmi = (dat3$C416/2.205)/(dat3$C417/39.37)^2
dat3$smoke = (dat3$C67== 1)
## dat4
dat4$bmi= (dat4$D401/2.205) / (dat4$D402/39.37)^2
dat4$smoke = (dat4$D093 == 1)
## dat5
dat5$bmi = (dat5$E024/2.205) / (dat5$E025/39.37)^2
dat5$smoke = (dat5$E319 == 1)
## dat6
dat6$bmi = (dat6$F007/2.205) / (dat6$F008/39.37)^2
dat6$smoke = (dat6$F288 == 1)
## dat7
dat7$bmi = (dat7$G440/2.205) / (dat7$G441/39.37)^2
dat7$smoke = (dat7$G116 == 1)

## construct a dyad data of spouse 
mat = c()
mat = cbind(rep(i,7), spouse_network$shareid[1], spouse_network$sharealterid[1] , c(1:7))
for(i in 2:nrow(spouse_network)){
  id1 = spouse_network$shareid[i]
  id2 = spouse_network$sharealterid[i] 
  mat = rbind(mat, cbind(rep(i,7), id1, id2, c(1:7)))
}
colnames(mat) = c("dyad", "eid", "aid", "wave")
mat = as.data.frame(mat)

## construct longitudinal data using the phenotype + genotype data per one dyad
matdist = esibs = asibs = eage = aage = eeduc = aeduc = emale = amale = ebmi = abmi = esmoke = asmoke = 
  eera2 = aera2 = eera3 = aera3 = eyrbirth = ayrbirth = 
  egene1 = egene2 = egene3 = egene4 = agene1 = agene2 = agene3 = agene4 =  rep(NA, nrow(mat))
count = 0
for(i in 1:nrow(mat)){
  dummy = (spouse_network$shareid == mat$eid[(count+1)] & spouse_network$sharealterid == mat$aid[(count+1)])
  matdist[(count+1):(count+7)] = c(spouse_network$DISTMI1[dummy],spouse_network$DISTMI2[dummy],
                                   spouse_network$DISTMI3[dummy],spouse_network$DISTMI4[dummy],
                                   spouse_network$DISTMI5[dummy],spouse_network$DISTMI6[dummy],
                                   spouse_network$DISTMI7[dummy])
  esibs[(count+1):(count+7)] = sum(network_c2$ALTERTYPE[network_c2$shareid %in% mat$eid[(count+1)]] == "BROTHER" |
                                     network_c2$ALTERTYPE[network_c2$shareid %in% mat$eid[(count+1)]] == "SISTER")
  asibs[(count+1):(count+7)] = sum(network_c2$ALTERTYPE[network_c2$shareid %in% mat$aid[(count+1)]] == "BROTHER" |
                                     network_c2$ALTERTYPE[network_c2$shareid %in% mat$aid[(count+1)]] == "SISTER")
  
  emale[(count+1):(count+7)] = rep(ifelse(sum(agedat$shareid %in% mat$eid[count+1]) == 0, NA,
                                          agedat$sex[agedat$shareid %in% mat$eid[count+1]] == 1), 7)
  amale[(count+1):(count+7)] =  rep(ifelse(sum(agedat$shareid %in% mat$aid[count+1]) == 0, NA,
                                           agedat$sex[agedat$shareid %in% mat$aid[count+1]] == 1), 7)
  if(sum(agedat$shareid %in% mat$eid[count+1]) >0){
    eage[(count+1):(count+7)] = c(agedat$age1[agedat$shareid %in% mat$eid[count+1]],
                                  agedat$age2[agedat$shareid %in% mat$eid[count+1]],
                                  agedat$age3[agedat$shareid %in% mat$eid[count+1]],
                                  agedat$age4[agedat$shareid %in% mat$eid[count+1]],
                                  agedat$age5[agedat$shareid %in% mat$eid[count+1]],
                                  agedat$age6[agedat$shareid %in% mat$eid[count+1]],
                                  agedat$age7[agedat$shareid %in% mat$eid[count+1]])
    eera2[(count+1):(count+7)] = rep(agedat$birthera[agedat$shareid %in% mat$eid[count+1]] == 2,7)
    eera3[(count+1):(count+7)] = rep(agedat$birthera[agedat$shareid %in% mat$eid[count+1]] == 3,7)
    eyrbirth[(count+1):(count+7)] = rep(agedat$yrbirth[agedat$shareid %in% mat$eid[count+1]],7)
  }
  
  if(sum(agedat$shareid %in% mat$aid[count+1]) >0){
    aage[(count+1):(count+7)] = c(agedat$age1[agedat$shareid %in% mat$aid[count+1]],
                                  agedat$age2[agedat$shareid %in% mat$aid[count+1]],
                                  agedat$age3[agedat$shareid %in% mat$aid[count+1]],
                                  agedat$age4[agedat$shareid %in% mat$aid[count+1]],
                                  agedat$age5[agedat$shareid %in% mat$aid[count+1]],
                                  agedat$age6[agedat$shareid %in% mat$aid[count+1]],
                                  agedat$age7[agedat$shareid %in% mat$aid[count+1]])
    aera2[(count+1):(count+7)] = rep(agedat$birthera[agedat$shareid %in% mat$aid[count+1]] == 2,7)
    aera3[(count+1):(count+7)] = rep(agedat$birthera[agedat$shareid %in% mat$aid[count+1]] == 3,7)
    ayrbirth[(count+1):(count+7)] = rep(agedat$yrbirth[agedat$shareid %in% mat$aid[count+1]],7)
    
  }
  eeduc[(count+1):(count+7)] = rep(ifelse(length(dat2$B43[dat2$shareid %in% mat$eid[count+1]]) == 0, NA, dat2$B43[dat2$shareid %in% mat$eid[count+1]]),7)
  aeduc[(count+1):(count+7)] = rep(ifelse(length(dat2$B43[dat2$shareid %in% mat$aid[count+1]]) == 0, NA, dat2$B43[dat2$shareid %in% mat$aid[count+1]]),7)
  
  ebmi[count+1] = ifelse(length(dat1$bmi[dat1$shareid %in% mat$eid[count+1]]) == 0, NA, dat1$bmi[dat1$shareid %in% mat$eid[count+1]])
  ebmi[count+2] = ifelse(length(dat2$bmi[dat2$shareid %in% mat$eid[count+1]]) == 0, NA, dat2$bmi[dat2$shareid %in% mat$eid[count+1]])
  ebmi[count+3] = ifelse(length(dat3$bmi[dat3$shareid %in% mat$eid[count+1]]) == 0, NA, dat3$bmi[dat3$shareid %in% mat$eid[count+1]])
  ebmi[count+4] = ifelse(length(dat4$bmi[dat4$shareid %in% mat$eid[count+1]]) == 0, NA, dat4$bmi[dat4$shareid %in% mat$eid[count+1]])
  ebmi[count+5] = ifelse(length(dat5$bmi[dat5$shareid %in% mat$eid[count+1]]) == 0, NA, dat5$bmi[dat5$shareid %in% mat$eid[count+1]])
  ebmi[count+6] = ifelse(length(dat6$bmi[dat6$shareid %in% mat$eid[count+1]]) == 0, NA, dat6$bmi[dat6$shareid %in% mat$eid[count+1]])
  ebmi[count+7] = ifelse(length(dat7$bmi[dat7$shareid %in% mat$eid[count+1]]) == 0, NA, dat7$bmi[dat7$shareid %in% mat$eid[count+1]])
  
  
  abmi[count+1] = ifelse(length(dat1$bmi[dat1$shareid %in% mat$aid[count+1]]) == 0, NA, dat1$bmi[dat1$shareid %in% mat$aid[count+1]])
  abmi[count+2] = ifelse(length(dat2$bmi[dat2$shareid %in% mat$aid[count+1]]) == 0, NA, dat2$bmi[dat2$shareid %in% mat$aid[count+1]])
  abmi[count+3] = ifelse(length(dat3$bmi[dat3$shareid %in% mat$aid[count+1]]) == 0, NA, dat3$bmi[dat3$shareid %in% mat$aid[count+1]])
  abmi[count+4] = ifelse(length(dat4$bmi[dat4$shareid %in% mat$aid[count+1]]) == 0, NA, dat4$bmi[dat4$shareid %in% mat$aid[count+1]])
  abmi[count+5] = ifelse(length(dat5$bmi[dat5$shareid %in% mat$aid[count+1]]) == 0, NA, dat5$bmi[dat5$shareid %in% mat$aid[count+1]])
  abmi[count+6] = ifelse(length(dat6$bmi[dat6$shareid %in% mat$aid[count+1]]) == 0, NA, dat6$bmi[dat6$shareid %in% mat$aid[count+1]])
  abmi[count+7] = ifelse(length(dat7$bmi[dat7$shareid %in% mat$aid[count+1]]) == 0, NA, dat7$bmi[dat7$shareid %in% mat$aid[count+1]])
  
  
  
  esmoke[count+1] = ifelse(length(dat1$smoke[dat1$shareid %in% mat$eid[count+1]]) == 0, NA, dat1$smoke[dat1$shareid %in% mat$eid[count+1]])
  esmoke[count+2] = ifelse(length(dat2$smoke[dat2$shareid %in% mat$eid[count+1]]) == 0, NA, dat2$smoke[dat2$shareid %in% mat$eid[count+1]])
  esmoke[count+3] = ifelse(length(dat3$smoke[dat3$shareid %in% mat$eid[count+1]]) == 0, NA, dat3$smoke[dat3$shareid %in% mat$eid[count+1]])
  esmoke[count+4] = ifelse(length(dat4$smoke[dat4$shareid %in% mat$eid[count+1]]) == 0, NA, dat4$smoke[dat4$shareid %in% mat$eid[count+1]])
  esmoke[count+5] = ifelse(length(dat5$smoke[dat5$shareid %in% mat$eid[count+1]]) == 0, NA, dat5$smoke[dat5$shareid %in% mat$eid[count+1]])
  esmoke[count+6] = ifelse(length(dat6$smoke[dat6$shareid %in% mat$eid[count+1]]) == 0, NA, dat6$smoke[dat6$shareid %in% mat$eid[count+1]])
  esmoke[count+7] = ifelse(length(dat7$smoke[dat7$shareid %in% mat$eid[count+1]]) == 0, NA, dat7$smoke[dat7$shareid %in% mat$eid[count+1]])
  
  asmoke[count+1] = ifelse(length(dat1$smoke[dat1$shareid %in% mat$aid[count+1]]) == 0, NA, dat1$smoke[dat1$shareid %in% mat$aid[count+1]])
  asmoke[count+2] = ifelse(length(dat2$smoke[dat2$shareid %in% mat$aid[count+1]]) == 0, NA, dat2$smoke[dat2$shareid %in% mat$aid[count+1]])
  asmoke[count+3] = ifelse(length(dat3$smoke[dat3$shareid %in% mat$aid[count+1]]) == 0, NA, dat3$smoke[dat3$shareid %in% mat$aid[count+1]])
  asmoke[count+4] = ifelse(length(dat4$smoke[dat4$shareid %in% mat$aid[count+1]]) == 0, NA, dat4$smoke[dat4$shareid %in% mat$aid[count+1]])
  asmoke[count+5] = ifelse(length(dat5$smoke[dat5$shareid %in% mat$aid[count+1]]) == 0, NA, dat5$smoke[dat5$shareid %in% mat$aid[count+1]])
  asmoke[count+6] = ifelse(length(dat6$smoke[dat6$shareid %in% mat$aid[count+1]]) == 0, NA, dat6$smoke[dat6$shareid %in% mat$aid[count+1]])
  asmoke[count+7] = ifelse(length(dat7$smoke[dat7$shareid %in% mat$aid[count+1]]) == 0, NA, dat7$smoke[dat7$shareid %in% mat$aid[count+1]])
  
  
  # FTO : (AA-reference, AT-gene1, TT-gene2)
  FTO.index = as.numeric(as.character(Care.c1.rs9939609$Subject_ID))
  if( sum(FTO.index %in% mat$eid[count+1]) > 0){
    egene1[(count+1):(count+7)] = rep(as.character(Care.c1.rs9939609$rs9939609[FTO.index %in% mat$eid[count+1]])[1] == "AT",7)
    egene2[(count+1):(count+7)] = rep(as.character(Care.c1.rs9939609$rs9939609[FTO.index %in% mat$eid[count+1]])[1] == "TT",7)
  }
  count = i*7
}

## part of the analyses below replicate the sample code available at https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12172
datc1 = data.frame(dyad = mat$dyad, eid = mat$eid, aid = mat$aid,  wave = mat$wave, ebmi = ebmi,
                   emale = emale, eage = eage, eera2 = eera2, eera3 = eera3, eyrbirth = eyrbirth,
                   esmoke = esmoke, esibs = esibs, eeduc = eeduc,
                   egene1 = egene1, egene2 = egene2, 
                   abmi = abmi, amale = amale, aage = aage, 
                   aera2 = aera2, aera3 = aera3, ayrbirth = ayrbirth, asmoke = asmoke,
                   asibs = asibs, aeduc = aeduc, 
                   agene1 = agene1, agene2 = agene2)
datc1 = as.data.frame(datc1)

# Generate male by age interactions
datc1$amaleage = datc1$amale*datc1$aage
datc1$emaleage = datc1$emale*datc1$eage

# Generate age, age^2 and gene interaction variables
datc1$agene1age = datc1$agene1*datc1$aage
datc1$agene2age = datc1$agene2*datc1$aage
datc1$aage2 = datc1$aage^2
datc1$agene1age2 = datc1$agene1age*datc1$aage
datc1$agene2age2 = datc1$agene2age*datc1$aage
datc1$egene1age = datc1$egene1*datc1$eage
datc1$egene2age = datc1$egene2*datc1$eage
datc1$eage2 = datc1$eage^2
datc1$egene1age2 = datc1$egene1age*datc1$eage
datc1$egene2age2 = datc1$egene2age*datc1$eage


datc1 = na.omit(datc1)
network = network_c1[(network_c1$SPELLBEGIN <= 12*32) & (network_c1$SPELLEND > 0),  ]

## pick (1,2)
datc1_12 = datc1[datc1$wave == 1 |datc1$wave == 2, ]
names(table(datc1_12$dyad))[which(table(datc1_12$dyad) == 2)]
datc1_12 = datc1_12[datc1_12$dyad %in% names(table(datc1_12$dyad))[which(table(datc1_12$dyad) == 2)], ]
sum(which(table(datc1_12$dyad) != 2))

eid_t =  ebmi_t = eage_t = emaleage_t = esmoke_t = esibs_t = eeduc_t = 
  egene1age_t = egene2age_t = abmi_t_1 = aage_t_1 = 
  amaleage_t_1 = asmoke_t_1 = asibs_t_1 = aeduc_t_1 = agene1age_t_1 = agene2age_t_1 = 
  eera2_t = aera2_t_1 = eera3_t = aera3_t_1 = eyrbirth_t = ayrbirth_t_1  = 
  emale_t = amale_t_1 = rep(NA, (nrow(datc1_12)/2))

for(i in 1:(nrow(datc1_12)/2)){
  eid_t[i] = datc1_12$eid[2*i]
  ebmi_t[i] = datc1_12$ebmi[2*i]
  eage_t[i] = datc1_12$eage[2*i]
  emaleage_t[i] = datc1_12$amaleage[2*i]
  esmoke_t[i] = datc1_12$esmoke[2*i]
  esibs_t[i] = datc1_12$esibs[2*i]
  eeduc_t[i] = datc1_12$eeduc[2*i]
  egene1age_t[i] = datc1_12$egene1age[2*i]
  egene2age_t[i] = datc1_12$egene2age[2*i]
  eera2_t[i] = datc1_12$eera2[2*i]
  eera3_t[i] = datc1_12$eera3[2*i]
  eyrbirth_t[i] = datc1_12$eyrbirth[2*i]
  emale_t[i] = datc1_12$emale[2*i]
  
  abmi_t_1[i] = datc1_12$abmi[2*i-1]
  aage_t_1[i] = datc1_12$aage[2*i-1]
  amaleage_t_1[i] = datc1_12$amaleage[2*i-1]
  asmoke_t_1[i] = datc1_12$asmoke[2*i-1]
  asibs_t_1[i] = datc1_12$asibs[2*i-1]
  aeduc_t_1[i] = datc1_12$aeduc[2*i-1]
  agene1age_t_1[i] = datc1_12$agene1age[2*i-1]
  agene2age_t_1[i] = datc1_12$agene2age[2*i-1]
  aera2_t_1[i] = datc1_12$aera2[2*i-1]
  aera3_t_1[i] = datc1_12$aera3[2*i-1]
  ayrbirth_t_1[i] = datc1_12$ayrbirth[2*i-1]
  amale_t_1[i] = datc1_12$amale[2*i-1]
}
tmp.dat = data.frame(eid_t = eid_t, ebmi_t = ebmi_t, eage_t = eage_t,
                     emaleage_t = emaleage_t, esmoke_t = esmoke_t,
                     esibs_t = esibs_t, eeduc_t = eeduc_t,
                     egene1age_t = egene1age_t, egene2age_t = egene2age_t, emale_t = emale_t,
                     abmi_t_1 = abmi_t_1,  aage_t_1 = aage_t_1,
                     amaleage_t_1 = amaleage_t_1, asmoke_t_1 = asmoke_t_1,
                     asibs_t_1 = asibs_t_1, aeduc_t_1 = aeduc_t_1,
                     agene1age_t_1 = agene1age_t_1, agene2age_t_1 = agene2age_t_1,
                     eera2_t = eera2_t, eera3_t = eera3_t, eyrbirth_t = eyrbirth_t,
                     aera2_t_1 = aera2_t_1, aera3_t_1 = aera3_t_1, ayrbirth_t_1 = ayrbirth_t_1,
                     amale_t_1 = amale_t_1)


result = ivreg(ebmi_t ~ abmi_t_1 + eage_t + emaleage_t + esmoke_t + esibs_t +
                 eeduc_t + egene1age_t + egene2age_t + eera2_t + eera3_t + eyrbirth_t + emale_t +
                 aage_t_1 + amaleage_t_1 + asmoke_t_1 +  
                 asibs_t_1 + aeduc_t_1 + aera2_t_1 + 
                 aera3_t_1 + ayrbirth_t_1 + amale_t_1 | .-abmi_t_1 + agene1age_t_1 + agene2age_t_1, data = tmp.dat)


Adj = matrix(0, nrow(tmp.dat), nrow(tmp.dat))
ids = tmp.dat$eid_t
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network$sharealterid[which(network$shareid == ids[i])])
  friends2 = which(ids %in% network$shareid[which(network$sharealterid == ids[i])])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i] = 1
}

result12.residual = make.permute.moran(Adj, result$residuals, np = 500)
result12.ebmi = make.permute.moran(Adj, tmp.dat$ebmi_t, np = 500)
result12.abmi = make.permute.moran(Adj, tmp.dat$abmi_t_1, np = 500)
result12.agene1age = make.permute.moran(Adj, tmp.dat$agene1age, np = 500)
result12.agene2age = make.permute.moran(Adj, tmp.dat$agene2age, np = 500)