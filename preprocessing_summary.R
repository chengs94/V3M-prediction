library(dplyr)
library(randomForest)
library(ROCR)
library(mice)
library(sjmisc)
library(glmnet)
library(Metrics)

#setwd("...")
old=read.csv("assess_composite_v6.csv") # this is where I store the composite outcomes that I coded; it's called "old" but the outcome is up-to-date
#lastver=read.csv("ASSESS4_Feb19_2020.csv")
newdat=read.csv("ASSESS6_Apr1_2020.csv")
old=old[,-1]

#setdiff(names(newdat),names(lastver))
#sum(as.character(newdat$subj)==as.character(old$subj))

compdat=cbind(newdat,old$time_iped,old$iped)
colnames(compdat)[(ncol(compdat)-1):ncol(compdat)]=c("time_iped","iped")

#convert to binary outcome
compdat$bin3yr[compdat$time_iped>3 & compdat$time_iped<=36 & compdat$iped=="Yes"]=1
compdat$bin3yr[compdat$time_iped>36 | compdat$time_iped==36 & compdat$iped=="No"]=0
compdat$bin3yr[compdat$time_iped>3 & compdat$time_iped<36 & compdat$iped=="No"]=NA

compdat$diffGFR=lastver$V3M_GFR-lastver$baseline_GFR
compdat$bin3yr=as.factor(compdat$bin3yr)

##### select columns ######
dat=as_tibble(compdat)
dat2=dat %>% select(# ID and outcome, 3 variables
  subj, match_id, bin3yr,
  # Demographics and clinical, 15 variables
  CHF, CLD, COPD, CVD, Hypertension, ICU, V0_AKIN_stage, V3M_GFR, age,
  diffGFR, #setting, 
  diabetes, gender, los_index, race, sepsis,
  # Post v3m lifestyle, 5 variables
  BMI_V3M, alcohol_V3M, diastolic_bp_V3M, systolic_bp_V3M, tobacco_V3M,
  # Post V3M medications, 10 variables
  drug_class_A_V3M, drug_class_A1_V3M, drug_class_A2_V3M, drug_class_D_V3M, 
  drug_class_E_V3M, drug_class_H_V3M, drug_class_K_V3M, drug_class_O_V3M, 
  drug_class_Q_V3M, drug_class_R_V3M,
  # Adult plasma biomarkers V3M, 19+9=28 variables
  CRP_V3M, FGF23_V3M, PTH_V3M, Phos_V3M, ST2_V3M, PBNP_NT_V3M, Troponin_V3M, 
  IFNg_V3M, IL2_V3M, IL4_V3M, IL6_V3M, IL8_V3M, IL10_V3M, IL13_V3M, IL12p70_V3M, 
  IL1b_V3M, TNF_RI_V3M, TNF_RII_V3M, TNFa_V3M,
  Ang1_V3M, Ang2_V3M, PIGF_V3M, Tie2_V3M, VEGF_V3M, VEGFc_V3M, VEGFd_V3M, 
  sFlt1_V3M, bFGF_V3M,
  # Adult urine biomarkers V3M, 10 variables
  uAlbumin_V3M, uCreatinine_V3M, uCystatinC_V3M, uIL18_V3M, uKIM1_V3M, uMCP1_V3M, 
  uNGAL_V3M, uOsmolality_V3M, uUMOD_V3M, uYKL40_V3M,
  # Chemistry and lipid panel biomarkers, 6 variables
  BUN_V3M, #HDLChol_V3M, LDLChol_V3M, 
  K_V3M, glucose_V3M, sodium_V3M, 
  totchol_V3M, triglyc_V3M
)

# processing of variables
dat2$CHF[dat2$CHF=="Unknown"]=NA
dat2$CLD[dat2$CLD=="Unknown"]=NA
dat2$COPD[dat2$COPD=="Unknown"]=NA
dat2$CVD[dat2$CVD=="Unknown"]=NA
dat2$Hypertension[dat2$Hypertension=="Unknown"]=NA
#dat2$uCystatinC_V3M[!is.na(dat2$uCystatinC_V3M) & dat2$uCystatinC_V3M<0]=NA

#for (j in which(names(dat2)=="drug_class_A_V3M"): which(names(dat2)=="drug_class_R_V3M"))
#  dat2[dat2[,j]=="",j]=NA
for (j in c(which(names(dat2)=="alcohol_V3M"), which(names(dat2)=="tobacco_V3M")))
  dat2[dat2[,j]=="",j]=NA
rm(j)

dat2=droplevels(dat2)

# set NA to "No" for tobacco_V3M and alcohol_V3M
dat2$tobacco_V3M[is.na(dat2$tobacco_V3M)]="No"
dat2$alcohol_V3M[is.na(dat2$alcohol_V3M)]="No"

# missing values
nrow(dat2)-nrow(na.omit(dat2)) # 281 missing
#nrow(dat2)-nrow(na.omit(dat2 %>% select(-Troponin_V3M))) # 309 missing
#sum(apply(X=dat2,MARGIN=1,FUN=function(vec) sum(is.na(vec))!=0) & !is.na(dat2$Troponin_V3M))
#sum(is.na(dat2$Troponin_V3M))
mis.ind=is.na(dat2[which(rowSums(is.na(dat2))!=0),which(colSums(is.na(dat2))!=0)])
mis.ind=cbind(mis.ind,rowSums(mis.ind))
mis.ind=mis.ind[order(mis.ind[,ncol(mis.ind)]),]
mis.ind=mis.ind[,-ncol(mis.ind)]
#mis.ind=is.na(dat2[which(rowSums(is.na(dat2))!=0),-c(1,2)])
par(mar=c(10,4,4,2))
heatmap(matrix(as.numeric(mis.ind),ncol=sum(colSums(is.na(dat2))!=0)), scale="none", 
        col=c("white","royalblue"),Colv = NA, Rowv = NA,
        labCol=names(dat2)[which(colSums(is.na(dat2))!=0)],labRow="")

# check event time of excluded patients (event by 3mon, 12 total)
dat[dat$time_iped<=3,] %>% select(ckd_incidence, time_ckd_incidence, 
                                  ckd_progression, time_ckd_progression,
                                  esrd, time_esrd, death, time_death)
# 10 dropouts by month 3; 
# 1 ESRD at month 2.9 and drop-out at month 81.1; 
# 1 ESRD at month 2.8 and death at month 13.7

# grouping of "NonWhite"
levels(dat2$race)=c(levels(dat2$race),"NonWhite")
dat2$race[dat2$race!="White"]="NonWhite"
dat2=droplevels(dat2)

# log transform (log2)
for (j in which(names(dat2)=="CRP_V3M"): which(names(dat2)=="triglyc_V3M")){
  dat2[!is.na(dat2[,j]),j]=log2(dat2[!is.na(dat2[,j]),j])
}
rm(j)

# treating event-free subjects after 30mon as non-events
compdat$bin3yr[compdat$time_iped>36 | compdat$time_iped>=30 & compdat$iped=="No"]=0
compdat$bin3yr[compdat$time_iped>3 & compdat$time_iped<30 & compdat$iped=="No"]=NA

######## missing outcome summary #########
temp=as.data.frame(newdat[is.na(dat$bin3yr),] %>% select(ckd_incidence, time_ckd_incidence, 
                                                         ckd_progression, time_ckd_progression,
                                                         esrd, time_esrd, death, time_death))
hist(dat[is.na(dat$bin3yr),]$time_esrd,xlab="time of censoring",main="Censoring time",
     breaks=seq(0,36,6),xlim=c(0,36),xaxt='n')
axis(side=1, at=seq(0,36,6), labels=seq(0,36,6))
#c(sum(dat[is.na(dat$bin3yr),]$time_esrd<=3),
#  sum(dat[is.na(dat$bin3yr),]$time_esrd>3 & dat[is.na(dat$bin3yr),]$time_esrd<=12),
#  sum(dat[is.na(dat$bin3yr),]$time_esrd>12 & dat[is.na(dat$bin3yr),]$time_esrd<=24),
#  sum(dat[is.na(dat$bin3yr),]$time_esrd>24 & dat[is.na(dat$bin3yr),]$time_esrd<=30),
#  sum(dat[is.na(dat$bin3yr),]$time_esrd>30 & dat[is.na(dat$bin3yr),]$time_esrd<=36))
#c(mean(dat[is.na(dat$bin3yr),]$time_esrd<=3),
#  mean(dat[is.na(dat$bin3yr),]$time_esrd>3 & dat[is.na(dat$bin3yr),]$time_esrd<=12),
#  mean(dat[is.na(dat$bin3yr),]$time_esrd>12 & dat[is.na(dat$bin3yr),]$time_esrd<=24),
#  mean(dat[is.na(dat$bin3yr),]$time_esrd>24 & dat[is.na(dat$bin3yr),]$time_esrd<=30),
#  mean(dat[is.na(dat$bin3yr),]$time_esrd>30 & dat[is.na(dat$bin3yr),]$time_esrd<=36))