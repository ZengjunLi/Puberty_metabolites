library(TwoSampleMR)
library(dplyr)
library(tidyr)
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO",force = TRUE)
library(MRPRESSO)
setwd("pathname/data")


#### Childhood BMI & birth weight SNPs ----

# #cbmi snp
# Downlord: http://egg-consortium.org/childhood-obesity-2019.html
# cbmi <- read.table("CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.txt.gz",header = T,na.strings=c("-","NA"))
# cbmi[,20] <- as.numeric(cbmi[,20])
# colnames(cbmi)
# cbmi2 <- cbmi %>% drop_na(EUR_P)
# cbmi3 <- subset(cbmi2, cbmi2$EUR_P < 5e-8)
# pos_rs <- read.table("pos_rs",header = TRUE)
# colnames(pos_rs)
# colnames(cbmi3)[1] <- "chr_pos"
# cbmi4 <- tidyr::unite(cbmi3, "ID", CHR, POS,sep = ":",remove=FALSE)
# cbmi4 <- cbmi4[,-(22:27)]
# cbmi4 <- cbmi4[,-(7:16)]
# cbmi_confounding <- left_join(cbmi4, pos_rs, by = 'ID')
# cbmi_confounding <- edit(cbmi_confounding)
# https://www.ncbi.nlm.nih.gov/
# 2:632592---rs12992672
# 2:645190---rs7564167
# 2:25159858---rs10865322
# 2:25191745---rs2118826
# 2:25205427---rs1982200
# 16:53814363---rs9972653
# 16:53822651---rs7185735
# 16:53837144---11075993
# cbmi_confounding$SNP[which(cbmi_confounding$ID=='2:632592')] <- 'rs12992672'
# cbmi_confounding$SNP[which(cbmi_confounding$ID=='2:645190')] <- 'rs7564167'
# cbmi_confounding$SNP[which(cbmi_confounding$ID=='2:25159858')] <- 'rs10865322'
# cbmi_confounding$SNP[which(cbmi_confounding$ID=='2:25191745')] <- 'rs2118826'
# cbmi_confounding$SNP[which(cbmi_confounding$ID=='2:25205427')] <- 'rs1982200'
# cbmi_confounding$SNP[which(cbmi_confounding$ID=='16:53814363')] <- 'rs9972653'
# cbmi_confounding$SNP[which(cbmi_confounding$ID=='16:53822651')] <- 'rs7185735'
# cbmi_confounding$SNP[which(cbmi_confounding$ID=='16:53837144')] <- 'rs11075993'
# colnames(cbmi_confounding)[12] <- "SNP"
# write.csv(cbmi_confounding,"cbmi_snp.csv", row.names = FALSE)
# 
# #birth weight snp
# Downlord: http://egg-consortium.org/birth-weight-2019.html
# bw= read.table("birthweightSNP/Fetal_Effect_European_meta_NG2019.txt.gz", header = T)
# colnames(bw)
# bw2 <- subset(bw, bw$p < 5e-8)
# colnames(bw2)[2] <- "SNP"
# write.csv(bw2,"bw_snp.csv", row.names = FALSE)
cbmi_snp <- read.csv("cbmi_snp.csv")
cbmi_snp2 <- as.data.frame(cbmi_snp[,12])
colnames(cbmi_snp2) <- "SNP"
cbmi_snp2$confounder_flag <- "cbmi"
bw2 <- read.csv("bw_snp.csv")
bw3 <- as.data.frame(bw2[,2])
colnames(bw3) <- "SNP"
bw3$confounder_flag <- "bw"

cbmi_bw_snp <- rbind(cbmi_snp2,bw3)
cbmi_bw_snp <- cbmi_bw_snp %>% distinct(SNP, .keep_all=TRUE)
write.csv(cbmi_bw_snp,"exposure_and_confounders/cbmi_bw.csv", row.names = FALSE)
#### END: Childhood BMI & birth weight 2samplemr SNPs ----

#### Puberty timing for two-sample-mr----
#puberty clumping
# Downlord: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5841952/
day389 <- read.csv("day389.csv",header = T)
colnames(day389)[10] <- "pval.exposure"

day389_clumped <- clump_data(day389)
day389_clumped[,4] <- sub(",", "", day389_clumped[,4], fixed = TRUE)
day389_clumped[,4] <- sub(",", "", day389_clumped[,4], fixed = TRUE)
day389_clumped[,4] <- as.numeric(day389_clumped[,4])

cbmi_bw_snp <- read.csv("cbmi_bw.csv")

day208 <- left_join(day389_clumped,cbmi_bw_snp,by="SNP")
table(day208$confounder_flag)
day201 <- day208[is.na(day208$confounder_flag),]

write.csv(day208,"exposure_and_confounders/day208.csv",row.names = FALSE)
write.csv(day201,"exposure_and_confounders/day201.csv",row.names = FALSE)
#### END: Puberty timing for two-sample-mr----

#### Metabolites for two-sample-mr----
#pt(snp)_metabolities
# Downlord: https://omicscience.org/apps/crossplatform/      (according to the pt(snp))
data=data.frame()

dir("./metabolitesSNP/metabolites_PT/")

for (i in dir("./metabolitesSNP/metabolites_PT/")){ 
  df=read.csv	(paste("./metabolitesSNP/metabolites_PT/",i,sep=""))
  data=rbind(data,df)
}
#pathway
pathway <- df[,-(3:41)]
colnames(pathway)[1] <- "outcome"
colnames(pathway)[2] <- "pathway"
write.csv(pathway,"pathway.csv",row.names = FALSE)
#metabolities
data <- data[,-(16:41)]
data <- data[,-(6:9)]
colnames(data)[5] <- "SNP"
colnames(data)[8] <- "EAF"
colnames(data)[9] <- "N"
data[,9] <- sub(",", "", unlist(data[,9]), fixed = TRUE)
colnames(data)[10] <- "Z_score"
colnames(data)[11] <- "P_value"
data[,c(8:11)] <- as.numeric(unlist(data[,c(8:11)]))
data$SNP[which(data$pos=='41,494,364')] <- 'rs4924538' 
data$SNP[which(data$pos=='87,319,950')] <- 'rs7465046' 
data$SNP[which(data$pos=='91,196,099')] <- 'rs7517629' 
for (j in 1:nrow(data)){
  if(data$EAF[j]>0.5){
    data$MAF[j]=1-data$EAF[j]
  }else {
    data$MAF[j]=data$EAF[j]
  }
}

# calculate BETA and SE
data$SE <- 1 / sqrt(2 * data$MAF * (1 - data$MAF) * (data$N + data$Z_score^2))
data$BETA <- data$Z_score / sqrt(2 * data$MAF * (1 - data$MAF) * (data$N + data$Z_score^2))

# export all outcome to one file
write.csv(data,"outcome/all_pt_metabolites.csv",row.names = FALSE)
#### END: Metabolites for two-sample-mr----




#### Adulthood BMI  for two-step-MR  ----
abmi_exposure_dat <- extract_instruments(c('ebi-a-GCST006368'))
abmi_exposure_dat_clumped <- clump_data(abmi_exposure_dat)
cbmi_bw_snp <- read.csv("cbmi_bw.csv")
abmi <- left_join(abmi_exposure_dat_clumped,cbmi_bw_snp,by="SNP")
table(abmi$confounder_flag)
abmi141_sen <- abmi[is.na(abmi$confounder_flag),]
write.csv(abmi,"exposure_and_confounders/abmi147.csv",row.names = FALSE)
write.csv(abmi141_sen,"exposure_and_confounders/abmi141_sen.csv",row.names = FALSE)
#### END: Adulthood BMI for two-step-MR  ----


#### Metabolites for two-step-MR ----
#Adulthood BMI(snp)_metabolities
# Downlord: https://omicscience.org/apps/crossplatform/      (according to the ABMI(snp))
abmi_met=data.frame()

dir("./metabolitesSNP/metabolites_ABMI/")

for (i in dir("./metabolitesSNP/metabolites_ABMI/")){ 
  df=read.csv	(paste("./metabolitesSNP/metabolites_ABMI/",i,sep=""))
  abmi_met=rbind(abmi_met,df)
}

#Metabolities
abmi_met <- abmi_met[,-(16:41)]
abmi_met <- abmi_met[,-(6:9)]
colnames(abmi_met)[5] <- "SNP"
colnames(abmi_met)[8] <- "EAF"
colnames(abmi_met)[9] <- "N"
abmi_met[,9] <- sub(",", "", unlist(abmi_met[,9]), fixed = TRUE)
colnames(abmi_met)[10] <- "Z_score"
colnames(abmi_met)[11] <- "P_value"
abmi_met[,c(8:11)] <- as.numeric(unlist(abmi_met[,c(8:11)]))
for (j in 1:nrow(abmi_met)){
  if(abmi_met$EAF[j]>0.5){
    abmi_met$MAF[j]=1-abmi_met$EAF[j]
  }else {
    abmi_met$MAF[j]=abmi_met$EAF[j]
  }
}

# calculate BETA and SE
abmi_met$SE <- 1 / sqrt(2 * abmi_met$MAF * (1 - abmi_met$MAF) * (abmi_met$N + abmi_met$Z_score^2))
abmi_met$BETA <- abmi_met$Z_score / sqrt(2 * abmi_met$MAF * (1 - abmi_met$MAF) * (abmi_met$N + abmi_met$Z_score^2))

# export all outcome to one file
write.csv(abmi_met,"outcome/all_abmi_metabolites.csv",row.names = FALSE)
#### END: Metabolites for two-step-MR ----



#### Reprogen GWAS ----
# Downlord: https://reprogen.org/
reprogen <- read.table("Menarche_1KG_NatGen2017_WebsiteUpload.txt",header = T)
reprogen2 <- dplyr::filter(reprogen, !grepl(':INDEL', Markername))
reprogen2[,1] <- sub("chr", "", unlist(reprogen2[,1]), fixed = TRUE)
write.table(reprogen2,file="reprogen2.txt",row.names = FALSE, quote=F)
write.table(reprogen2,file="D:/research/linux/reprogen2.txt",row.names = FALSE, quote=F)
## Handle in Linux to match the pos and rsnumber
#cd /mnt/d/research/linux/
#cat reprogen2.txt|head -n 15
#sed -i '1d' reprogen2.txt
#cat snp150_hg19.txt|head -n 15
#sed -i '1d' snp150_hg19.txt
#sort snp150_hg19.txt  -k 1 > file2.txt && cat file2.txt|head -n 15
#join -1 1 -2 1 -a1 -eNa -o "1.1 2.2 1.2 1.3 1.4 1.5 1.6" reprogen2.txt file2.txt > file3.txt
#sed -i '1 i pos SNP Allele1 Allele2 Effect Pvalue Minor_Allele' file3.txt
reprogen3 <-read.table("D:/research/linux/file3.txt",header=TRUE)
# calculate se
reprogen3$se <-abs(reprogen3$Effect
                   /qnorm(reprogen3$Pvalue/2,lower.tail=F))
# find the new exposure
reprogen_5e <- subset(reprogen3, reprogen3$Pvalue < 5e-8)
reprogen_Na <- subset(reprogen_5e, SNP=="Na")
# 203 SNPs were lack of Rsnumber
reprogen_5e$SNP[which(reprogen_5e$pos=='14:93880734')] <- 'rs76336754'
reprogen_5e$SNP[which(reprogen_5e$pos=='2:209630397')] <- 'rs76150399'
# The rest of 201 SNPs without Rsnumber were from chr23 
reprogen_All <- dplyr::filter(reprogen_5e, !grepl('Na', SNP))
colnames(reprogen_All)[6] <- "pval.exposure"
reprogen_clumped <- clump_data(reprogen_All)
write.table(reprogen3,"reprogen_full.txt",row.names = FALSE)
write.csv(reprogen_clumped,"exposure_and_confounders/reprogen_209.csv",row.names = FALSE)
#### END: Reprogen GWAS ----

#### Metabolities for MVMR----
# # Using rs11227247 as a proxy to replace rs10750766 r2=0.35 delete
pt_mvmr_met=data.frame()
dir("./metabolitesSNP/metabolites_MVMR_reprogen/")

for (i in dir("./metabolitesSNP/metabolites_MVMR_reprogen/")){ 
  df=read.csv	(paste("./metabolitesSNP/metabolites_MVMR_reprogen/",i,sep=""))
  pt_mvmr_met=rbind(pt_mvmr_met,df)
}
pt_mvmr_met$rsid[which(pt_mvmr_met$pos=='130,774,674')] <- 'rs3111740'
pt_mvmr_met <- pt_mvmr_met[,-(16:41)]
pt_mvmr_met <- pt_mvmr_met[,-(6:9)]
colnames(pt_mvmr_met)[5] <- "SNP"
colnames(pt_mvmr_met)[8] <- "EAF"
colnames(pt_mvmr_met)[9] <- "N"
pt_mvmr_met[,9] <- sub(",", "", unlist(pt_mvmr_met[,9]), fixed = TRUE)
colnames(pt_mvmr_met)[10] <- "Z_score"
colnames(pt_mvmr_met)[11] <- "P_value"
pt_mvmr_met[,c(8:11)] <- as.numeric(unlist(pt_mvmr_met[,c(8:11)]))
for (j in 1:nrow(pt_mvmr_met)){
  if(pt_mvmr_met$EAF[j]>0.5){
    pt_mvmr_met$MAF[j]=1-pt_mvmr_met$EAF[j]
  }else {
    pt_mvmr_met$MAF[j]=pt_mvmr_met$EAF[j]
  }
}

# calculate BETA and SE
pt_mvmr_met$SE <- 1 / sqrt(2 * pt_mvmr_met$MAF * (1 - pt_mvmr_met$MAF) * (pt_mvmr_met$N + pt_mvmr_met$Z_score^2))
pt_mvmr_met$BETA <- pt_mvmr_met$Z_score / sqrt(2 * pt_mvmr_met$MAF * (1 - pt_mvmr_met$MAF) * (pt_mvmr_met$N + pt_mvmr_met$Z_score^2))

# export all outcome to one file
all_abmi_met <- read.csv("outcome/all_abmi_metabolites.csv")
all_mvmr_met <- rbind(pt_mvmr_met,all_abmi_met)
write.csv(pt_mvmr_met,"outcome/pt_mvmr_metabolites.csv",row.names = FALSE)
write.csv(all_mvmr_met,"outcome/all_mvmr_metabolites.csv",row.names = FALSE)
# seperate  
met_mvmr <- read.csv("outcome/all_mvmr_metabolites.csv")
met_mvmr[,1] <- sub(":", "_", met_mvmr[,1], fixed = TRUE)
full <- unique(met_mvmr$met_full)
for (i in full) {
  file_name = paste(i,'.csv',sep = '')
  file_path = paste('D:/research/outcome/mvmr_met_sep/',file_name,sep = '')
  file_full <- met_mvmr %>% filter(met_full == i)
  write.csv(file_full,file_path)
}

#### END: Metabolities  for MVMR----

#### Extract MVMR exposure ----
## determine eaf
reprogen209 <- read.csv("exposure_and_confounders/reprogen_209.csv")
reprogen209_in_abmigwas <- extract_outcome_data(reprogen209$SNP, c('ebi-a-GCST006368'))
reprogen209_2 <- left_join(reprogen209,reprogen209_in_abmigwas,by="SNP")
# Using rs4717904 as a proxy to replace rs2267812 r2=0.76 delete
reprogen209_2$Allele1[which(reprogen209_2$Allele1 =='a')] <- 'A'
reprogen209_2$Allele1[which(reprogen209_2$Allele1 =='g')] <- 'G'
reprogen209_2$Allele1[which(reprogen209_2$Allele1 =='c')] <- 'C'
reprogen209_2$Allele1[which(reprogen209_2$Allele1 =='t')] <- 'T'
reprogen209_2$Allele2[which(reprogen209_2$Allele2 =='a')] <- 'A'
reprogen209_2$Allele2[which(reprogen209_2$Allele2 =='g')] <- 'G'
reprogen209_2$Allele2[which(reprogen209_2$Allele2 =='c')] <- 'C'
reprogen209_2$Allele2[which(reprogen209_2$Allele2 =='t')] <- 'T'
reprogen209_2$Minor_Allele[which(reprogen209_2$Minor_Allele =='a')] <- 'A'
reprogen209_2$Minor_Allele[which(reprogen209_2$Minor_Allele =='g')] <- 'G'
reprogen209_2$Minor_Allele[which(reprogen209_2$Minor_Allele =='c')] <- 'C'
reprogen209_2$Minor_Allele[which(reprogen209_2$Minor_Allele =='t')] <- 'T'
# ncbi 	rs2267812: C=0.21742	A=0.78258
reprogen209_2$eaf.outcome[which(reprogen209_2$SNP=='rs2744688')] <- '0.78258'
del <- which(reprogen209_2$SNP=="rs10750766")
repeogen208 <- reprogen209_2[-del,]
reprogen207 <- na.omit(repeogen208)
reprogen207[,c(16)] <- as.numeric(unlist(reprogen207[,c(16)]))
for (j in 1:nrow(reprogen207)){
  if(reprogen207$eaf.outcome[j]>0.5){
    reprogen207$MAF[j]=1-reprogen207$eaf.outcome[j]
  }else reprogen207$MAF[j]=reprogen207$eaf.outcome[j]
}
reprogen207$nMAF <- 1-reprogen207$MAF
reprogen207$iden <- reprogen207$Allele1 == reprogen207$Minor_Allele
for (j in 1:nrow(reprogen207)){
  if(reprogen207$iden[j]== "TRUE" ){
    reprogen207$eaf.pt[j]=reprogen207$MAF[j]
  }else {
    reprogen207$eaf.pt[j]=reprogen207$nMAF[j]
  }
}
# reprogen_pt
reprogen_pt <- reprogen207[,-(9:27)]
reprogen_pt <- reprogen_pt[,-7]
reprogen_pt <- reprogen_pt[,-1]
colnames(reprogen_pt)[4] <- "beta"
colnames(reprogen_pt)[5] <- "p"
colnames(reprogen_pt)[7] <- "eaf"
reprogen_pt$exposure <- "pt"
reprogen_pt$N <- "179117"
write.csv(reprogen_pt,"exposure_and_confounders/reprogen_pt.csv",row.names = FALSE)
# pt_in_abmigwas
pt_in_abmigwas <- reprogen207[,-(20:28)]
pt_in_abmigwas <- pt_in_abmigwas[,-(3:11)]
pt_in_abmigwas <- pt_in_abmigwas[,-1]
colnames(pt_in_abmigwas)[2] <- "beta"
colnames(pt_in_abmigwas)[3] <- "se"
colnames(pt_in_abmigwas)[4] <- "N"
colnames(pt_in_abmigwas)[5] <- "p"
colnames(pt_in_abmigwas)[6] <- "eaf"
colnames(pt_in_abmigwas)[7] <- "Allele1"
colnames(pt_in_abmigwas)[8] <- "Allele2"
colnames(pt_in_abmigwas)[9] <- "exposure"
write.csv(pt_in_abmigwas,"exposure_and_confounders/pt_in_abmigwas.csv",row.names = FALSE)
# abmi_in_abmigwas
abmi_in_abmigwas <- read.csv("exposure_and_confounders/abmi147.csv")
abmi_in_abmigwas <- abmi_in_abmigwas[,-(13:16)]
abmi_in_abmigwas <- abmi_in_abmigwas[,-7]
abmi_in_abmigwas <- abmi_in_abmigwas[,-5]
abmi_in_abmigwas <- abmi_in_abmigwas[,-2]
colnames(abmi_in_abmigwas)[1] <- "N"
colnames(abmi_in_abmigwas)[3] <- "se"
colnames(abmi_in_abmigwas)[2] <- "beta"
colnames(abmi_in_abmigwas)[4] <- "p"
colnames(abmi_in_abmigwas)[8] <- "eaf"
colnames(abmi_in_abmigwas)[6] <- "Allele1"
colnames(abmi_in_abmigwas)[7] <- "Allele2"
colnames(abmi_in_abmigwas)[9] <- "exposure"
write.csv(abmi_in_abmigwas,"exposure_and_confounders/abmi_in_abmigwas.csv",row.names = FALSE)
# 
repeogen_full <- read.table("reprogen_full.txt",header=TRUE)
abmi <- read.csv("exposure_and_confounders/abmi_in_abmigwas.csv")
reprogen_abmi <- left_join(abmi,repeogen_full,by="SNP")
reprogen_abmi$Allele1.y[which(reprogen_abmi$Allele1.y =='a')] <- 'A'
reprogen_abmi$Allele1.y[which(reprogen_abmi$Allele1.y =='g')] <- 'G'
reprogen_abmi$Allele1.y[which(reprogen_abmi$Allele1.y =='c')] <- 'C'
reprogen_abmi$Allele1.y[which(reprogen_abmi$Allele1.y =='t')] <- 'T'
reprogen_abmi$Allele2.y[which(reprogen_abmi$Allele2.y =='a')] <- 'A'
reprogen_abmi$Allele2.y[which(reprogen_abmi$Allele2.y =='g')] <- 'G'
reprogen_abmi$Allele2.y[which(reprogen_abmi$Allele2.y =='c')] <- 'C'
reprogen_abmi$Allele2.y[which(reprogen_abmi$Allele2.y =='t')] <- 'T'
reprogen_abmi$Minor_Allele[which(reprogen_abmi$Minor_Allele =='a')] <- 'A'
reprogen_abmi$Minor_Allele[which(reprogen_abmi$Minor_Allele =='g')] <- 'G'
reprogen_abmi$Minor_Allele[which(reprogen_abmi$Minor_Allele =='c')] <- 'C'
reprogen_abmi$Minor_Allele[which(reprogen_abmi$Minor_Allele =='t')] <- 'T'
# determine eaf
reprogen_abmi[,c(8)] <- as.numeric(unlist(reprogen_abmi[,c(8)]))
for (j in 1:nrow(reprogen_abmi)){
  if(reprogen_abmi$eaf[j]>0.5){
    reprogen_abmi$MAF[j]=1-reprogen_abmi$eaf[j]
  }else reprogen_abmi$MAF[j]=reprogen_abmi$eaf[j]
}
reprogen_abmi$nMAF <- 1-reprogen_abmi$MAF
reprogen_abmi$iden <- reprogen_abmi$Allele1.y == reprogen_abmi$Minor_Allele
for (j in 1:nrow(reprogen_abmi)){
  if(reprogen_abmi$iden[j]== "TRUE" ){
    reprogen_abmi$eaf.2[j]=reprogen_abmi$MAF[j]
  }else {
    reprogen_abmi$eaf.2[j]=reprogen_abmi$nMAF[j]
  }
}
# reprogen_abmi
reprogen_abmi <- reprogen_abmi[,-(17:19)]
reprogen_abmi <- reprogen_abmi[,-15]
reprogen_abmi <- reprogen_abmi[,-(9:10)]
reprogen_abmi <- reprogen_abmi[,-(6:8)]
reprogen_abmi <- reprogen_abmi[,-(1:4)]
reprogen_abmi$exposure <- "pt"
reprogen_abmi$N <- "179117"
colnames(reprogen_abmi)[2] <- "Allele1"
colnames(reprogen_abmi)[3] <- "Allele2"
colnames(reprogen_abmi)[4] <- "beta"
colnames(reprogen_abmi)[5] <- "p"
colnames(reprogen_abmi)[6] <- "se"
colnames(reprogen_abmi)[7] <- "eaf"
write.csv(reprogen_abmi,"exposure_and_confounders/reprogen_abmi.csv",row.names = FALSE)
# pt_mvmr
reprogen_pt <- read.csv("exposure_and_confounders/reprogen_pt.csv")
reprogen_abmi <- read.csv("exposure_and_confounders/reprogen_abmi.csv")
all_pt_mvmr <- rbind(reprogen_pt,reprogen_abmi)
all_pt_mvmr <- all_pt_mvmr %>% distinct(SNP, .keep_all=TRUE)

cbmi_bw_snp <- read.csv("cbmi_bw.csv")
pt_MVMR_sen <- left_join(all_pt_mvmr,cbmi_bw_snp,by="SNP")
table(pt_MVMR_sen$confounder_flag)
pt_MVMR_sen2 <- pt_MVMR_sen[is.na(pt_MVMR_sen$confounder_flag),]

write.csv(all_pt_mvmr,"exposure_and_confounders/all_pt_mvmr.csv", row.names = FALSE)
write.csv(pt_MVMR_sen2,"exposure_and_confounders/sen_pt_MVMR.csv", row.names = FALSE)
# abmi_mvmr
pt_in_abmigwas <- read.csv("exposure_and_confounders/pt_in_abmigwas.csv")
abmi_in_abmigwas <- read.csv("exposure_and_confounders/abmi_in_abmigwas.csv")
colnames(pt_in_abmigwas)[9] <- "exposure"
all_abmi_mvmr <- rbind(pt_in_abmigwas,abmi_in_abmigwas)
all_abmi_mvmr <- all_abmi_mvmr %>% distinct(SNP, .keep_all=TRUE)

cbmi_bw_snp <- read.csv("cbmi_bw.csv")
abmi_MVMR_sen <- left_join(all_abmi_mvmr,cbmi_bw_snp,by="SNP")
table(abmi_MVMR_sen$confounder_flag)
abmi_MVMR_sen2 <- abmi_MVMR_sen[is.na(abmi_MVMR_sen$confounder_flag),]
write.csv(all_abmi_mvmr,"exposure_and_confounders/all_abmi_mvmr.csv", row.names = FALSE)
write.csv(abmi_MVMR_sen2,"exposure_and_confounders/sen_abmi_MVMR.csv", row.names = FALSE)
#### END: Extractc MVMR exposure ----