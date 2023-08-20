library(TwoSampleMR)
library(dplyr)
library(tidyr)
setwd("pathname/data")


#### two-step-MR----
#two-step-MR Step1
#two-step-MR Step1    exposure:ABMI

twostep1_abmi <- read_exposure_data(
  filename = "exposure_and_confounders/abmi147.csv", 
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  pval_col = "pval.exposure",
  se_col = "se.exposure",
  samplesize_col = "samplesize.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure" )

#two-step-MR Step1  outcome: Metabolites
twostep1_met <- read.csv("outcome/all_abmi_metabolites.csv")
colnames(twostep1_met)
twostep1_met <- read_outcome_data(
  filename = "outcome/all_abmi_metabolites.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "EAF" ,
  pval_col = "P_value",
  phenotype_col = "met_full"
) 

# # Using rs656326 as a proxy to replace rs1561277 r2=0.72 delete

twostep1h <- harmonise_data(twostep1_abmi, twostep1_met)
table(twostep2h$outcome)
twostep1h2 <- twostep2h[twostep1h$mr_keep,]
table(twostep1h2$outcome)
pathway <- read.csv("pathway.csv")

#two-step-MR Step1 MR
twostep1_res <- mr(twostep1h, method_list=c('mr_ivw_mre'))
twostep1_res <- merge(twostep1_res,pathway,by="outcome")
twostep1_res$lower = twostep1_res$b-1.96*twostep1_res$se
twostep1_res$upper = twostep1_res$b+1.96*twostep1_res$se
write.csv(twostep1_res,"result/step1res.csv",row.names = FALSE)


#two-step-MR Step2   ABMI on Metabolites
#two-step-MR Step2   exposure: PT
twostep2_exposure <- read_exposure_data(
  filename = "exposure_and_confounders/day208.csv", 
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA..y.allele.",
  pval_col = "pval.exposure",
  se_col = "SE",
  samplesize_col = "N",
  effect_allele_col = "Effect.allele",
  other_allele_col = "Other.allele",
  eaf_col = "EAF" 
)

#two-step-MR Step2  outcome: Adulthood BMI
twostep2_outcome <- extract_outcome_data(twostep2_exposure$SNP, c('ebi-a-GCST006368'))
# Using rs4717904 as a proxy to replace rs2267812 r2=0.76 delete
#two-step-MR Step2  harmonise
step2h <- harmonise_data(twostep2_exposure, twostep2_outcome)
step2h2 <- step2h[step2h$mr_keep,]
table(step2h2$outcome)

#two-step-MR Step2 MR
step2res <- mr(step2h2, method_list=c('mr_ivw_mre',
                                      'mr_egger_regression',
                                      'mr_weighted_median',
                                      'mr_weighted_mode' ))

step2res$lower = step2res$b-1.96*step2res$se
step2res$upper = step2res$b+1.96*step2res$se


write.csv(step2res,"result/step2res.csv",row.names = FALSE)

#### END: two-step-MR ----

#### two-step-MR (Excluding 6 snps) ----
#two-step-MR Step1
#two-step-MR Step1    exposure:ABMI

twostep1_abmi_sen <- read_exposure_data(
  filename = "exposure_and_confounders/abmi141_sen.csv", 
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  pval_col = "pval.exposure",
  se_col = "se.exposure",
  samplesize_col = "samplesize.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure" )

#two-step-MR Step1  outcome: Metabolites
twostep1_met <- read.csv("outcome/all_abmi_metabolites.csv")
colnames(twostep1_met)
twostep1_met <- read_outcome_data(
  filename = "outcome/all_abmi_metabolites.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "EAF" ,
  pval_col = "P_value",
  phenotype_col = "met_full"
) 

# # Using rs656326 as a proxy to replace rs1561277 r2=0.72 delete
twostep1h_sen <- harmonise_data(twostep1_abmi_sen, twostep1_met)
table(twostep2h_sen$outcome)
twostep1h2_sen <- twostep1h_sen[twostep1h_sen$mr_keep,]
table(twostep1h2_sen$outcome)
pathway <- read.csv("pathway.csv")

#two-step-MR Step1 MR
twostep1_res_sen <- mr(twostep1h_sen, method_list=c('mr_ivw_mre'))
twostep1_res_sen <- merge(twostep1_res_sen,pathway,by="outcome")
twostep1_res_sen$lower = twostep1_res_sen$b-1.96*twostep1_res_sen$se
twostep1_res_sen$upper = twostep1_res_sen$b+1.96*twostep1_res_sen$se


write.csv(twostep1_res_sen,"result/step1res_sen.csv",row.names = FALSE)

#two-step-MR Step2   ABMI on Metabolites
#two-step-MR Step2   exposure: PT

twostep2_exposure_sensitivity <- read_exposure_data(
  filename = "exposure_and_confounders/day201.csv", 
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA..y.allele.",
  pval_col = "pval.exposure",
  se_col = "SE",
  samplesize_col = "N",
  effect_allele_col = "Effect.allele",
  other_allele_col = "Other.allele",
  eaf_col = "EAF" 
)
#two-step-MR Step2  outcome: Adulthood BMI
twostep2_outcome <- extract_outcome_data(twostep2_exposure_sensitivity$SNP, c('ebi-a-GCST006368'))
step2h_sen <- harmonise_data(twostep2_exposure_sensitivity, twostep2_outcome)
step2h2_sen <- step2h_sen[step2h_sen$mr_keep,]
table(step2h2_sen$outcome)

#two-step-MR Step2 MR
step2res_sen <- mr(step2h2_sen, method_list=c('mr_ivw_mre',
                                              'mr_egger_regression',
                                              'mr_weighted_median',
                                              'mr_weighted_mode' ))

step2res_sen$lower = step2res_sen$b-1.96*step2res_sen$se
step2res_sen$upper = step2res_sen$b+1.96*step2res_sen$se

write.csv(step2res_sen,"result/step2res_sen.csv",row.names = FALSE)


#### END: two-step-MR (Excluding 6 snps)----


#### mediated----
step1res <- read.csv("result/step1res.csv")
step1res_sen <- read.csv("result/step1res_sen.csv")
step1res <- subset(step1res,step1res$method== "Inverse variance weighted (multiplicative random effects)")
step1res_sen <- subset(step1res_sen,step1res_sen$method== "Inverse variance weighted (multiplicative random effects)") <- read.csv("result/step2res.csv")
twostep2_res <- read.csv("result/step2res.csv")
twostep2_res_sen <- read.csv("result/step2res_sen.csv")
step2res <- subset(twostep2_res,twostep2_res$method== "Inverse variance weighted (multiplicative random effects)")
step2res_sen <- subset(twostep2_res_sen,twostep2_res_sen$method== "Inverse variance weighted (multiplicative random effects)")

res_2smr <- read.csv("result/pt_met_2smr_results(ivw).csv")
res_2smr_sen <- read.csv("result/pt_met_2smr_results(ivw)_sensitivity.csv")


step1res$b2=step2res$b
step1res$se2=step2res$se
step1res$p2=step2res$p
step1res$nsnp2=step2res_sen$nsnp
step1res_sen$b2=step2res_sen$b
step1res_sen$se2=step2res_sen$se
step1res_sen$p2=step2res_sen$p
step1res_sen$nsnp2=step1res_sen$nsnp
colnames(step1res)[6] <- "nsnp1"
colnames(step1res)[7] <- "b1"
colnames(step1res)[8] <- "se1"
colnames(step1res)[9] <- "p1"
colnames(step1res_sen)[6] <- "nsnp1"
colnames(step1res_sen)[7] <- "b1"
colnames(step1res_sen)[8] <- "se1"
colnames(step1res_sen)[9] <- "p1"
step1res <- step1res[,-(11:12)]
step1res <- step1res[,-(2:5)]
step1res_sen <- step1res_sen[,-(11:12)]
step1res_sen <- step1res_sen[,-(2:5)]
res_2smr <- res_2smr[,-(10:13)]
res_2smr <- res_2smr[,-(2:5)]
res_2smr_sen <- res_2smr_sen[,-(10:12)]
res_2smr_sen <- res_2smr_sen[,-(2:5)]

twostep <- left_join(step1res,res_2smr,by ='outcome')
twostep_sen <- left_join(step1res_sen,res_2smr_sen,by ='outcome')
twostep$b3 <- twostep$b1*twostep$b2
twostep_sen$b3 <- twostep_sen$b1*twostep_sen$b2
twostep$se3 = sqrt(twostep$b1*twostep$b1*twostep$se2*twostep$se2+twostep$b2*twostep$b2*twostep$se1*twostep$se1)
twostep_sen$se3 = sqrt(twostep_sen$b1*twostep_sen$b1*twostep_sen$se2*twostep_sen$se2+twostep_sen$b2*twostep_sen$b2*twostep_sen$se1*twostep_sen$se1)
twostep$Z3 = twostep$b3/twostep$se3
twostep_sen$Z3 = twostep_sen$b3/twostep_sen$se3
twostep$p3 <-  2*pnorm(q=abs(twostep$Z3), lower.tail=FALSE)
twostep_sen$p3 <-  2*pnorm(q=abs(twostep_sen$Z3), lower.tail=FALSE)
twostep <- subset(twostep,twostep$p2 <0.05)
twostep_sen <- subset(twostep_sen,twostep_sen$p2 <0.05)
twostep$Proportion_mediated = twostep$b3/twostep$b
for (j in 1:nrow(twostep)){
  if(twostep$Proportion_mediated[j]>0){
    twostep$Directionally_consistent[j]= 'Y'
  }else twostep$Directionally_consistent[j]= 'N'
}
for (j in 1:nrow(twostep)){
  if(twostep$Proportion_mediated[j] < 0 ){
    twostep$Proportion_mediated[j]= 'NA'}
}

twostep_sen$Proportion_mediated = twostep_sen$b3/twostep_sen$b
for (j in 1:nrow(twostep_sen)){
  if(twostep_sen$Proportion_mediated[j]>0){
    twostep_sen$Directionally_consistent[j]= 'Y'
  }else twostep_sen$Directionally_consistent[j]= 'N'
}
for (j in 1:nrow(twostep_sen)){
  if(twostep_sen$Proportion_mediated[j] < 0 ){
    twostep_sen$Proportion_mediated[j]= 'NA'}
}

write.csv(twostep,"result/meditated.csv",row.names = FALSE)
write.csv(twostep_sen,"result/meditated_sen.csv",row.names = FALSE)
#### END: mediated ----
