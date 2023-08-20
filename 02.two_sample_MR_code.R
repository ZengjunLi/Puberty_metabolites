library(TwoSampleMR)
library(dplyr)
library(tidyr)
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO",force = TRUE)
library(MRPRESSO)
setwd("pathname/data")

#### 2SMR PT on Met ----
#exposure
day208 <- read.csv("exposure_and_confounders/day208.csv")
colnames(day208)

exposure_dat <- read_exposure_data(
  filename = "day208.csv", 
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
#outcome
met <- read.csv("outcome/all_pt_metabolites.csv")
colnames(met)

outcome_dat <- read_outcome_data(
  filename = "outcome/all_pt_metabolites.csv",
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


# # Using rs11227247 as a proxy to replace rs10750766 r2=0.35 delete

h <- harmonise_data(exposure_dat, outcome_dat)
table(h$outcome)
h2 <- h[h$mr_keep,]
table(h2$outcome)
write.csv(h2,"exposure_and_confounders/presso/h.csv",row.names = FALSE)

# MR
res <- mr(h2, method_list=c('mr_ivw_mre',
                            'mr_egger_regression',
                            'mr_weighted_median',
                            'mr_weighted_mode'))
pathway <- read.csv("pathway.csv")
res <- merge(res,pathway,by="outcome")
res_2smr_pt_met_ivw <- res[res$method == "Inverse variance weighted (multiplicative random effects)",]
res_2smr_pt_met_ivw <- merge(res_2smr_pt_met_ivw,pathway,by="outcome")
res_2smr_pt_met_ivw$lower = res_2smr_pt_met_ivw$b-1.96*res_2smr_pt_met_ivw$se
res_2smr_pt_met_ivw$upper = res_2smr_pt_met_ivw$b+1.96*res_2smr_pt_met_ivw$se


#heterogeneity   pleiotropy
hete_met <- mr_heterogeneity(h2)
hete_met <- merge(hete_met,pathway,by="outcome")
pleio_met <- mr_pleiotropy_test(h2)
pleio_met <- merge(pleio_met,pathway,by="outcome")

write.csv(res,"result/pt_met_2smr_results.csv",row.names = FALSE)
write.csv(res_2smr_pt_met_ivw,"result/pt_met_2smr_results(ivw).csv",row.names = FALSE)
write.csv(hete_met,"result/pt_met_2smr_heterogeneity.csv",row.names = FALSE)
write.csv(pleio_met,"result/pt_met_2smr_pleiotropy.csv",row.names = FALSE)

#### End: 2SMR PT on Met ----

#### MR PRESSO ----
h <- read.csv("exposure_and_confounders/presso/h.csv")

first_category_name = list.files("C:/Users/admin/Desktop/metabolites2/") 
h[,18] <- sub(":", "_", h[,18], fixed = TRUE)
full <- unique(h$outcome)
for (circle in full) {
  file_name = paste(circle,'.csv',sep = '')
  file_path = paste('exposure_and_confounders/presso/presso_metabolites/',file_name,sep = '')
  file_full <- h %>% filter(outcome == circle)
  write.csv(file_full,file_path)
}



first_category_name = list.files("exposure_and_confounders/presso/presso_metabolites/") 
for (i in 1 : 174) {
  meta <- read.csv(paste0('exposure_and_confounders/presso/presso_metabolites/',
                          first_category_name[i]))
  
  presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure",
                      SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                      OUTLIERtest = TRUE,DISTORTIONtest = TRUE, 
                      data = meta, 
                      NbDistribution = 300,  
                      SignifThreshold = 0.05)
  presso2 <- data.frame(new=c(1,2))
  presso2$MRAnalysis = presso[["Main MR results"]][["MR Analysis"]]
  presso2$causalestimate = presso[["Main MR results"]][["Causal Estimate"]]
  presso2$sd =presso[["Main MR results"]][["Sd"]]
  presso2$Pvalue = presso[["Main MR results"]][["P-value"]]
  presso2$GlobalTestRSSobs = presso[["MR-PRESSO results"]][["Global Test"]][["RSSobs"]]
  presso2$GlobalTestPvalue = presso[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
  presso2$outcome = first_category_name[i]
  if (i == 1){
    pressoall <- presso2
  }else{
    pressoall <- rbind(pressoall,presso2)  
  }
  
}

pressoall[,8] <- sub("_", ":", pressoall[,8], fixed = TRUE)
write.csv(pressoall,"result/pressoall.csv", row.names = FALSE)
#### End: MR PRESSO----

#### 2SMR PT on Met (Excluding 7 snps) ----
day201 <- read.csv("exposure_and_confounders/day201.csv")
colnames(day201)

exposure_dat_sensitivity <- read_exposure_data(
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

#outcome
met <- read.csv("outcome/all_pt_metabolites.csv")
colnames(met)

outcome_dat <- read_outcome_data(
  filename = "outcome/all_pt_metabolites.csv",
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


# # Using rs11227247 as a proxy to replace rs10750766 r2=0.35 delete
h_sen <- harmonise_data(exposure_dat_sensitivity, outcome_dat)
table(h_sen$outcome)
h_sen2 <- h_sen[h_sen$mr_keep,]
table(h_sen2$outcome)

res_sen_2smr <- mr(h_sen2, method_list=c('mr_ivw_mre'))
res_sen_2smr <- merge(res_sen_2smr,pathway,by="outcome")
res_sen_2smr$lower = res_sen_2smr$b-1.96*res_sen_2smr$se
res_sen_2smr$upper = res_sen_2smr$b+1.96*res_sen_2smr$se
res_sen_2smr <- merge(res_sen_2smr,pathway,by="outcome")

#heterogeneity   pleiotropy
hete_met2 <- mr_heterogeneity(h_sen2)
hete_met2 <- merge(hete_met2,pathway,by="outcome")
pleio_met2 <- mr_pleiotropy_test(h_sen2)
pleio_met2 <- merge(pleio_met2,pathway,by="outcome")

write.csv(hete_met2,"result/pt_met_2smr_sensitivity_heterogeneity.csv",row.names = FALSE)

write.csv(pleio_met2,"result/pt_met_2smr_sensitivity_pleiotropy.csv",row.names = FALSE)

write.csv(res_sen_2smr,"result/pt_met_2smr_sensitivity_IVW.csv", row.names = FALSE)
#### End: 2SMR PT on Met (Excluding 7 snps)----





