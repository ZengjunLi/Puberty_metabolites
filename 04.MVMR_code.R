library(TwoSampleMR)
setwd("pathname/data")
library(dplyr)


#### MVMR ----
# exposure
filenames_exposure=c("exposure_and_confounders/all_pt_mvmr.csv",
                     "exposure_and_confounders/all_abmi_mvmr.csv")
exposure_dat <- mv_extract_exposures_local(
  filenames_exposure, 
  sep = ",",
  phenotype = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  pval_col = "pval",
  se_col = "se",
  samplesize_col = "N",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "eaf" ,  pval_threshold = 5e-08,
  clump_r2 = 0.001,
  clump_kb = 10000,
  harmonise_strictness = 2
)

# outcome
first_category_name = list.files("D:/research/outcome/mvmr_met_sep/") 
for (i in 1 : 174) {
  data_ <- read_outcome_data(
    filename = paste0('D:/research/outcome/mvmr_met_sep/',first_category_name[i]),
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
# MR
  h <- mv_harmonise_data(exposure_dat, data_, harmonise_strictness = 2)
  result <- data.frame(mv_multiple(h))
  
  if (i == 1){
    metab <- result
    
  }else{
    metab <- rbind(metab,result)  
    
  }
}

metab[,4] <- sub("_", ":", metab[,4], fixed = TRUE)
colnames(metab)[4] <- "outcome"
# pathway
pathway <- read.csv("pathway.csv")
metab_full <- merge(metab,pathway,by="outcome")
metab_full$lower = metab_full$result.b-1.96*metab_full$result.se
metab_full$upper = metab_full$result.b+1.96*metab_full$result.se
write.csv(metab_full,"result/mvmr_res.csv", row.names = FALSE)
#### END: MVMR ----



#### MVMR sensitivity----
# exposure
filenames_exposure2=c("exposure_and_confounders/sen_pt_mvmr.csv",
                      "exposure_and_confounders/sen_abmi_mvmr.csv")
sen_exposure_dat <- mv_extract_exposures_local(
  filenames_exposure2, 
  sep = ",",
  phenotype = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  pval_col = "pval",
  se_col = "se",
  samplesize_col = "N",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "eaf" ,  pval_threshold = 5e-08,
  clump_r2 = 0.001,
  clump_kb = 10000,
  harmonise_strictness = 2
)
# outcome
first_category_name = list.files("D:/research/outcome/mvmr_met_sep/") 
for (i in 1 : 174) {
  data_ <- read_outcome_data(
    filename = paste0('D:/research/outcome/mvmr_met_sep/',first_category_name[i]),
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
# MR
  h_sen <- mv_harmonise_data(sen_exposure_dat, data_, harmonise_strictness = 2)
  result_sen <- data.frame(mv_multiple(h_sen))
  
  if (i == 1){
    metab_sen <- result_sen
    
  }else{
    metab_sen <- rbind(metab_sen,result_sen)  
    
  }
}


metab_sen[,4] <- sub("_", ":", metab_sen[,4], fixed = TRUE)
colnames(metab_sen)[4] <- "outcome"
# pathway

metab_sen_full <- merge(metab_sen,pathway,by="outcome")
metab_sen_full$lower = metab_sen_full$result.b-1.96*metab_sen_full$result.se
metab_sen_full$upper = metab_sen_full$result.b+1.96*metab_sen_full$result.se

write.csv(metab_sen_full,"result/mvmr_sen_res.csv", row.names = FALSE)

#### END: MVMR ----
