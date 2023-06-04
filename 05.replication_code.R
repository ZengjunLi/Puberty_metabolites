install.packages("wesanderson")
remotes::install_github("sritchie73/ukbnmr")
install.packages("patchwork")
library(patchwork)
library(TwoSampleMR)
library(tidyverse)
library(wesanderson)
library(ukbnmr)
library(data.table)
setwd("pathname/data")
library(scales)
show_col(hue_pal()(6))


#### Two-sample MR : 9 amini acids replication analysis ----
# IEU Open GWAS: https://gwas.mrcieu.ac.uk/
pt <- "ieu-a-1095"
confounders <- c("ukb-b-13378", "ukb-b-4650", "ukb-a-34", "ieu-a-1083", "ebi-a-GCST90002409")
aa.gwas.id <- paste0("met-d-",nmr_info %>% filter(Group == "Amino acids") %>% select(Biomarker) %>% unlist)
# replication Two-sample MR 
exp_dat <- extract_instruments(aa.gwas.id)
write.csv(exp_dat,"replication/aa_tophits.csv",row.names = FALSE)
exp_dat <- extract_instruments(pt)
out_dat <- extract_outcome_data(exp_dat$SNP, aa.gwas.id)
dat3 <- harmonise_data(exp_dat, out_dat)
res3 <- mr(dat3)
res3$lower = res3$b-1.96*res3$se
res3$upper = res3$b+1.96*res3$se
res3 <- split_outcome(res3)
res3[which(res3$outcome == "Total concentration of branched-chain amino acids (leucine + isoleucine + valine)"), ]$outcome <- c("Total BCAA")
write.csv(res3,"replication/2smr_pt_aa_res.csv",row.names = FALSE)

#### End: Two-sample MR : 9 amini acids replication analysis ----

#### MVMR : 9 amini acids replication analysis ----
# IEU Open GWAS: https://gwas.mrcieu.ac.uk/
exp_dat <- mv_extract_exposures(c("ieu-a-1095", "ieu-a-835"))
mvmr_res <- list()
for (i in 1:10) {
  aa.gwas=aa.gwas.id[[i]]
  out_dat <- extract_outcome_data(exp_dat$SNP, aa.gwas)
  dat3 <- mv_harmonise_data(exp_dat, out_dat)
  res3 <- mv_multiple(dat3)
  mvmr_res[[i]] <- res3$result
  
}

mvmr_res <- rbindlist(mvmr_res)
mvmr_res$lower = mvmr_res$b-1.96*mvmr_res$se
mvmr_res$upper = mvmr_res$b+1.96*mvmr_res$se
mvmr_res <- split_exposure(split_outcome(mvmr_res))
mvmr_res[which(mvmr_res$outcome == "Total concentration of branched-chain amino acids (leucine + isoleucine + valine)"), ]$outcome <- c("Total BCAA")
write.csv(mvmr_res,"replication/mvmr_pt_aa_res.csv",row.names = FALSE)
#### End: MVMR : 9 amini acids replication analysis ----