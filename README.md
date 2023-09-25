# Puberty on plasma metabolites: a two-sample MR analysis
<br>
This repository accompanies the paper titled "Lifelong impacts of puberty timing on human plasma metabolic profiles."

This study employed MR (Mendelian randomization) analysis to explore the causal effects of puberty timing on plasma metabolites.

## Environment details

This study was carried out using the following language: R.4.1.3. Some of the data processing work was carried out using Linux.

## Analysis overview

We used a two-sample Mendelian randomization (MR) analysis to assess the total effect of puberty timing on 174 plasma metabolites. We utilized a two-step MR and MVMR (multivariable Mendelian randomization) to estimate the direct and indirect effects considering adulthood adiposity as a mediator. Single nucleotide polymorphisms (SNPs) related to confounding factors, birth weight, and childhood adiposity, were further excluded, which provided evidence for the robustness of our results. Replication analyses were also conducted using UK Biobank data to examine the effects of puberty timing on nine amino acids.

## Data sources

The genetic instruments for puberty timing were obtained from the corresponding literature and https://reprogen.org/data_download.html. Data for plasma metabolites mGWAS can be found at www.omicscience.org. Full summary data on adulthood body mass index were obtained from the MRC Integrative Epidemiology Unit OpenGWAS project (https://gwas.mrcieu.ac.uk/). Full summary data on birth weight and childhood BMI were obtained from Early Growth Genetics Consortium (egg-consortium.org). The replication analyses were conducted using full summary data from UK Biobank (https://www.ukbiobank.ac.uk/).

## Analysis steps 

1. Data pre-processing for puberty timing and plasma metabolites. <br>
2. Two-sample MR: exposure:: puberty timing; outcome: plasma metabolites. <br>
3. Two-step MR: exposure:: puberty timing; outcome: plasma metabolites; mediator: adulthood BMI. <br>
4. MVMR: exposure:: puberty timing and adulthood BMI; outcome: plasma metabolites. <br>
5. Replication analyses using data from UK Biobank. <br>
