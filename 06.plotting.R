install.packages("ggplot2")
library(ggplot2)
install.packages("devtools")
install.packages("cli")
install.packages("lifecycle")
devtools::install_github("lionel-/ggstance")
devtools::install_github("thomasp85/ggforce")
library(dplyr)
library(tidyr)
setwd("D:/research")


####-- Two-sample MR ####
result_com <- read.csv("plot/plot4/2smr/polt.csv")
xname <- expression(paste("Standardized mean differences in metabolites per year later onset of puberty (95% CI)"))

# edit the low and up ci in a new dataset 
result_com2 <- result_com
result_com2$lower <- result_com2$beta-1.96*result_com2$se
result_com2$upper <- result_com2$beta+1.96*result_com2$se


result_com2$pathway[which(result_com2$pathway =='Lyso-phosphatidylcholines')] <- 'LPC'
result_com2$pathway[which(result_com2$pathway =='Amino acids')] <- 'AA'
result_com2$pathway[which(result_com2$pathway =='Phosphatidylcholines')] <- 'PC'
result_com2$pathway[which(result_com2$pathway =='Biogenic amines')] <- 'BA'
result_com2$pathway[which(result_com2$pathway =='Acylcarnitines')] <- 'AC'

result_com2$method[which(result_com2$method =='MR-Presso')] <- 'MR-PRESSO'
result_com2$method[which(result_com2$method =='IVW (confounders controlled)')] <- 'IVW (excluding 7 SNPs)'
result_com2$Method <- result_com2$method
result_com2$approach <- factor(result_com2$group, levels = c("IVW","Sensitivity analysis"))


result_com2 <- result_com2 %>% filter(outcome == "lysoPC a C20:4" | outcome=="lysoPC a C18:2"|
outcome =="lysoPC a C18:1"|outcome =="lysoPC a C18:0"|outcome =="lysoPC a C17:0"|
  outcome =="lysoPC a C16:0"|outcome =="trans-Hydroxyproline"|outcome =="Sarcosine"|
  outcome =="Kynurenine"|outcome == "Tyrosine"|outcome =="Tryptophan"|outcome =="Threonine"|
  outcome =="Serine"|outcome =="Phenylalanine"|outcome =="Methionine"|outcome =="Lysine"|
  outcome =="Carnitine"|outcome =="Aspartate"|outcome =="Asparagine"|outcome =="Valerylcarnitine"|
  outcome =="Tigylcarnitine"|outcome =="Nonaylcarnitine"|outcome =="Hydroxytetradecenoylcarnitine"|
  outcome =="Hydroxyhexadecenoylcarnitine"|
  outcome =="Dodecenoylcarnitine"|outcome =="Dodecanedioylcarnitine"|outcome =="Carnitine")

result_com2$method[which(result_com2$method =='MR-Presso')] <- 'MR-PRESSO'
result_com2$method[which(result_com2$method =='IVW (confounders controlled)')] <- 'IVW (excluding 7 SNPs)'


tsmr <-  ggplot2::ggplot(data=result_com2, aes(y=outcome, x=beta, xmin=lower, xmax=upper)) +
  geom_errorbar(aes(xmin=lower, xmax=upper,col= Method),width=0.1,cex=1, position = ggstance::position_dodgev(height = 0.5))+
  
  geom_pointrange(aes(col= Method), size=0.5, position = ggstance::position_dodgev(height = 0.5) )+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.4,
    colour = "black"
  ) +
  
  scale_x_continuous(name=xname, breaks = c(-0.4,-0.3,-0.2, -0.1, 0, 0.1, 0.2, 0.3,0.4 )) +
  ggplot2::coord_cartesian(xlim = c(-0.4, 0.4)) +
  
  ggforce::facet_grid_paginate (
    facets = pathway ~ approach,
    scales = "free_y",
    space = "free",
    switch = "y"
  ) +
  
  theme_bw() +

  
  #thematic stuff
  theme(strip.text.y = element_text(size = 35, face = "bold"),
        strip.text.x = element_text(size = 30, face = "bold"),
        
        legend.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 30),
        
        axis.text.y = element_text(size = 35),
        axis.text.x = element_text(size = 30),
        
        axis.title.y = element_blank(),
        # axis.title.y = element_text(size = 25, face = "bold" ),
        axis.title.x = element_text(size = 30, face = "bold" ),
        
        panel.grid.minor.x = element_blank())
colors_6<-c(	"#3cb346","#00abf0","#d75427",
             "#2e409a","#942d8d","#eeb401")

tsmr+
  scale_colour_manual(values = colors_6)



ggsave(file="C:/Users/admin/Desktop/mr/material and methods/20230804/ld2.tiff", width = 35, height = 40, dpi = 100)
ggsave(file="C:/Users/admin/Desktop/mr/material and methods/20230829/20230909/FIG_PDF/FIG4.PDF", width = 35, height = 40, dpi = 600)
ggsave(file="plot/plot4/2smr/ld.tiff", width = 35, height = 25, dpi = 100)
####-- End: two-sample MR ####

####-- MVMR ####
result_mvmr <- read.csv("plot/plot4/MVMR/mvmr_plot.csv")
result_mvmr <- result_mvmr %>% filter(outcome == "lysoPC a C20:4" | outcome=="PC aa C40:4"|
                                        outcome =="trans-Hydroxyproline"|outcome =="Ornithine"|
                                        outcome =="Tigylcarnitine"|outcome =="Lysine"|
                                        outcome =="Octadecenoylcarnitine"|
                                        outcome =="Hydroxytetradecenoylcarnitine"|outcome =="Dodecanedioylcarnitine"|outcome =="Carnitine")


result_mvmr$result.exposure[which(result_mvmr$result.exposure =='Body mass index || id:ebi-a-GCST006368')] <- 'Adulthood BMI'
result_mvmr$result.exposure[which(result_mvmr$result.exposure =='pt')] <- 'Puberty timing'
result_mvmr$group[which(result_mvmr$group =='MVMR (confounders controlled)')] <- 'MVMR (excluding 8 SNPs)'

result_mvmr$pathway[which(result_mvmr$pathway =='Lyso-phosphatidylcholines')] <- 'LPC'
result_mvmr$pathway[which(result_mvmr$pathway =='Amino acids')] <- 'AA'
result_mvmr$pathway[which(result_mvmr$pathway =='Phosphatidylcholines')] <- 'PC'
result_mvmr$pathway[which(result_mvmr$pathway =='Biogenic amines')] <- 'BA'
result_mvmr$pathway[which(result_mvmr$pathway =='Acylcarnitines')] <- 'AC'

xname <- expression(paste("Standardized mean differences in metabolites per year later onset of puberty (95% CI)"))

result_mvmr$approach <- factor(result_mvmr$result.exposure, levels = c("Puberty timing","Adulthood BMI"))
colnames(result_mvmr)[10] <- "Method"

mvmr <-  ggplot2::ggplot(data=result_mvmr, aes(y=outcome, x=result.b, xmin=lower, xmax=upper)) +
  geom_errorbar(aes(xmin=lower, xmax=upper,col= Method),width=0.1,cex=1, position = ggstance::position_dodgev(height = 0.5))+
  
  geom_pointrange(aes(col= Method), size=0.5, position = ggstance::position_dodgev(height = 0.5) )+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.4,
    colour = "black"
  ) +
  
  scale_x_continuous(name=xname, breaks = c(-0.4,-0.3,-0.2, -0.1, 0, 0.1, 0.2, 0.3,0.4 )) +
  ggplot2::coord_cartesian(xlim = c(-0.4, 0.4)) +
  
  ggforce::facet_grid_paginate (
    facets = pathway ~ approach,
    scales = "free_y",
    space = "free",
    switch = "y"
  ) +
  theme_bw() +
  theme(strip.text.y = element_text(size = 20, face = "bold"),
        strip.text.x = element_text(size = 23, face = "bold"),
        
        legend.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 23),
        
        axis.text.y = element_text(size = 25),
        axis.text.x = element_text(size = 23),
        
        axis.title.y = element_blank(),
        # axis.title.y = element_text(size = 25, face = "bold" ),
        axis.title.x = element_text(size = 30, face = "bold" ),
        
        panel.grid.minor.x = element_blank())
mvmr
colors_5<-c("#4da0a0","#9b3a74")
mvmr+scale_colour_manual(values = colors_5)
ggsave(file="plot/plot4/MVMR/mvmr.tiff", width = 30, height = 10, dpi = 100)
ggsave(file="C:/Users/admin/Desktop/mr/material and methods/20230829/20230909/FIG_PDF/FIG5.PDF", width = 30, height = 10, dpi = 600)
####-- End:MVMR ####

####-- Replication analysis ####
library(patchwork)
library(TwoSampleMR)
library(tidyverse)
library(wesanderson)
library(ukbnmr)
library(data.table)
setwd("D:/research")
library(scales)
show_col(hue_pal()(6))
# 2smr
pt_2smr <- read.csv("result/pt_met_2smr_results(ivw).csv")
MVMR <- read.csv("result/mvmr_res.csv")

pt_2smr <- filter(pt_2smr, pathway.x == "Amino acids")
pt_2smr$group <- "Two sample MR"
pt_2smr <- pt_2smr[,-(2:5)]

pt_2smr <- pt_2smr[,-7]
colnames(pt_2smr)[6] <- "pathway"
MVMR <- filter(MVMR, pathway == "Amino acids")
MVMR <- filter(MVMR, result.exposure == "pt")
MVMR$group <- "MVMR"
MVMR <- MVMR[,-3]
colnames(MVMR)[4] <- "b"
colnames(MVMR)[5] <- "se"
colnames(MVMR)[6] <- "pval"

res3 <- read.csv("replication/2smr_pt_aa_res.csv")
for (j in 1:nrow(res3)){
  if(res3$pval[j] < 0.05 ){
    res3$outcome[j]= paste('*',res3$outcome[j],sep = '')
  }else {
    res3$outcome[j]= res3$outcome[j]
  }
}


res3 <- filter(res3, method == "Inverse variance weighted")
res3 <- res3[,-(4:5)]
res3 <- res3[,-(1:2)]
res3$group <- "Two sample MR (replication in UKB)"
res3$pathway <- "Amino acids"

mvmr_res <- read.csv("replication/mvmr_pt_aa_res.csv")
mvmr_res <- filter(mvmr_res, exposure == "Age at menarche")
mvmr_res <- mvmr_res[,-(1:3)]
mvmr_res$group <- "MVMR (replication in UKB)"
mvmr_res$pathway <- "Amino acids"

data <- rbind(pt_2smr,MVMR,res3,mvmr_res)
write.csv(data, "F:/data.csv", row.names = FALSE)
data <- read.csv("F:/data.csv")




data <- data %>% filter(outcome=="Valine"|
                          outcome =="Tyrosine"|outcome =="Phenylalanine"|
                          outcome =="Leucine"|outcome =="Isoleucine"|
                          outcome =="Histidine"|
                          outcome =="Glycine"|outcome =="Glutamine"|outcome =="Alanine")
data$approach[which(data$group =='MVMR (replication in UKB)')] <- 'MVMR'
data$approach[which(data$group =='Two sample MR')] <- 'Two sample MR'
data$approach[which(data$group =='Two sample MR (replication in UKB)')] <- 'Two sample MR'
data$approach[which(data$group =='MVMR')] <- 'MVMR'


xname <- expression(paste("Standardized mean differences in metabolites per year later onset of puberty (95% CI)"))
data$approach <- factor(data$approach, levels = c("Two sample MR","MVMR"))
rep <-  ggplot2::ggplot(data=data, aes(y=outcome, x=b, xmin=lower, xmax=upper)) +
  geom_errorbar(aes(xmin=lower, xmax=upper,col= group),width=0.1,cex=1, position = ggstance::position_dodgev(height = 0.5))+
  
  geom_pointrange(aes(col= group), size=0.5, position = ggstance::position_dodgev(height = 0.5) )+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.4,
    colour = "black"
  ) +
  
  scale_x_continuous(name=xname, breaks = c(-0.1, 0, 0.1)) +
  ggplot2::coord_cartesian(xlim = c(-0.1, 0.1)) +
  ggforce::facet_grid_paginate (
    facets = pathway ~ approach,
    scales = "free_y",
    space = "free",
    switch = "y"
  ) +
  theme_bw() +
  theme(strip.text.y = element_text(size = 15, face = "bold"),
        strip.text.x = element_text(size = 15, face = "bold"),
        
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15),
        
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        
        axis.title.y = element_blank(),
        # axis.title.y = element_text(size = 25, face = "bold" ),
        axis.title.x = element_text(size = 18, face = "bold" ),
        
        panel.grid.minor.x = element_blank())

rep

colors_4<-c("#4da0a0","#9b3a74","#156077","#f46f20")

rep+scale_colour_manual(values = colors_4)

rep

ggsave( "plot/plot4/replication/Fig 5. Forest plots comparing the effect estimates of puberty timing on 9 amino acids in replication analyses..tiff",
        width = 30, height =13, units = "cm", dpi = 300)

ggsave( "C:/Users/admin/Desktop/mr/material and methods/20230829/20230909/FIG_PDF/Fig 6. Forest plots comparing the effect estimates of puberty timing on 9 amino acids in replication analyses..PDF",
        width = 30, height =13, units = "cm", dpi = 600)


####--####
# Install devtools
library(usethis)
install.packages("devtools")
library(devtools)

# Install EpiViz directly from GitHub
devtools::install_github("mattlee821/EpiViz/R_package")
library(EpiViz)
library(circlize) 
library(dplyr)
library(epivizr)
library(BiocManager)
install.packages("devtools")

install.packages("stringi")

install.packages("usethis")
library(usethis)
library(devtools)
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap",force = T)
devtools::install_github("mattlee821/EpiViz/R_package")
library(EpiViz)
#### tutorial data ----
dt <- EpiViz::EpiViz_data1 %>% arrange(subclass, label)
dt2 <- dt %>% arrange(subclass, label) %>%
  distinct(subclass) %>%
  mutate(subclass.number = row_number(), legend = paste0(row_number(),".",subclass)) 
dt3 <- dt %>% group_by(subclass) %>%
  summarise(n= n()) %>%
  mutate(subclass.distance = 1/(n+1))
dt <- dt %>% left_join(dt2, by = c("subclass" = "subclass")) %>%
  left_join(dt3, by = c("subclass" = "subclass")) %>%
  mutate(label.number = row_number())

dt$subclass.number <- factor(dt$subclass.number)

a <- 1
b <- 1
for (i in 1:123) {
  if (dt$subclass.number[i] == a){
    dt$subclass.order[i] = dt$label.number[i] - b + 1
  }else {
    a <- dt$subclass.number[i]
    b <- dt$label.number[i]
    dt$subclass.order[i] = dt$label.number[i] - b + 1
  }
}

color <- colorRampPalette(c("#86E3CE", "#D0E6A5", "#FFDD94", "#FA897B", "#CCABD8"))(28) # change to number of sector
# set reference line
ref <- 0 
#### END: tutorial data ----

#### PT & Met. project 2SMR ----
category <- read.csv("D:/Research/pathway.csv")
colnames(category) <- c("label","subclass")

pt_met_2smr <- read.csv("D:/research/result/pt_met_2smr_results(ivw).csv")
dat <- pt_met_2smr[pt_met_2smr$method =="Inverse variance weighted (multiplicative random effects)",c(1,7:9)]
colnames(dat) <- c("label","beta","se","pvalue")
dat$lower_confidence_interval <- dat$beta - 1.96*dat$se
dat$upper_confidence_interval <- dat$beta + 1.96*dat$se
dat <- left_join(dat, category, by = "label")
dat <- dat %>% arrange(subclass, label)


dat2 <- dat %>% arrange(subclass, label) %>%
  distinct(subclass) %>%
  mutate(subclass.number = row_number(), legend = paste0(row_number(),".",subclass)) 
dat3 <- dat %>% group_by(subclass) %>%
  summarise(n= n()) %>%
  mutate(subclass.distance = 1/(n+1))

dat <- dat %>% left_join(dat2, by = c("subclass" = "subclass")) %>%
  left_join(dat3, by = c("subclass" = "subclass")) %>%
  mutate(label.number = row_number())

dat$subclass.number <- factor(dat$subclass.number)

a <- 1
b <- 1
for (i in 1:174) {
  if (dat$subclass.number[i] == a){
    dat$subclass.order[i] = dat$label.number[i] - b + 1
  }else {
    a <- dat$subclass.number[i]
    b <- dat$label.number[i]
    dat$subclass.order[i] = dat$label.number[i] - b + 1
  }
}

dt <- dat
dt2 <- dat2
dt3 <- dat3

color <- colorRampPalette(c("#86E3CE", "#D0E6A5", "#FFDD94", "#FA897B", "#CCABD8"))(7) # change to number of sector
# set reference line
ref <- 0 
#### END: PT & Met. project 2SMR ----

pdf("D:/research/result/pt__met_2smr.pdf",
    width = 20, height = 15)

graphics::par(mar = c(0.5, 0.5, 0.5, 0.5) * 20, 
              cex = 1, xpd = NA)
# ?circos.par
circos.par(start.degree = 90, 
           gap.degree = c(rep(1,6), 15), # change to 1 to number of sector - 1
           cell.padding = c(0.02, 0, 0.02, 0),
           points.overflow.warning = FALSE,
           track.margin = c(0.012, 0.012),
           track.height = 0.2
)
circos.initialize(factors = dt$subclass.number, xlim = c(0,1), 
                  sector.width = dt3$n)

# track1 - secter & label
circos.track(factors = dt$subclass.number, ylim = c(0,1), 
             track.height = 0.1, 
             bg.col = color,
             bg.border = NA,
             
             panel.fun = function(x, y) {
               # sector names inside the track     
               circos.text(CELL_META$xcenter, 
                           CELL_META$ycenter,
                           CELL_META$sector.index, facing = "bending.inside", niceFacing = TRUE,
                           adj = c(0.5, 0.25), cex = 0.8) 
             }
)

# labels outside the track
circos.trackText(dt$subclass.number,
                 dt$subclass.order*dt$subclass.distance,
                 dt$pvalue*0 + 2,
                 labels = dt$label, facing = "clockwise", niceFacing = TRUE,
                 adj = c(0, 0.5), cex = 0.8, col = ifelse(dt$lower_confidence_interval > ref, "#66CCCC",
                                                          ifelse(dt$upper_confidence_interval < ref, "#B22222", "black"))
)

# track2 - point 
circos.track(factors = dt$subclass.number, 
             ylim = c(min(dt$lower_confidence_interval) - 0.01,max(dt$upper_confidence_interval) + 0.01), 
             track.height = 0.5, 
             bg.col = "white",
             bg.border = NA)

# beta
circos.trackPoints(dt$subclass.number,
                   dt$subclass.order*dt$subclass.distance, 
                   dt$beta, 
                   cex = 0.7, pch = 20,
                   col = ifelse(dt$lower_confidence_interval > ref, "#66CCCC",
                                ifelse(dt$upper_confidence_interval < ref, "#B22222", "grey")) 
)

# errorbar
for (i in 1:nrow(dt)) {
  dt1 = subset(dt, dt$subclass.number == i)
  # ?circos.segments
  circos.segments(x0 = dt1$subclass.order*dt1$subclass.distance, 
                  x1 = dt1$subclass.order*dt1$subclass.distance, 
                  y0 = dt1$subclass.order * 0 + dt1$lower_confidence_interval, 
                  y1 = dt1$subclass.order * 0 + dt1$upper_confidence_interval, 
                  col = ifelse(dt1$lower_confidence_interval > ref, "#66CCCC",
                               ifelse(dt1$upper_confidence_interval < ref, "#B22222", "grey")),
                  lwd = 1, sector.index = i)
}

# reference line
circos.trackLines(dt$subclass.number,
                  dt$subclass.order*dt$subclass.distance,
                  dt$subclass.order * 0 - ref,
                  col = "black",
                  lwd = 0.1, lty = 1)

# circos.yaxis
circos.yaxis(side = "left", 
             at = c(round(min(dt$lower_confidence_interval),2), 0, round(max(dt$upper_confidence_interval),2)), 
             labels.cex = 0.8, sector.index = 1, track.index = 2, 
             labels.col = "black", col = "black")
circos.trackText(1, 0, 0,
                 "(β)",
                 cex = 0.8,
                 adj = c(1.1,2))

# legend
# legend(x = 1.1, y = -0.8, pch = 15, col = c("#66CCCC","#B22222"), 
#        legend = c("green beta XXXX","red beta XXXX"), cex = 0.8, 
#        box.col = "black", ncol = 1, text.col = "black",
#        title = "  beta", title.col = "black", title.adj = 0)

legend(x = 1.1, y = -1.1, pch = 15, col = color, legend = dt2$legend, cex = 0.8,
       box.col = "black", ncol = 2, text.col = "black",
       title = "  Metabolites", title.col = "black", title.adj = 0)

circos.clear()

dev.off()

