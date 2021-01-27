###################################################################################
## Title: Calculate mean total community abundance for insectivores
## User: Lisbeth Hordley 
## email: l.hordley@pgr.reading.ac.uk
## Date: January 2020
###################################################################################
options(scipen=999)

rm(list = ls())
library(dplyr)

## read in data
bbs_data1 <- read.csv("../Data/BBS_2004_2018_invert_interpol_repeat_det.csv", header=TRUE)
length(unique(bbs_data$ENGLISH_NAME))

# calculate mean total community abundance

## total abundance for each site, year and iteration
BBS_tot <- bbs_data %>% group_by(GRIDREF,YEAR,i) %>% summarise(tot_abund = sum(TOT_det))
## calculate the mean abund over time for each site and iteration
BBS_tot2 <- BBS_tot %>% group_by(GRIDREF, i) %>% summarise(mean_abund=mean(tot_abund))
## calculate mean abundance across repitions
BBS_proxy <- BBS_tot2 %>% group_by(GRIDREF) %>% summarise(mean_abund=mean(mean_abund))
## save file
write.csv(BBS_proxy, file="../Data/Invert_mean_abund_det1.csv", row.names=FALSE)
