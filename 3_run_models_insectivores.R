###################################################################################
## Title: Run linear models (abund ~ FDis) for insectivores
## User: Lisbeth Hordley 
## email: l.hordley@pgr.reading.ac.uk
## Date: January 2020
###################################################################################

options(scipen=999)
library(dplyr)
library(purrr)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ape)
library(foreach)
library(doParallel)
rm(list = ls())

## read in data
FDis_abund <- read.csv("../Results/Insectivores/FDis_abund_det2_gow.csv", header=TRUE)
FDis_PA <- read.csv("../Results/Insectivores/FDis_PA_det2_gow.csv", header=TRUE)
mean_abund <- read.csv("../Data/Invert_mean_abund_det2.csv", header=TRUE)
BBS_site_dist <- read.csv("../Data/BBS_site_dist.csv", header=TRUE)

## create unique code to run each model on
## e.g. 6_4_2_0.5 (trait number: 6, trait number combination: 4, proportion combination: 2, prop. continuous traits: 0.5)
FDis_abund$unique_code = paste(FDis_abund$number_comb, FDis_abund$mix_comb, FDis_abund$p_cont, sep="_")
## same for PA data
FDis_PA$unique_code = paste(FDis_PA$number_comb, FDis_PA$mix_comb, FDis_PA$p_cont, sep="_")
## each unique code has 200 rows i.e. 200 sites

## create 100km blocking factor in each data frame to use as random factor
FDis_abund$GRIDREF_100km <- substr(FDis_abund$GRIDREF, 1, 2) ## take first two characters from gridref
FDis_PA$GRIDREF_100km <- substr(FDis_PA$GRIDREF, 1, 2) ## take first two characters from gridref
mean_abund$GRIDREF_100km <- substr(mean_abund$GRIDREF, 1, 2) ## take first two characters from gridref
mean_abund <- merge(mean_abund, BBS_site_dist, by="GRIDREF") ## so mean_abund has lat lon values for morans i test
## create distance matrix for morans i test
dists = as.matrix(dist(mean_abund[,c(4,5)])) ## using lat and lon values
dists_inv = 1/dists
diag(dists_inv) = 0

## check normality of response variable (total mean community abundance)
hist(mean_abund$mean_abund) ## right skew
mean_abund$log_mean_abund <- log(mean_abund$mean_abund) 
hist(mean_abund$log_mean_abund)

################################
######### Abundance ############
################################

## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
## change code below from %dopar% to %do% to switch from parallel to non-parallel
start_time <- Sys.time()
results_table_final <- foreach(i=unique(FDis_abund$unique_code),  .combine=rbind, .packages=c('lme4', 'dplyr', 'purrr', 'MuMIn', 'ape')) %dopar% {
  print(i)
  FDis_temp <- FDis_abund[FDis_abund$unique_code==i,]
  abund_fdis <- list(mean_abund, FDis_temp) %>% reduce(left_join, by = c("GRIDREF", "GRIDREF_100km")) ## join fdis and mean abund data together
  mod <- lmer(log_mean_abund ~ FDis + (1|GRIDREF_100km), data=abund_fdis) ## run model
  ## create dataframe with model residuals and site lat & lon values
  df = data.frame(res = residuals(mod), lat = abund_fdis$latitude, lon = abund_fdis$longitude)
  x <- Moran.I(df$res, dists_inv) ## calculate Moran's I
  ## save R2 and Moran's I p-value (to check for spatial autocorrelation)
  data.frame(R2=r.squaredGLMM(mod)[1], moran_p_value=x$p.value,
             trait_combination=i)
}
end_time <- Sys.time()
end_time - start_time ## nonparallel = ~5 hours, paralell = ~50 mins
## save file
write.csv(results_table_final, file="../Results/Insectivores/R2_abund_seed_det2_gow_lmer.csv", row.names=FALSE)

#######################################
######### Presence/Absence ############
#######################################

## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

start_time <- Sys.time()
results_table_final2 <- foreach (i=unique(FDis_PA$unique_code),  .combine=rbind, .packages=c('lme4', 'dplyr', 'purrr', 'MuMIn', 'ape')) %dopar% {
  print(i)
  FDis_temp <- FDis_PA[FDis_PA$unique_code==i,]
  pa_fdis <- list(mean_abund, FDis_temp) %>% reduce(left_join, by = c("GRIDREF", "GRIDREF_100km")) ## join fdis and mean abund data together
  mod <- lmer(log_mean_abund ~ FDis + (1|GRIDREF_100km), data=pa_fdis) ## run model
  ## create dataframe with model residuals and site lat & lon values
  df = data.frame(res = residuals(mod), lat = pa_fdis$latitude, lon = pa_fdis$longitude)
  x <- Moran.I(df$res, dists_inv) ## calculate Moran's I
  ## save R2 and Moran's I p-value (to check for spatial autocorrelation)
  ## save results
  data.frame(R2=r.squaredGLMM(mod)[1], moran_p_value=x$p.value,
             trait_combination=i)
}
end_time <- Sys.time()
end_time - start_time 

## save file
write.csv(results_table_final2, file="../Results/Insectivores/R2_PA_seed_det2_gow_lmer.csv", row.names=FALSE)












#### old code without blocking factor
##########################################################################################
## linear model morans i test

################################
######### Abundance ############
################################

moran_results_final <- NULL
## loops to calculate linear model mean_abund ~ FDis for each group of trait combination
for(i in unique(FDis_abund$unique_code)){
  FDis_temp <- FDis_abund[FDis_abund$unique_code==i,]
  abund_fdis <- list(mean_abund, FDis_temp) %>% reduce(left_join, by = c("GRIDREF", "GRIDREF_100km"))
  mod <- lm(log_mean_abund ~ FDis, data=abund_fdis)
  df = data.frame(res = residuals(mod), lat = abund_fdis$latitude, lon = abund_fdis$longitude)
  # coordinates(df)<-c('lon','lat')
  # var.mod<-variogram(res~1,data=df)
  # plot(var.mod) ## plot variogram if needed
  x <- Moran.I(df$res, dists_inv)
  moran_results <- data.frame(p_value=x$p.value,
                              trait_combination=i) 
  moran_results_final <- rbind(moran_results_final, moran_results)
}
## save file
write.csv(moran_results_final, file="../Results/Insectivores/moran_results_abund_noPCA_lm.csv", row.names=FALSE)

#######################################
######### Presence/Absence ############
#######################################
moran_results_final2 <- NULL
## loops to calculate linear model mean_abund ~ FDis for each group of trait combination
for(i in unique(FDis_PA$unique_code)){
  print(i)
  FDis_temp2 <- FDis_PA[FDis_PA$unique_code==i,]
  pa_fdis <- list(mean_abund, FDis_temp2) %>% reduce(left_join, by = c("GRIDREF", "GRIDREF_100km"))
  mod2 <- lm(log_mean_abund ~ FDis, data=pa_fdis)
  df = data.frame(res = residuals(mod2), lat = pa_fdis$latitude, lon = pa_fdis$longitude)
  # coordinates(df)<-c('lon','lat')
  # var.mod<-variogram(res~1,data=df)
  # plot(var.mod) ## plot variogram if needed
  x <- Moran.I(df$res, dists_inv)
  moran_results <- data.frame(p_value=x$p.value,
                              trait_combination=i) 
  moran_results_final2 <- rbind(moran_results_final2, moran_results)
}

## save file
write.csv(moran_results_final2, file="../Results/Insectivores/moran_results_PA_noPCA_lm.csv", row.names=FALSE)






###################################################################################################################
options(scipen=999)
library(purrr)
library(dplyr)

rm(list = ls())

## read in data
FDis_abund <- read.csv("../Results/Insectivores/FDis_abund_det2.csv", header=TRUE)
FDis_PA <- read.csv("../Results/Insectivores/FDis_PA_det2.csv", header=TRUE)
mean_abund <- read.csv("../Data/Invert_mean_abund_det2.csv", header=TRUE)
BBS_site_dist <- read.csv("../Data/BBS_site_dist.csv", header=TRUE)

## create unique code e.g. 6_4_2_0.5 (trait number: 6, trait number combination: 4, proportion combination: 2, prop. continuous traits: 0.5)
FDis_abund$unique_code = paste(FDis_abund$number_comb, FDis_abund$mix_comb, FDis_abund$p_cont, sep="_")
## same for PA data
FDis_PA$unique_code = paste(FDis_PA$number_comb, FDis_PA$mix_comb, FDis_PA$p_cont, sep="_")
## each unique code has 199 rows i.e. 199 sites

## create 100km blocking factor in each data frame
FDis_abund$GRIDREF_100km <- substr(FDis_abund$GRIDREF, 1, 2)
FDis_PA$GRIDREF_100km <- substr(FDis_PA$GRIDREF, 1, 2)
mean_abund$GRIDREF_100km <- substr(mean_abund$GRIDREF, 1, 2)
mean_abund <- merge(mean_abund, BBS_site_dist, by="GRIDREF") ## so mean_abund has lat lon values for morans i test
## create distance matrix for morans i test
dists = as.matrix(dist(mean_abund[,c(4,5)]))
dists_inv = 1/dists
diag(dists_inv) = 0

## check normality of response variable
hist(mean_abund$mean_abund) ## right skew
mean_abund$log_mean_abund <- log(mean_abund$mean_abund) 
hist(mean_abund$log_mean_abund)

################################
######### Abundance ############
################################
library(ape)
results_table_final <- NULL
moran_results_final <- NULL
## loops to calculate linear model mean_abund ~ FDis for each group of trait combination
for(i in unique(FDis_abund$unique_code[1])){
  FDis_temp <- FDis_abund[FDis_abund$unique_code==i,]
  abund_fdis <- list(mean_abund, FDis_temp) %>% reduce(left_join, by = c("GRIDREF", "GRIDREF_100km"))
  mod <- lmer(log_mean_abund ~ FDis + (1|GRIDREF_100km), data=abund_fdis)
  df = data.frame(res = residuals(mod), lat = abund_fdis$latitude, lon = abund_fdis$longitude)
  # coordinates(df)<-c('lon','lat')
  # var.mod<-variogram(res~1,data=df)
  # plot(var.mod) ## plot variogram if needed
  x <- Moran.I(df$res, dists_inv)
  moran_results <- data.frame(p_value=x$p.value,
                              trait_combination=i) 
  results_table_temp <- data.frame(R2=r.squaredGLMM(mod)[1] , 
                                   trait_combination=i) ## marginal R2 i.e. only the variation explained by FDis, not the blocking factor too
  moran_results_final <- rbind(moran_results_final, moran_results)
  results_table_final <- rbind(results_table_final, results_table_temp)
}
## save file
write.csv(results_table_final, file="../Results/Insectivores/R2_abund_det2_lmer.csv", row.names=FALSE)
write.csv(moran_results_final, file="../Results/Insectivores/moran_results_abund_lmer.csv", row.names=FALSE)

#######################################
######### Presence/Absence ############
#######################################
results_table_final2 <- NULL
moran_results_final2 <- NULL
## loops to calculate linear model mean_abund ~ FDis for each group of trait combination
for(i in unique(FDis_PA$unique_code)){
  FDis_temp2 <- FDis_PA[FDis_PA$unique_code==i,]
  pa_fdis <- list(mean_abund, FDis_temp2) %>% reduce(left_join, by = c("GRIDREF", "GRIDREF_100km"))
  mod2 <- lmer(log_mean_abund ~ FDis + (1|GRIDREF_100km), data=pa_fdis)
  df = data.frame(res = residuals(mod2), lat = pa_fdis$latitude, lon = pa_fdis$longitude)
  # coordinates(df)<-c('lon','lat')
  # var.mod<-variogram(res~1,data=df)
  # plot(var.mod) ## plot variogram if needed
  x <- Moran.I(df$res, dists_inv)
  moran_results <- data.frame(p_value=x$p.value,
                              trait_combination=i) 
  results_table_temp <- data.frame(R2=r.squaredGLMM(mod2)[1] , 
                                   trait_combination=i) ## marginal R2 i.e. only the variation explained by FDis, not the blocking factor too
  moran_results_final2 <- rbind(moran_results_final2, moran_results)
  results_table_final2 <- rbind(results_table_final2, results_table_temp)
}

## save file
write.csv(results_table_final2, file="../Results/Insectivores/R2_PA_det2_lmer.csv", row.names=FALSE)
write.csv(moran_results_final2, file="../Results/Insectivores/moran_results_PA_lmer.csv", row.names=FALSE)

##########################################################################################
## linear model morans i test

################################
######### Abundance ############
################################

moran_results_final <- NULL
## loops to calculate linear model mean_abund ~ FDis for each group of trait combination
for(i in unique(FDis_abund$unique_code)){
  FDis_temp <- FDis_abund[FDis_abund$unique_code==i,]
  abund_fdis <- list(mean_abund, FDis_temp) %>% reduce(left_join, by = c("GRIDREF", "GRIDREF_100km"))
  mod <- lm(log_mean_abund ~ FDis, data=abund_fdis)
  df = data.frame(res = residuals(mod), lat = abund_fdis$latitude, lon = abund_fdis$longitude)
  # coordinates(df)<-c('lon','lat')
  # var.mod<-variogram(res~1,data=df)
  # plot(var.mod) ## plot variogram if needed
  x <- Moran.I(df$res, dists_inv)
  moran_results <- data.frame(p_value=x$p.value,
                              trait_combination=i) 
  moran_results_final <- rbind(moran_results_final, moran_results)
}
## save file
write.csv(moran_results_final, file="../Results/Insectivores/moran_results_abund_lm.csv", row.names=FALSE)

#######################################
######### Presence/Absence ############
#######################################
moran_results_final2 <- NULL
## loops to calculate linear model mean_abund ~ FDis for each group of trait combination
for(i in unique(FDis_PA$unique_code)){
  print(i)
  FDis_temp2 <- FDis_PA[FDis_PA$unique_code==i,]
  pa_fdis <- list(mean_abund, FDis_temp2) %>% reduce(left_join, by = c("GRIDREF", "GRIDREF_100km"))
  mod2 <- lm(log_mean_abund ~ FDis, data=pa_fdis)
  df = data.frame(res = residuals(mod2), lat = pa_fdis$latitude, lon = pa_fdis$longitude)
  # coordinates(df)<-c('lon','lat')
  # var.mod<-variogram(res~1,data=df)
  # plot(var.mod) ## plot variogram if needed
  x <- Moran.I(df$res, dists_inv)
  moran_results <- data.frame(p_value=x$p.value,
                              trait_combination=i) 
  moran_results_final2 <- rbind(moran_results_final2, moran_results)
}

## save file
write.csv(moran_results_final2, file="../Results/Insectivores/moran_results_PA_lm.csv", row.names=FALSE)

############# check all moran results - does spatial autocorrelation disappear when the 100km blocking factor is added?
## i.e. compare between linear and mixed effects model moran results

## abundance no PCA
moran_results_abund_lm <- read.csv("../Results/Insectivores/moran_results_abund_noPCA_lm.csv",header=TRUE)
nrow(moran_results_abund_lm[moran_results_abund_lm$p_value<0.05, ])/nrow(moran_results_abund_lm)*100 ## 0% significant results
moran_results_abund_lmer <- read.csv("../Results/Insectivores/moran_results_abund_noPCA_lmer.csv",header=TRUE)
nrow(moran_results_abund_lmer[moran_results_abund_lmer$p_value<0.05, ])/nrow(moran_results_abund_lmer)*100 
## still 0% results show spatial autocorrelation

## P/A no PCA
moran_results_pa_lm <- read.csv("../Results/Insectivores/moran_results_PA_noPCA_lm.csv",header=TRUE)
nrow(moran_results_pa_lm[moran_results_pa_lm$p_value<0.05, ])/nrow(moran_results_pa_lm)*100 ## 0% significant results
moran_results_pa_lmer <- read.csv("../Results/Insectivores/moran_results_PA_noPCA_lmer.csv",header=TRUE)
nrow(moran_results_pa_lmer[moran_results_pa_lmer$p_value<0.05, ])/nrow(moran_results_pa_lmer)*100 
## still 0% results show spatial autocorrelation

## abundance with PCA
moran_results_abund_lm <- read.csv("../Results/Insectivores/moran_results_abund_lm.csv",header=TRUE)
nrow(moran_results_abund_lm[moran_results_abund_lm$p_value<0.05, ])/nrow(moran_results_abund_lm)*100 ## 0% significant results
moran_results_abund_lmer <- read.csv("../Results/Insectivores/moran_results_abund_lmer.csv",header=TRUE)
nrow(moran_results_abund_lmer[moran_results_abund_lmer$p_value<0.05, ])/nrow(moran_results_abund_lmer)*100 
## still 0% results show spatial autocorrelation

## P/A with PCA
moran_results_pa_lm <- read.csv("../Results/Insectivores/moran_results_PA_lm.csv",header=TRUE)
nrow(moran_results_pa_lm[moran_results_pa_lm$p_value<0.05, ])/nrow(moran_results_pa_lm)*100 ## 0% significant results
moran_results_pa_lmer <- read.csv("../Results/Insectivores/moran_results_PA_lmer.csv",header=TRUE)
nrow(moran_results_pa_lmer[moran_results_pa_lmer$p_value<0.05, ])/nrow(moran_results_pa_lmer)*100 
## still 0% results show spatial autocorrelation

## looks like presence-absence results might need a 50km blocking factor instead
