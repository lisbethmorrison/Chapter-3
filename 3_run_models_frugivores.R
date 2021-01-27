###################################################################################
## Title: Run linear models (abund ~ FDis) for frugivores
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
FDis_abund <- read.csv("../Results/Frugivores/FDis_abund_det2_gow.csv", header=TRUE)
FDis_PA <- read.csv("../Results/Frugivores/FDis_PA_det2_gow.csv", header=TRUE)
mean_abund <- read.csv("../Data/Seed_mean_abund_det2.csv", header=TRUE)
BBS_site_dist <- read.csv("../Data/BBS_site_dist.csv", header=TRUE)

## remove site with only 1 species
BBS_site_dist <- BBS_site_dist[!(BBS_site_dist$GRIDREF=="SU0451"),]
FDis_abund <- FDis_abund[!(FDis_abund$GRIDREF=="SU0451"),] 
FDis_PA <- FDis_PA[!(FDis_PA$GRIDREF=="SU0451"),]
## 199 sites now

## create unique code e.g. 6_4_2_0.5 (trait number: 6, trait number combination: 4, proportion combination: 2, prop. continuous traits: 0.5)
FDis_abund$unique_code = paste(FDis_abund$number_comb, FDis_abund$mix_comb, FDis_abund$p_cont, sep="_")
## same for PA data
FDis_PA$unique_code = paste(FDis_PA$number_comb, FDis_PA$mix_comb, FDis_PA$p_cont, sep="_")
## each unique code has 199 rows i.e. 199 sites

## create 100km blocking factor in each data frame to use as random factor
FDis_abund$GRIDREF_100km <- substr(FDis_abund$GRIDREF, 1, 2) ## take first two characters from gridref
FDis_PA$GRIDREF_100km <- substr(FDis_PA$GRIDREF, 1, 2)
mean_abund$GRIDREF_100km <- substr(mean_abund$GRIDREF, 1, 2)
mean_abund <- merge(mean_abund, BBS_site_dist, by="GRIDREF") ## so mean_abund has lat lon values for morans i test
## create distance matrix for morans i test
dists = as.matrix(dist(mean_abund[,c(4,5)])) ## using lat and lon values
dists_inv = 1/dists
diag(dists_inv) = 0

## check normality of response variable
hist(mean_abund$mean_abund) ## right skew
mean_abund$log_mean_abund <- log(mean_abund$mean_abund) ## log transform looks best
hist(mean_abund$log_mean_abund)

################################
######### Abundance ############
################################

## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

## loops to calculate mixed effects model mean_abund ~ FDis for each group of trait combination
start_time <- Sys.time()
results_table_final <- foreach (i=unique(FDis_abund$unique_code),  .combine=rbind, .packages=c('lme4', 'dplyr', 'purrr', 'MuMIn', 'ape')) %dopar% {
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
end_time - start_time 

## save file
write.csv(results_table_final, file="../Results/Frugivores/R2_abund_seed_det2_gow_lmer.csv", row.names=FALSE)

#######################################
######### Presence/Absence ############
#######################################

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

## loops to calculate linear model mean_abund ~ FDis for each group of trait combination
start_time <- Sys.time()
results_table_final2 <- foreach (i=unique(FDis_abund$unique_code),  .combine=rbind, .packages=c('lme4', 'dplyr', 'purrr', 'MuMIn', 'ape')) %dopar% {
  print(i)
  FDis_temp <- FDis_PA[FDis_PA$unique_code==i,]
  pa_fdis <- list(mean_abund, FDis_temp) %>% reduce(left_join, by = c("GRIDREF", "GRIDREF_100km")) ## join fdis and mean abund data together
  mod2 <- lmer(log_mean_abund ~ FDis + (1|GRIDREF_100km), data=pa_fdis) ## run model
  ## create dataframe with model residuals and site lat & lon values
  df = data.frame(res = residuals(mod2), lat = pa_fdis$latitude, lon = pa_fdis$longitude)
  x <- Moran.I(df$res, dists_inv) ## calculate Moran's I
  ## save R2 and Moran's I p-value (to check for spatial autocorrelation)
  data.frame(R2=r.squaredGLMM(mod2)[1], moran_p_value=x$p.value,
             trait_combination=i)
}
end_time <- Sys.time()
end_time - start_time 
write.csv(results_table_final2, file="../Results/Frugivores/R2_PA_seed_det2_gow_lmer.csv", row.names=FALSE)




############# check lmer moran results - does spatial autocorrelation disappear when the 100km blocking factor is added?

## Frugivore abundance
moran_results_abund_lmer <- read.csv("../Results/Frugivores/R2_abund_seed_det2_gow_lmer.csv",header=TRUE)
nrow(moran_results_abund_lmer[moran_results_abund_lmer$moran_p_value<0.05, ])/nrow(moran_results_abund_lmer)*100 
## only 0.07% results show spatial autocorrelation

## Frugivore P/A 
moran_results_pa_lmer <- read.csv("../Results/Frugivores/R2_PA_seed_det2_gow_lmer.csv",header=TRUE)
nrow(moran_results_pa_lmer[moran_results_pa_lmer$moran_p_value<0.05, ])/nrow(moran_results_pa_lmer)*100 
## still 67% results show spatial autocorrelation

## Insectivore abundance
moran_results_abund_lmer <- read.csv("../Results/Insectivores/R2_abund_seed_det2_gow_lmer.csv",header=TRUE)
nrow(moran_results_abund_lmer[moran_results_abund_lmer$moran_p_value<0.05, ])/nrow(moran_results_abund_lmer)*100 
## only 0.07% results show spatial autocorrelation

## Insectivore P/A
moran_results_pa_lmer <- read.csv("../Results/Insectivores/R2_PA_seed_det2_gow_lmer.csv",header=TRUE)
nrow(moran_results_pa_lmer[moran_results_pa_lmer$moran_p_value<0.05, ])/nrow(moran_results_pa_lmer)*100 
## only 0.07% results show spatial autocorrelation

## looks like frugivore presence-absence results might need a 50km blocking factor instead?
























## old code
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
write.csv(moran_results_final, file="../Results/Frugivores/moran_results_abund_gow_lm.csv", row.names=FALSE)

#######################################
######### Presence/Absence ############
#######################################
moran_results_final2 <- NULL
## loops to calculate linear model mean_abund ~ FDis for each group of trait combination
for(i in unique(FDis_PA$unique_code)){
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
write.csv(moran_results_final2, file="../Results/Frugivores/moran_results_PA_gow_lm.csv", row.names=FALSE)




























# ## add abundance_weighted column
# FDis_abund$abundance_weighted <- "Y"
# FDis_PA$abundance_weighted <- "N"
# 
# ## remove site with only 1 species
# FDis_abund <- FDis_abund[!(FDis_abund$GRIDREF=="SU0451"),] ## 199 sites now
# FDis_abund <- droplevels(FDis_abund)
# FDis_PA <- FDis_PA[!(FDis_PA$GRIDREF=="SU0451"),] ## 199 sites now
# FDis_PA <- droplevels(FDis_PA)
# 
# FDis_summary <- rbind(FDis_abund, FDis_PA) ## add the two FDis datasets together
# FDis_summary <- merge(FDis_summary, mean_abund, by="GRIDREF") ## add in mean abundance
# ## change/remove some columns
# # FDis_summary <- FDis_summary[,-6] ## remove PC axes variation
# names(FDis_summary)[4] <- "prop_cont_comb"
# names(FDis_summary)[5] <- "prop_cont"
# ## merge in trait data combination info
# FDis_summary2 <- merge(FDis_summary, trait_combinations, by=c("prop_cont", "number_comb", "prop_cont_comb"))
# ## create trait number and trait number combination columns
# FDis_summary2$trait_number <- sub("\\_.*", "", FDis_summary2$number_comb) ## keep first values
# FDis_summary2$trait_number_comb <- sub(".*_|^[^_]*$", "", FDis_summary2$number_comb) ## keep last values
# FDis_summary2 <- FDis_summary2[,-2]
# ## create unique code column
# FDis_summary2$unique_code = paste(FDis_summary2$trait_number, FDis_summary2$trait_number_comb,
#                                   FDis_summary2$prop_cont, FDis_summary2$prop_cont_comb, sep="_")
# ## reorder columns
# FDis_summary2 <- FDis_summary2[c(3,6,4,5,16,17,1,2,18,7:15)]
# ## save dataset
# write.csv(FDis_summary2, file="../Results/Frugivores/FDis_summary_gow.csv", row.names=FALSE)
# 
# 
# ## subset trait number = 9 and prop = 1 OR 0.8888889
# FDis_abund_subset1 <- FDis_abund[FDis_abund$number_comb=="9_1" & FDis_abund$p_cont=="1",]
# FDis_abund_subset1 <- FDis_abund_subset1[,c(1,2)]
# names(FDis_abund_subset1)[2] <- "FDis_9_1"
# 
# FDis_abund_subset2 <- FDis_abund[FDis_abund$p_cont=="0.888888888888889",]
# FDis_abund_subset2 <- FDis_abund_subset2[,c(1,2)]
# names(FDis_abund_subset2)[2] <- "FDis_9_0.89"
# 
# FDis_abund_subset <- merge(FDis_abund_subset1, FDis_abund_subset2, by="GRIDREF")
# ggplot(FDis_abund_subset, aes(FDis_9_1, FDis_9_0.89)) +
#   geom_point() +
#   theme_classic()
# 
# 
# FDis_abund_subset2 <- FDis_abund_subset2 %>% group_by(GRIDREF) %>% summarize(FDis_9_0.89 = mean(FDis_9_0.89, na.rm=TRUE))
# FDis_abund_subset <- merge(FDis_abund_subset1, FDis_abund_subset2, by="GRIDREF")
# ggplot(FDis_abund_subset, aes(FDis_9_1, FDis_9_0.89)) +
#   geom_point() +
#   theme_classic()

# mean_abund <- merge(mean_abund_old, mean_abund_new, by="GRIDREF")
# plot(mean_abund$mean_abund.x, mean_abund$mean_abund.y)
# mean_abund$abund_diff <- mean_abund$mean_abund.x - mean_abund$mean_abund.y
# 
# FDis_abund_old <- FDis_abund_old[,-c(3:6)]
# FDis_abund_new <- FDis_abund_new[,-c(3:6)]
# 
# FDis_abund <- merge(FDis_abund_old, FDis_abund_new, by=c("GRIDREF", "unique_code"))
# FDis_abund$FDis_diff <- FDis_abund$FDis.x - FDis_abund$FDis.y
# plot(FDis_abund$FDis.x, FDis_abund$FDis.y)
# 
# FDis_PA_old <- FDis_PA_old[,-c(3:6)]
# FDis_PA_new <- FDis_PA_new[,-c(3:6)]
# 
# FDis_PA <- merge(FDis_PA_old, FDis_PA_new, by=c("GRIDREF", "unique_code"))
# FDis_PA$FDis_diff <- FDis_PA$FDis.x - FDis_PA$FDis.y