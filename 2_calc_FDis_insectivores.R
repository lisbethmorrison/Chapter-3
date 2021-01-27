###################################################################################
## Title: Calculate FDis for each combination of traits: insectivores
## User: Lisbeth Hordley 
## email: l.hordley@pgr.reading.ac.uk
## Date: January 2020
###################################################################################
options(scipen=999)

rm(list = ls())
library(dplyr)
library(tidyr)
library(FD)
library(factoextra)
library(reshape2)
library(purrr)
library(stringr)
library(tidytext)

## read in data
bbs_data <- read.csv("../Data/BBS_2004_2018_invert_interpol_repeat_det2.csv", header=TRUE)
trait_data <- read.csv("../Data/BBS_invert_trait_data.csv", header=TRUE)

## Step 1:
## calculate mean relative abundance at each site
## and presence/absence at each site 
## to weight FDis metric

## put BBS data into a dataframe with species as columns and sites as rows filled with average abundance across iterations and year
BBS_av_abund <- bbs_data %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(mean_abund=mean(TOT_det)) %>% 
  spread(ENGLISH_NAME, mean_abund, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_av_abund[,(1)]
## add in numbers 1:200
site_match$site <- 1:200
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_av_abund <- BBS_av_abund[,-(1)]

## calculate P/A at each BBS site
BBS_PA <- bbs_data %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(PA=mean(TOT_det))
BBS_PA$PA <- 1
BBS_PA <- BBS_PA %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  spread(ENGLISH_NAME, PA, drop = FALSE, fill = 0) 
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_PA <- BBS_PA[,-(1)]

## Step 2:
## calculate FDis for each combination of traits

## prep trait data
trait_data <- trait_data[,c(1,6:15)] ## keep effect and both traits (10 total)
rownames(trait_data) <- trait_data[, 1]
trait_data <- trait_data[, -1] ## move species names to row names

## function to calculate all combinations of trait numbers to loop through
make_combinations <- function(x) {
  
  l <- length(x)
  mylist <- lapply(1:l, function(y) {
    combn(x, y, simplify = FALSE)
  })
  mylist
  
}

results <- make_combinations(colnames(trait_data)) ## use above function to calculate all combinations of traits for each number of traits
trait_combinations <- melt(results) ## convert list into dataframe
## reorder columns
trait_combinations <- trait_combinations[,c(1,3,2)]
## rename columns
colnames(trait_combinations) <- c("trait", "trait_number", "trait_combination")
## new column with unique combination per number of trait
trait_combinations$number_combination = paste(trait_combinations$trait_number, trait_combinations$trait_combination, sep="_")
## remove trait number of 1 (only investigating trait number 2-10)
trait_combinations<-trait_combinations[!(trait_combinations$trait_number==1),]

#####################################################################################################
######################################### ABUNDANCE #################################################
#####################################################################################################

## create vector for FDis values to fill
FDis_final<- vector(mode = "list", length = 11805600) # preallocate size

## loop through each trait combination
start_time <- Sys.time()
count<-1
for(i in unique(trait_combinations$number_combination)){
  print(i)
  trait_temp <- trait_combinations[trait_combinations$number_combination==i,] ## only look at traits of interest (i)
  trait_cont <- trait_data[colnames(trait_data) %in% trait_temp$trait] ## subset main trait data file to traits of interest
  
  trait_cat <- as.data.frame(lapply(trait_cont, cut, breaks=2, labels=c("low", "high"))) ## create categorical data
  colnames(trait_cont) <- paste(colnames(trait_cont), "cont", sep = "_") ## add continuous suffix
  colnames(trait_cat) <- paste(colnames(trait_cat), "cat", sep = "_") ## add categorical suffix 
  trait_data_final <- cbind(trait_cont, trait_cat) ## dataframe with all continuous and categorical traits, used later
  
  results2 <- do.call(expand.grid,
                      data.frame(rbind(names(trait_cont),names(trait_cat)))) ## create all combinations of cont/cat traits
  cat_con_comb <- reshape(results2, direction="long", 
                          varying=list(colnames(results2)), v.names = c("trait", "id"), timevar=NULL) ## change format
  cat_con_comb <-cat_con_comb[order(cat_con_comb$id),] ## order by ID i.e. the groups of cat/con combinations
  
  ## loop through each ID
  for(j in unique(cat_con_comb$id)){
    trait_temp2 <- cat_con_comb[cat_con_comb$id==j,] ## only look at traits of interest (j)
    trait_temp3 <- trait_data_final[colnames(trait_data_final) %in% trait_temp2$trait] ## subset main trait data file to traits of interest
    
    ## calculate gower distance matrix 
    dist <- gowdis(trait_temp3)
    fdis_temp <- dbFD(dist, BBS_av_abund, calc.FRic = FALSE, calc.FDiv = FALSE, 
                      calc.FGR = FALSE, calc.CWM = FALSE) ## calculate FDis only
    ## calculate p_cont - the porportion of traits that are continuous
    x <- as.data.frame(colnames(trait_temp3))
    x2 <- x %>%
      filter(grepl("cont", `colnames(trait_temp3)`))
    y <- length(x2$`colnames(trait_temp3)`)/length(colnames(trait_temp3)) ## calculate p_cont
    ## save all results into a dataframe
    fdis_temp <- data.frame(GRIDREF=site_match$GRIDREF, FDis=fdis_temp$FDis, number_comb=i, mix_comb=j, p_cont=y)
    ## add results to final list
    FDis_final[[count]] <- fdis_temp
    count <- count + 1
    
  }
}
end_time <- Sys.time()
end_time - start_time 

FDis_final_df <- do.call("rbind", FDis_final) ## convert final list into dataframe 
## save file
write.csv(FDis_final_df, file="../Results/Insectivores/FDis_abund_det2_gow.csv", row.names=FALSE) ## 2-10 traits using abundance data


############################################################################################################
######################################### PRESENCE-ABSENCE #################################################
############################################################################################################

## create vector for FDis values to fill
FDis_final2<- vector(mode = "list", length = 11805600) # preallocate size

## loop through each trait combination
start_time <- Sys.time()
count<-1
for(i in unique(trait_combinations$number_combination)){
  print(i)
  trait_temp <- trait_combinations[trait_combinations$number_combination==i,] ## only look at traits of interest (i)
  trait_cont <- trait_data[colnames(trait_data) %in% trait_temp$trait] ## subset main trait data file to traits of interest
  
  trait_cat <- as.data.frame(lapply(trait_cont, cut, breaks=2, labels=c("low", "high"))) ## create categorical data
  colnames(trait_cont) <- paste(colnames(trait_cont), "cont", sep = "_") ## add cont suffix
  colnames(trait_cat) <- paste(colnames(trait_cat), "cat", sep = "_") ## add cat suffix
  trait_data_final <- cbind(trait_cont, trait_cat) ## dataframe with all continuous and categorical traits, used later
  
  results2 <- do.call(expand.grid,
                      data.frame(rbind(names(trait_cont),names(trait_cat)))) ## create all combinations of cont/cat traits
  cat_con_comb <- reshape(results2, direction="long", 
                          varying=list(colnames(results2)), v.names = c("trait", "id"), timevar=NULL) ## change format
  cat_con_comb <-cat_con_comb[order(cat_con_comb$id),] ## order by ID i.e. the groups of trait combinations
  
  ## loop through each ID
  for(j in unique(cat_con_comb$id)){
    trait_temp2 <- cat_con_comb[cat_con_comb$id==j,] ## only look at traits of interest (j)
    trait_temp3 <- trait_data_final[colnames(trait_data_final) %in% trait_temp2$trait] ## subset main trait data file to traits of interest
    
    ## calculate gowdis distance matrix 
    dist <- gowdis(trait_temp3)
    fdis_temp <- dbFD(dist, BBS_PA, calc.FRic = FALSE, calc.FDiv = FALSE, 
                      calc.FGR = FALSE, calc.CWM = FALSE) ## calculate FDis only)
    ## calculate p_cont - the porportion of traits that are continuous
    x <- as.data.frame(colnames(trait_temp3))
    x2 <- x %>%
      filter(grepl("cont", `colnames(trait_temp3)`))
    y <- length(x2$`colnames(trait_temp3)`)/length(colnames(trait_temp3))
    ## save all results into a dataframe
    fdis_temp <- data.frame(GRIDREF=site_match$GRIDREF, FDis=fdis_temp$FDis,number_comb=i, mix_comb=j, p_cont=y)
    ## add results to final list
    FDis_final2[[count]] <- fdis_temp
    count <- count + 1      
  }
}
end_time <- Sys.time()
end_time - start_time 

FDis_final_df2 <- do.call("rbind", FDis_final2) ## convert final list into dataframe 
## save file
write.csv(FDis_final_df2, file="../Results/Insectivores/FDis_PA_det2_gow.csv", row.names=FALSE) ## 2-10 traits using P/A data

