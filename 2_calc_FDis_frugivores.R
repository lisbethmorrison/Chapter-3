###################################################################################
## Title: Calculate FDis for each combination of traits: frugivores
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
bbs_data <- read.csv("../Data/BBS_2004_2018_seed_interpol_repeat_det2.csv", header=TRUE)
trait_data <- read.csv("../Data/BBS_seed_trait_data.csv", header=TRUE)
###

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
## for abundance-weighted and presence/absence separately

## prep trait data
trait_data <- trait_data[-10,] ## remove Crane - no detectability data
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
end_time - start_time ## new = 1.3399 days mins old = ~2.3 days

FDis_final_df <- do.call("rbind", FDis_final) ## convert final list into dataframe 
## save file
write.csv(FDis_final_df, file="../Results/Frugivores/FDis_abund_det2_gow.csv", row.names=FALSE) ## 2-10 traits using abundance data


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
write.csv(FDis_final_df2, file="../Results/Frugivores/FDis_PA_det2_gow.csv", row.names=FALSE) ## 2-10 traits using P/A data

















options(scipen=999)

rm(list = ls())
library(dplyr)
library(tidyr)
library(FD)
library(factoextra)
library(reshape2)

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
trait_data <- trait_data[,c(1,6:15)] ## keep effect and both traits (9 total)
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
## remove trait number of 1 (only investigating trait number 2-9)
trait_combinations<-trait_combinations[!(trait_combinations$trait_number==1),]

#####################################################################################################
######################################### ABUNDANCE #################################################
#####################################################################################################

start_time <- Sys.time()

FDis_cont_final <- NULL
FDis_mix_final <- NULL
FDis_final <- NULL
species <- rownames(trait_data)

## loop through each trait combination
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
    ## need to figure out way to say if all columns are continuous, then do PCA, otherwise just gowdis/FDis calc
    ## also how to get prop of cont value out
    ## number of columns with _cont/total number of columns
    trait_temp2 <- cat_con_comb[cat_con_comb$id==j,] ## only look at traits of interest (j)
    trait_temp3 <- trait_data_final[colnames(trait_data_final) %in% trait_temp2$trait] ## subset main trait data file to traits of interest
    
    
    if(all(sapply(trait_temp3,is.numeric))){ ## if all columns are numeric --> do PCA before calculating FDis
      
      ## centre and scale traits
      # trait_temp3 <- scale(trait_temp3)
      # ## Run PCA on traits
      # pca_traits <- prcomp(trait_temp3) 
      # summary(pca_traits) ## for now select axes which explain min. 85% of variation
      # ## i.e. cumulative proportion>=85
      # cumsum(pca_traits$sdev^2 / sum(pca_traits$sdev^2))
      # s <- summary(prcomp(trait_temp3))
      # cum_prop <- as.data.frame(s$importance[3, ])
      # cum_prop$axes <- row.names(cum_prop)
      # ## if there are 4 or more PCA axes, only look at first 4
      # if(nrow(cum_prop)>=4){
      #   cum_prop <- cum_prop[1:4,]
      # } else { ## else, take all axes
      #   cum_prop <- cum_prop
      # }
      # x <- predict(pca_traits) ## save scores from PCA
      # x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
      # pca_temp <- x[colnames(x) %in% cum_prop$axes] ## subset main PCA axes to axes of interest
      # row.names(pca_temp) <- species
      # ## use these axes as traits
      
      ## calculate gowdis distance matrix for all types of traits
      dist <- gowdis(trait_temp3)
      fdis_temp <- dbFD(dist, BBS_av_abund, calc.FRic = FALSE, calc.FDiv = FALSE)
      x <- as.data.frame(colnames(trait_temp3))
      x2 <- x %>%
        filter(grepl("cont", `colnames(trait_temp3)`))
      y <- length(x2$`colnames(trait_temp3)`)/length(colnames(trait_temp3))
      
      # z2 <- tail(cum_prop, n=1) ## last row of PCA axes
      # z3 <- paste(z2$`s$importance[3, ]`, z2$axes, sep="_") ## prop variance explained _ PC axes number
      
      fdis_temp <- data.frame(GRIDREF=site_match$GRIDREF, FDis=fdis_temp$FDis,number_comb=i, mix_comb=j, p_cont=y)
      FDis_cont_final <- rbind(fdis_temp, FDis_cont_final)
      
    } else { ## if traits are not all continuous, calc FDis on raw data (or do gowdis first??)
      
      ## calculate gowdis distance matrix for all types of traits
      dist <- gowdis(trait_temp3)
      fdis_temp <- dbFD(dist, BBS_av_abund, calc.FRic = FALSE, calc.FDiv = FALSE)
      
      x <- as.data.frame(colnames(trait_temp3))
      x2 <- x %>%
        filter(grepl("cont", `colnames(trait_temp3)`))
      y <- length(x2$`colnames(trait_temp3)`)/length(colnames(trait_temp3))
      
      fdis_temp <- data.frame(GRIDREF=site_match$GRIDREF, FDis=fdis_temp$FDis, number_comb=i, mix_comb=j, p_cont=y)
      FDis_mix_final <- rbind(fdis_temp, FDis_mix_final)
      
    }
    
    
  }
  
  FDis_final <- rbind(FDis_cont_final, FDis_mix_final)
  
}

end_time <- Sys.time()
end_time - start_time

write.csv(FDis_final, file="../Results/Insectivores/FDis_abund_det2_gow.csv", row.names=FALSE) ## 2-9 traits using abundance data


############################################################################################################
######################################### PRESENCE-ABSENCE #################################################
############################################################################################################
## loop to calculate FDis for each trait combination
start_time <- Sys.time()

FDis_cont_final <- NULL
FDis_mix_final <- NULL
FDis_final <- NULL
species <- rownames(trait_data)

## loop through each trait combination
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
    ## need to figure out way to say if all columns are continuous, then do PCA, otherwise just gowdis/FDis calc
    ## also how to get prop of cont value out
    ## number of columns with _cont/total number of columns
    trait_temp2 <- cat_con_comb[cat_con_comb$id==j,] ## only look at traits of interest (j)
    trait_temp3 <- trait_data_final[colnames(trait_data_final) %in% trait_temp2$trait] ## subset main trait data file to traits of interest
    
    
    if(all(sapply(trait_temp3,is.numeric))){ ## if all columns are numeric --> do PCA before calculating FDis
      
      ## centre and scale traits
      # trait_temp3 <- scale(trait_temp3)
      # ## Run PCA on traits
      # pca_traits <- prcomp(trait_temp3) 
      # summary(pca_traits) ## for now select axes which explain min. 85% of variation
      # ## i.e. cumulative proportion>=85
      # cumsum(pca_traits$sdev^2 / sum(pca_traits$sdev^2))
      # s <- summary(prcomp(trait_temp3))
      # cum_prop <- as.data.frame(s$importance[3, ])
      # cum_prop$axes <- row.names(cum_prop)
      # ## if there are 4 or more PCA axes, only look at first 4
      # if(nrow(cum_prop)>=4){
      #   cum_prop <- cum_prop[1:4,]
      # } else { ## else, take all axes
      #   cum_prop <- cum_prop
      # }
      # x <- predict(pca_traits) ## save scores from PCA
      # x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
      # pca_temp <- x[colnames(x) %in% cum_prop$axes] ## subset main PCA axes to axes of interest
      # row.names(pca_temp) <- species
      # ## use these axes as traits
      
      ## calculate gowdis distance matrix for all types of traits
      dist <- gowdis(trait_temp3)
      fdis_temp <- dbFD(dist, BBS_PA, calc.FRic = FALSE, calc.FDiv = FALSE)
      
      x <- as.data.frame(colnames(trait_temp3))
      x2 <- x %>%
        filter(grepl("cont", `colnames(trait_temp3)`))
      y <- length(x2$`colnames(trait_temp3)`)/length(colnames(trait_temp3))
      
      # z2 <- tail(cum_prop, n=1) ## last row of PCA axes
      # z3 <- paste(z2$`s$importance[3, ]`, z2$axes, sep="_") ## prop variance explained _ PC axes number
      
      fdis_temp <- data.frame(GRIDREF=site_match$GRIDREF, FDis=fdis_temp$FDis,number_comb=i, mix_comb=j, p_cont=y)
      FDis_cont_final <- rbind(fdis_temp, FDis_cont_final)
      
    } else { ## if traits are not all continuous, calc FDis on raw data (or do gowdis first??)
      
      ## calculate gowdis distance matrix for all types of traits
      dist <- gowdis(trait_temp3)
      fdis_temp <- dbFD(dist, BBS_PA, calc.FRic = FALSE, calc.FDiv = FALSE)
      
      x <- as.data.frame(colnames(trait_temp3))
      x2 <- x %>%
        filter(grepl("cont", `colnames(trait_temp3)`))
      y <- length(x2$`colnames(trait_temp3)`)/length(colnames(trait_temp3))
      
      fdis_temp <- data.frame(GRIDREF=site_match$GRIDREF, FDis=fdis_temp$FDis, number_comb=i, mix_comb=j, p_cont=y)
      FDis_mix_final <- rbind(fdis_temp, FDis_mix_final)
      
    }
    
    
  }
  
  FDis_final <- rbind(FDis_cont_final, FDis_mix_final)
  
}

write.csv(FDis_final, file="../Results/Insectivores/FDis_PA_det2_gow.csv", row.names=FALSE) ## 2, 6 and 10 traits using abundance data

end_time <- Sys.time()
end_time - start_time




























#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################

### Calculate FDis for each trait separately
## Use for R2 from single-trait FDis in trait-dropout analysis

trait <- colnames(trait_data)
FDis_trait_dropout <- NULL

for(i in trait){print(i)
  trait_temp <- trait_data[!colnames(trait_data) %in% i]
  fdis_temp <- dbFD(trait_temp, BBS_av_abund, calc.FRic = FALSE, calc.FDiv = FALSE)
  fdis_temp <- data.frame(GRIDREF=site_match$GRIDREF, FDis=fdis_temp$FDis, missing_trait=i)
  FDis_trait_dropout <- rbind(fdis_temp, FDis_trait_dropout)
  }
write.csv(FDis_trait_dropout, file="../Results/Frugivores/FDis_trait_dropout_abund_det2.csv", row.names=FALSE)

trait <- colnames(trait_data)
FDis_trait_dropout <- NULL

for(i in trait){print(i)
  trait_temp <- trait_data[!colnames(trait_data) %in% i]
  fdis_temp <- dbFD(trait_temp, BBS_PA, calc.FRic = FALSE, calc.FDiv = FALSE)
  fdis_temp <- data.frame(GRIDREF=site_match$GRIDREF, FDis=fdis_temp$FDis, missing_trait=i)
  FDis_trait_dropout <- rbind(fdis_temp, FDis_trait_dropout)
}
write.csv(FDis_trait_dropout, file="../Results/Frugivores/FDis_trait_dropout_PA_det2.csv", row.names=FALSE)


