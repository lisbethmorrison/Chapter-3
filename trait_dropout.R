###################################################################################
## Title: Difference in R2 when traits are removed for frugivores
## User: Lisbeth Hordley 
## email: l.hordley@pgr.reading.ac.uk
## Date: September 2020
###################################################################################

options(scipen=999)
rm(list = ls())

library(dplyr)
library(tidyr)
library(reshape2)
library(purrr)
library(stringr)
library(tidytext)
library(reshape2)
library(ggplot2)

## read in data
trait_data <- read.csv("../Data/BBS_invert_trait_data.csv", header=TRUE) 
## can use seed or invert dataset - it's only the trait names used here

## prep trait data
trait_data <- trait_data[,c(1,6:13,15)] ## keep effect and both traits (9 total)
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
write.csv(trait_combinations, file="../Data/trait_combinations.csv", row.names=FALSE)

## only interested in trait number 8 & 9
trait_combinations <- trait_combinations[trait_combinations$trait_number==8 | trait_combinations$trait_number==9,]

x <- colnames(trait_data) ## list of traits 1:9
## create column missing_trait which lists all traits which are missing for each trait combination
results2 <- trait_combinations %>% group_by(number_combination) %>% 
  summarise(missing_trait = toString(setdiff(x, trait))) 
## assign each trait 0/1
## 1 = trait missing
## 0 = trait used
dtm <- results2 %>%
  unnest_tokens("word", "missing_trait", token = "regex", pattern = ",") %>% 
  mutate(word = str_trim(word)) %>%
  count(number_combination, word) %>% 
  pivot_wider(names_from = "word", values_from = "n", values_fill = list(n = 0))
## add in combination 9_1 which has no missing traits
add <- c("9_1",0,0,0,0,0,0,0,0,0)
dtm <- rbind(dtm, add) ## bind this to dtm

## find rows with ONE non-matching difference (i.e. where there is a 0 and 1 for ONE trait)
row.names(dtm) <- dtm$number_combination ## move trait combintaion to rownames
dtm[1] <- NULL
dtm[,1:9] <- sapply(dtm[,1:9],as.numeric) ## make 0 and 1 numeric
## calculate absolute sum of the number of times a difference between rows is found
final <- as.data.frame(setNames(combn(1:nrow(dtm), 2, FUN = function(i) sum(abs(dtm[i[1],] - dtm[i[2],]))), 
                                combn(rownames(dtm), 2, toString))) 
## tidy up dataframe
final$comb_match <- row.names(final) ## move comb_match (the matching trait combinations) to a column
row.names(final) <- 1:nrow(final)
colnames(final)[1] <- "difference"
## filter to only look at difference = 1 (i.e. where there is a mismatch for ONE trait)
final<-final[(final$difference==1),] ## 2295 rows
## split up comb_match column
final$trait_comb1 <- gsub(",.*$", "", final$comb_match)
final$trait_comb2 <- sub('.*,\\s*', '', final$comb_match)
## remove first two columns (not needed)
final <- final[,-c(1,2)]

## merge final with results2 which has lists of missing traits (use this to find which trait is unique between matches)
final2 <- merge(final, results2, by.x="trait_comb1", by.y="number_combination", all=TRUE)
colnames(final2)[3] <- "missing_trait1"
final3 <- merge(final2, results2, by.x="trait_comb2", by.y="number_combination", all.x=TRUE)
colnames(final3)[4] <- "missing_trait2"
final3 <- na.omit(final3) ## 2295 rows still

## match up traits between two columns and return non-matching trait
## this is the unique trait between matches - the one that influences the R2
final4 <- final3 %>%
  mutate(unique_trait = map2_chr(strsplit(as.character(missing_trait1), ",\\s+"), 
                                 strsplit(as.character(missing_trait2), ",\\s+"), 
                                 ~ str_c(setdiff(.x, .y), collapse=", ")))
## re-order columns
final4 <- final4[,c(2,1,3,4,5)] 

## save file 
write.csv(final4, file="../Data/trait_comb_matches.csv", row.names=FALSE)

#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
rm(list = ls())

###################### FRUGIVORES ######################
## read in data
model_results_abund <- read.csv("../Results/Frugivores/R2_abund_seed_det2_noPCA_lmer.csv", header=TRUE)
model_results_PA <- read.csv("../Results/Frugivores/R2_PA_seed_det2_noPCA_lmer.csv", header=TRUE)
trait_comb_matches <- read.csv("../Data/trait_comb_matches.csv", header=TRUE)

## tidy up trait_comb_matches
trait_comb_matches$trait_comb1 = paste(trait_comb_matches$trait_comb1, "_1_1", sep="")
trait_comb_matches$trait_comb2 = paste(trait_comb_matches$trait_comb2, "_1_1", sep="")
trait_comb_matches <- trait_comb_matches[-c(4,5)]
colnames(trait_comb_matches)[3] <- "missing_trait"

## merge the two together
## abundance 
trait_matches_abund <- merge(trait_comb_matches, model_results_abund, by.x="trait_comb1", by.y="trait_combination", all.x=TRUE)
colnames(trait_matches_abund)[4] <- "R2_trait_comb1"
trait_matches_abund <- merge(trait_matches_abund, model_results_abund, by.x="trait_comb2", by.y="trait_combination", all.x=TRUE)
colnames(trait_matches_abund)[5] <- "R2_trait_comb2"
## presence-absence
trait_matches_pa <- merge(trait_comb_matches, model_results_PA, by.x="trait_comb1", by.y="trait_combination", all.x=TRUE)
colnames(trait_matches_pa)[4] <- "R2_trait_comb1"
trait_matches_pa <- merge(trait_matches_pa, model_results_PA, by.x="trait_comb2", by.y="trait_combination", all.x=TRUE)
colnames(trait_matches_pa)[5] <- "R2_trait_comb2"

## calculate R2 difference
trait_matches_abund$R2_difference <- trait_matches_abund$R2_trait_comb1 - trait_matches_abund$R2_trait_comb2
trait_matches_pa$R2_difference <- trait_matches_pa$R2_trait_comb1 - trait_matches_pa$R2_trait_comb2

write.csv(trait_matches_abund, file="../Results/Frugivores/trait_match_R2_abund_det2_noPCA_lmer.csv", row.names=FALSE)
write.csv(trait_matches_pa, file="../Results/Frugivores/trait_match_R2_pa_det2_noPCA_lmer.csv", row.names=FALSE)

###################### INSECTIVORES ######################

## read in data
model_results_abund <- read.csv("../Results/Insectivores/R2_abund_det2_noPCA_lmer.csv", header=TRUE)
model_results_PA <- read.csv("../Results/Insectivores/R2_PA_det2_noPCA_lmer.csv", header=TRUE)
trait_comb_matches <- read.csv("../Data/trait_comb_matches.csv", header=TRUE)

## tidy up trait_comb_matches
trait_comb_matches$trait_comb1 = paste(trait_comb_matches$trait_comb1, "_1_1", sep="")
trait_comb_matches$trait_comb2 = paste(trait_comb_matches$trait_comb2, "_1_1", sep="")
trait_comb_matches <- trait_comb_matches[-c(4,5)]
colnames(trait_comb_matches)[3] <- "missing_trait"

## merge the two together
## abundance 
trait_matches_abund <- merge(trait_comb_matches, model_results_abund, by.x="trait_comb1", by.y="trait_combination", all.x=TRUE)
colnames(trait_matches_abund)[4] <- "R2_trait_comb1"
trait_matches_abund <- merge(trait_matches_abund, model_results_abund, by.x="trait_comb2", by.y="trait_combination", all.x=TRUE)
colnames(trait_matches_abund)[5] <- "R2_trait_comb2"
## presence-absence
trait_matches_pa <- merge(trait_comb_matches, model_results_PA, by.x="trait_comb1", by.y="trait_combination", all.x=TRUE)
colnames(trait_matches_pa)[4] <- "R2_trait_comb1"
trait_matches_pa <- merge(trait_matches_pa, model_results_PA, by.x="trait_comb2", by.y="trait_combination", all.x=TRUE)
colnames(trait_matches_pa)[5] <- "R2_trait_comb2"

## calculate R2 difference
trait_matches_abund$R2_difference <- trait_matches_abund$R2_trait_comb1 - trait_matches_abund$R2_trait_comb2
trait_matches_pa$R2_difference <- trait_matches_pa$R2_trait_comb1 - trait_matches_pa$R2_trait_comb2

write.csv(trait_matches_abund, file="../Results/Insectivores/trait_match_R2_abund_det2_noPCA_lmer.csv", row.names=FALSE)
write.csv(trait_matches_pa, file="../Results/Insectivores/trait_match_R2_pa_det2_noPCA_lmer.csv", row.names=FALSE)


########################################################################################################

## Calculate average correlation coefficient between unique trait and rest of traits 
## for each trait combination
rm(list = ls())

trait_combinations <- read.csv(file="../Data/trait_combinations.csv", header=TRUE)
trait_matches <- read.csv("../Data/trait_comb_matches.csv", header=TRUE)

## only interested in trait number 8 & 9
trait_combinations <- trait_combinations[trait_combinations$trait_number==8 | trait_combinations$trait_number==9,]

################# FRUGIVORES ###############
######### calculate average correlation for each unique trait
trait_data <- read.csv("../Data/BBS_seed_trait_data.csv", header=TRUE) 

trait_data <- trait_data[-10,] ## remove crane
rownames(trait_data) <- trait_data$ENGLISH_NAME
trait_data <- trait_data[, -1] ## move species names to row names

av_cor_frug <- NULL
av_cor <- NULL
for(i in unique(trait_combinations$number_combination)){
  print(i)
  comb_temp <- trait_combinations[trait_combinations$number_combination==i,] ## only look at traits of interest (i)
  trait_temp <- trait_data[colnames(trait_data) %in% comb_temp$trait] ## subset main trait data file to traits of interest
  cor_matrix <- cor(trait_temp)
  cor_matrix_long <- melt(cor_matrix)
  ## remove correlations = 1 (correlations between the same trait)
  cor_matrix_long <- cor_matrix_long[!(cor_matrix_long$value==1),]
  
  unique_trait <- trait_matches$unique_trait[trait_matches$trait_comb2==i]
  for (j in unique_trait){
    temp_cor <- cor_matrix_long[cor_matrix_long$Var2==j,]
    ## take absolute average correlation (so e.g. -0.9 will still count as high correlation)
    average_cor <- mean(abs(temp_cor$value))
    av_cor_temp <- data.frame(trait_comb2=i, unique_trait=j, av_cor=average_cor)
    av_cor_frug <- rbind(av_cor_temp, av_cor_frug)
  }
}
write.csv(av_cor_frug, file="../Results/Frugivores/Average_correlation.csv", row.names=FALSE)

################# INSECTIVORES ###############
######### calculate average correlation for each unique trait
trait_data <- read.csv("../Data/BBS_invert_trait_data.csv", header=TRUE) 

rownames(trait_data) <- trait_data$ENGLISH_NAME
trait_data <- trait_data[, -1] ## move species names to row names

av_cor_insect <- NULL
av_cor <- NULL
for(i in unique(trait_combinations$number_combination)){
  print(i)
  comb_temp <- trait_combinations[trait_combinations$number_combination==i,] ## only look at traits of interest (i)
  trait_temp <- trait_data[colnames(trait_data) %in% comb_temp$trait] ## subset main trait data file to traits of interest
  cor_matrix <- cor(trait_temp)
  cor_matrix_long <- melt(cor_matrix)
  ## remove correlations = 1 (correlations between the same trait)
  cor_matrix_long <- cor_matrix_long[!(cor_matrix_long$value==1),]
  
  unique_trait <- trait_matches$unique_trait[trait_matches$trait_comb2==i]
  for (j in unique_trait){
    temp_cor <- cor_matrix_long[cor_matrix_long$Var2==j,]
    ## take absolute average correlation (so e.g. -0.9 will still count as high correlation)
    average_cor <- mean(abs(temp_cor$value))
    av_cor_temp <- data.frame(trait_comb2=i, unique_trait=j, av_cor=average_cor)
    av_cor_insect <- rbind(av_cor_temp, av_cor_insect)
  }
}
write.csv(av_cor_insect, file="../Results/Insectivores/Average_correlation.csv", row.names=FALSE)


###################################################################################################
### MODELS AND GRAPHS 
rm(list = ls())

########### FRUGIVORES
## read in data
trait_matches_abund <- read.csv("../Results/Frugivores/trait_match_R2_abund_det2_lmer.csv", header=TRUE)
trait_matches_PA <- read.csv("../Results/Frugivores/trait_match_R2_PA_det2_noPCA_lmer.csv", header=TRUE)
av_cor_frug <- read.csv("../Results/Frugivores/Average_correlation.csv", header=TRUE)

av_cor_frug$trait_comb2 = paste(av_cor_frug$trait_comb2, "_1_1", sep="")
colnames(av_cor_frug)[2] <- "missing_trait"

## merge correlation scores
trait_matches_abund <- merge(trait_matches_abund, av_cor_frug, by=c("trait_comb2", "missing_trait"), all=TRUE)
trait_matches_PA <- merge(trait_matches_PA, av_cor_frug, by=c("trait_comb2", "missing_trait"), all=TRUE)

## plot difference in R2 for each trait
lookup <- c("Bill_Depth" = "Bill depth", "Bill_TotalCulmen" = "Bill length",
            "Bill_Width" = "Bill width", "Body_Mass" = "Body mass", "Gape_width" = "Gape width",
            "Hand.Wing.Index..Claramunt.2011." = "Hand-wing index", "Kipp.s_Distance" = 
              "Kipp's distance", "Tarsus_Length" = "Tarsus length", "Wing_Chord" = "Wing length")
trait_matches_abund <- trait_matches_abund %>% mutate(missing_trait=lookup[as.character(missing_trait)])
trait_matches_PA <- trait_matches_PA %>% mutate(missing_trait=lookup[as.character(missing_trait)])

## order traits in decreasing R2_difference
trait_matches_abund <- trait_matches_abund[order(trait_matches_abund$R2_difference, decreasing=TRUE),]
trait_matches_abund$missing_trait <- as.character(trait_matches_abund$missing_trait)
trait_matches_abund$missing_trait <- factor(trait_matches_abund$missing_trait, levels=unique(trait_matches_abund$missing_trait))
trait_matches_PA <- trait_matches_PA[order(trait_matches_PA$R2_difference, decreasing=TRUE),]
trait_matches_PA$missing_trait <- as.character(trait_matches_PA$missing_trait)
trait_matches_PA$missing_trait <- factor(trait_matches_PA$missing_trait, levels=unique(trait_matches_PA$missing_trait))

trait_differences_abund_frug <- ggplot(data=trait_matches_abund, aes(x=missing_trait, y=R2_difference)) + 
  geom_point(size=3) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  labs(x="Removed trait", y =expression("Difference in"~R^2))+
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
trait_differences_abund_frug

trait_differences_pa_frug <- ggplot(data=trait_matches_PA, aes(x=missing_trait, y=R2_difference)) + 
  geom_point(size=3) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  labs(x="Removed trait", y =expression("Difference in"~R^2))+
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
trait_differences_pa_frug

## check normality
par(mfrow=c(1,1))
plot(hist(trait_matches_abund$R2_difference))
trait_matches_abund$R2_difference_transf <- log(trait_matches_abund$R2_difference + 10)
plot(hist(trait_matches_abund$R2_difference_transf))

qqnorm(trait_matches_abund$R2_difference_transf)
qqline(trait_matches_abund$R2_difference_transf)

plot(hist(trait_matches_PA$R2_difference))
trait_matches_PA$R2_difference_transf <- log(trait_matches_PA$R2_difference + 10)
plot(hist(trait_matches_PA$R2_difference_transf))
qqnorm(trait_matches_PA$R2_difference_transf)
qqline(trait_matches_PA$R2_difference_transf) 

########## Run models
model <- lm(R2_difference ~ av_cor, data=trait_matches_abund)
summary(model) ## non-significant
par(mfrow=c(2,2))
plot(model)

model2 <- lm(R2_difference ~ av_cor, data=trait_matches_PA)
summary(model2) ## non-significant
par(mfrow=c(2,2))
plot(model2)

### non-parametric regression as R2 difference is non-normal
library(mblm)
abund_seed <- mblm(R2_difference ~ av_cor, data=trait_matches_abund)
pa_seed <- mblm(R2_difference ~ av_cor, data=trait_matches_PA)

########### INSECTIVORES
## read in data
trait_matches_abund <- read.csv("../Results/Insectivores/trait_match_R2_abund_det2_noPCA_lmer.csv", header=TRUE)
trait_matches_PA <- read.csv("../Results/Insectivores/trait_match_R2_PA_det2_noPCA_lmer.csv", header=TRUE)
av_cor_frug <- read.csv("../Results/Insectivores/Average_correlation.csv", header=TRUE)

av_cor_frug$trait_comb2 = paste(av_cor_frug$trait_comb2, "_1_1", sep="")
colnames(av_cor_frug)[2] <- "missing_trait"

## merge correlation scores
trait_matches_abund <- merge(trait_matches_abund, av_cor_frug, by=c("trait_comb2", "missing_trait"), all=TRUE)
trait_matches_PA <- merge(trait_matches_PA, av_cor_frug, by=c("trait_comb2", "missing_trait"), all=TRUE)
## subset data
trait_matches_abund <- trait_matches_abund[trait_matches_abund$trait_comb2 == "9_1_1_1", ]
trait_matches_PA <- trait_matches_PA[trait_matches_PA$trait_comb2 == "9_1_1_1", ]

## plot difference in R2 for each trait
lookup <- c("Bill_Depth" = "Bill depth", "Bill_TotalCulmen" = "Bill length",
            "Bill_Width" = "Bill width", "Body_Mass" = "Body mass", "Gape_width" = "Gape width",
            "Hand.Wing.Index..Claramunt.2011." = "Hand-wing index", "Kipp.s_Distance" = 
              "Kipp's distance", "Tarsus_Length" = "Tarsus length", "Wing_Chord" = "Wing length")
trait_matches_abund <- trait_matches_abund %>% mutate(missing_trait=lookup[as.character(missing_trait)])
trait_matches_PA <- trait_matches_PA %>% mutate(missing_trait=lookup[as.character(missing_trait)])
## order traits in decreasing R2_difference
trait_matches_abund <- trait_matches_abund[order(trait_matches_abund$R2_difference, decreasing=TRUE),]
trait_matches_abund$missing_trait <- as.character(trait_matches_abund$missing_trait)
trait_matches_abund$missing_trait <- factor(trait_matches_abund$missing_trait, levels=unique(trait_matches_abund$missing_trait))
trait_matches_PA <- trait_matches_PA[order(trait_matches_PA$R2_difference, decreasing=TRUE),]
trait_matches_PA$missing_trait <- as.character(trait_matches_PA$missing_trait)
trait_matches_PA$missing_trait <- factor(trait_matches_PA$missing_trait, levels=unique(trait_matches_PA$missing_trait))

trait_differences_abund_insect <- ggplot(data=trait_matches_abund, aes(x=missing_trait, y=R2_difference)) + 
  geom_point(size=3) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  labs(x="Removed trait", y =expression("Difference in"~R^2))+
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
trait_differences_abund_insect

trait_differences_pa_insect <- ggplot(data=trait_matches_PA, aes(x=missing_trait, y=R2_difference)) + 
  geom_point(size=3) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  labs(x="Removed trait", y =expression("Difference in"~R^2))+
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
trait_differences_pa_insect

### put graphs together and save
library(cowplot)
trait_differences <- plot_grid(trait_differences_abund_frug, trait_differences_pa_frug, trait_differences_abund_insect, trait_differences_pa_insect, 
                               labels=c("(a)", "(b)", "(c)", "(d)"), ncol = 2, nrow = 2, hjust=0)
trait_differences
ggsave(file="../Graphs/Final/FigureS7.png", trait_differences, height=12, width=12) 

## check normality
par(mfrow=c(1,1))
plot(hist(trait_matches_abund$R2_difference))
trait_matches_abund$R2_difference_transf <- (trait_matches_abund$R2_difference)^1/3
plot(hist(trait_matches_abund$R2_difference_transf))
qqnorm(trait_matches_abund$R2_difference)
qqline(trait_matches_abund$R2_difference) 
## can't transform = use non-parametric regression/correlation

plot(hist(trait_matches_PA$R2_difference)) 
qqnorm(trait_matches_PA$R2_difference)
qqline(trait_matches_PA$R2_difference)
## can't transform = use non-parametric regression/correlation

########## Run models
library(mblm)
abund_insect <- lm(R2_difference ~ av_cor, data=trait_matches_abund)
summary(abund_insect) ## just non-significant (p=0.055)
par(mfrow=c(2,2))
plot(abund_insect)
cor1 <- cor.test(trait_matches_abund$R2_difference, trait_matches_abund$av_cor,  method="pearson")
cor1 ## non-significant
cor2 <- cor.test(trait_matches_abund$R2_difference, trait_matches_abund$av_cor,  method="spearman")
cor2 ## non-significant

PA_insect <- lm(R2_difference ~ av_cor,  data=trait_matches_PA)
par(mfrow=c(2,2))
plot(PA_insect)
summary(PA_insect) ## non-significant
cor1 <- cor.test(trait_matches_PA$R2_difference, trait_matches_PA$av_cor,  method="pearson")
cor1 ## non-significant
cor2 <- cor.test(trait_matches_PA$R2_difference, trait_matches_PA$av_cor,  method="spearman")
cor2 ## non-significant

