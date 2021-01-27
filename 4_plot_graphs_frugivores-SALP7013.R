###################################################################################
## Title: Plot model results for frugivores
## User: Lisbeth Hordley 
## email: l.hordley@pgr.reading.ac.uk
## Date: January 2020
###################################################################################

options(scipen=999)
library(ggplot2)
library(dplyr)
library(ggeffects)

rm(list = ls())

## read in data
model_results_abund <- read.csv("../Results/Frugivores/R2_abund_seed_det2.csv", header=TRUE)
model_results_PA <- read.csv("../Results/Frugivores/R2_PA_seed_det2.csv", header=TRUE)

## create trait number and proportion continuous columns
model_results_abund$trait_number <- sub("\\_.*", "", model_results_abund$trait_combination) ## keep first values
model_results_abund$prop_cont <- sub(".*_|^[^_]*$", "", model_results_abund$trait_combination) ## keep last values
model_results_abund$trait_number <- as.numeric(model_results_abund$trait_number)
model_results_abund$prop_cont <- as.numeric(model_results_abund$prop_cont)
model_results_abund <- model_results_abund[,-3]

model_results_PA$trait_number <- sub("\\_.*", "", model_results_PA$trait_combination) ## keep first values
model_results_PA$prop_cont <- sub(".*_|^[^_]*$", "", model_results_PA$trait_combination) ## keep last values
model_results_PA$trait_number <- as.numeric(model_results_PA$trait_number)
model_results_PA$prop_cont <- as.numeric(model_results_PA$prop_cont)
model_results_PA <- model_results_PA[,-3]


####### PLOT HEATMAP OF R2 AND SD ##########
#### ABUNDANCE #####
## take average and SE of R2 per trait group
mean_r2_abund <- model_results_abund %>% group_by(trait_number, prop_cont) %>% 
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))

## R2 values
R2_abund <- ggplot(mean_r2_abund, aes(prop_cont, trait_number, fill= R2_mean)) + 
  geom_tile(width=0.1) +
  labs(y="Number of traits", x="Proportion of continuous traits") + 
  scale_fill_viridis_c(name =expression(""~R^2)) +
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
  scale_y_continuous(breaks=seq(from=2,to=9,by=1)) + 
  theme_classic()
R2_abund
ggsave(file="../Graphs/Frugivores_R2_abund_det2.png", R2_abund, height=7, width=10) 

### Standard Deviation 
mean_r2_abund2 <- na.omit(mean_r2_abund) ## remove NAs where variation cannot be calculated (only one combination e.g. all 9 traits)
SD_abund <- ggplot(mean_r2_abund2, aes(prop_cont, trait_number, fill= R2_sd)) + 
  geom_tile(width=0.1) +
  scale_fill_viridis_c(name ="SD") +
  labs(y="Number of traits", x="Proportion of continuous traits") +
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
  scale_y_continuous(breaks=seq(from=2,to=9,by=1)) +
  theme_classic()
SD_abund
ggsave(file="../Graphs/Frugivores_SD_abund_det2.png", SD_abund, height=7, width=10) 

# ## Model slope + variation plots (for supp. material)
# ## plot heatmap of slope values (leaves gaps)
# slope_abund <- ggplot(mean_r2_abund, aes(prop_cont, trait_number, fill= slope_mean)) + 
#   geom_tile(width=0.1) +
#   labs(y="Number of traits", x="Proportion of continuous traits") + 
#   scale_fill_viridis_c(name ="Slope") +
#   scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
#   scale_y_continuous(breaks=seq(from=2,to=9,by=1)) + 
#   theme_classic()
# slope_abund
# ggsave(file="../Graphs/Frugivores_slope_abund.png", slope_abund, height=7, width=10) 
# 
# slope_var_abund <- ggplot(mean_r2_abund2, aes(prop_cont, trait_number, fill= slope_sd)) + 
#   geom_tile(width=0.1) +
#   labs(y="Number of traits", x="Proportion of continuous traits") + 
#   scale_fill_viridis_c(name ="SD") +
#   scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
#   scale_y_continuous(breaks=seq(from=2,to=9,by=1)) + 
#   theme_classic()
# slope_var_abund
# ggsave(file="../Graphs/Frugivores_slope_var_abund.png", slope_var_abund, height=7, width=10) 

#### PRESENCE-ABSENCE #####

## take average and SE of R2 per trait group
mean_r2_PA <- model_results_PA %>% group_by(trait_number, prop_cont) %>% 
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))

## R2 values
R2_PA <- ggplot(mean_r2_PA, aes(prop_cont, trait_number, fill= R2_mean)) + 
  geom_tile(width=0.1) +
  labs(y="Number of traits", x="Proportion of continuous traits") + 
  scale_fill_viridis_c(name =expression(""~R^2)) +
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
  scale_y_continuous(breaks=seq(from=2,to=9,by=1)) + 
  theme_classic()
R2_PA
ggsave(file="../Graphs/Frugivores_R2_PA_det2.png", R2_PA, height=7, width=10) 

## Standard Deviation
mean_r2_PA2 <- na.omit(mean_r2_PA)
SD_PA <- ggplot(mean_r2_PA2, aes(prop_cont, trait_number, fill= R2_sd)) + 
  geom_tile(width=0.1) +  
  scale_fill_viridis_c(name ="SD") +
  labs(y="Number of traits", x="Proportion of continuous traits") +
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
  scale_y_continuous(breaks=seq(from=2,to=9,by=1)) +
  theme_classic()
SD_PA
ggsave(file="../Graphs/Frugivores_SD_PA_det2.png", SD_PA, height=7, width=10) 

# ## Model slope + variation plots (for supp. material)
# ## plot heatmap of slope values (leaves gaps)
# slope_PA <- ggplot(mean_r2_PA, aes(prop_cont, trait_number, fill= slope_mean)) + 
#   geom_tile(width=0.1) +
#   labs(y="Number of traits", x="Proportion of continuous traits") + 
#   scale_fill_viridis_c(name ="Slope") +
#   scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
#   scale_y_continuous(breaks=seq(from=2,to=9,by=1)) + 
#   theme_classic()
# slope_PA
# ggsave(file="../Graphs/Frugivores_slope_PA.png", slope_PA, height=7, width=10) 
# 
# slope_var_PA <- ggplot(mean_r2_PA2, aes(prop_cont, trait_number, fill= slope_sd)) + 
#   geom_tile(width=0.1) +
#   labs(y="Number of traits", x="Proportion of continuous traits") + 
#   scale_fill_viridis_c(name ="SD") +
#   scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
#   scale_y_continuous(breaks=seq(from=2,to=9,by=1)) + 
#   theme_classic()
# slope_var_PA
# ggsave(file="../Graphs/Frugivores_slope_var_PA.png", slope_var_PA, height=7, width=10) 

####################################################################################################
####################################################################################################

## More tests to explain variation in R2 values

############################
##### abundance models ##### 
############################

## split dataframe into those with continuous trait results and mix/categorical trait results
mean_r2_abund_cont <- mean_r2_abund[mean_r2_abund$prop_cont == "1", ]
mean_r2_abund_mix <- mean_r2_abund[!mean_r2_abund$prop_cont == "1", ]

## Model 1: R2 mean ~ no. continuous traits (polynomial model)
## check normality
plot(hist(mean_r2_abund_cont$R2_mean)) ## looks good
qqnorm(mean_r2_abund_cont$R2_mean)
qqline(mean_r2_abund_cont$R2_mean)

frug_abund1 <- lm(R2_mean ~ trait_number, data=mean_r2_abund_cont)
summary(frug_abund1) ## non-significant
frug_abund2 <- lm(R2_mean ~ poly(trait_number,2), data=mean_r2_abund_cont)
summary(frug_abund2) ## significant - use this one
frug_abund3 <- lm(R2_mean ~ poly(trait_number,3), data=mean_r2_abund_cont)
summary(frug_abund3) ## non-significant
frug_abund4 <- lm(R2_mean ~ poly(trait_number,4), data=mean_r2_abund_cont)
summary(frug_abund4) ## non-significant
AIC(frug_abund1,frug_abund2,frug_abund3,frug_abund4) ## 4 has lowest AIC, but 2 makes most sense
## but actually linear and squared v similar

par(mfrow=c(2,2))
plot(frug_abund2) ## looks ok

## plot graph
mydf <- ggpredict(frug_abund2, terms = "trait_number")
mydf2 <- merge(mydf, mean_r2_abund_cont, by.x="x", by.y="trait_number")

frug_abund_cont_traits <- ggplot(mydf2, aes(x = x, y = predicted)) +
  geom_smooth(aes(x = x, y = predicted), color="black", stat="identity") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha=0.1) +
  geom_point(aes(x = x, y = R2_mean), size=2, colour="darkgrey") +
  geom_errorbar(aes(ymin=R2_mean-R2_sd, ymax=R2_mean+R2_sd), width=.2, colour="darkgrey") +
  labs(x="Number of traits", y =expression("Mean"~R^2))+
  theme_classic()
frug_abund_cont_traits 

## Model 2: R2 ~ no. traits * prop_cont of mixed traits
## does the R2 depend on prop continuous and no traits?

## check distribution
plot(hist(mean_r2_abund_mix$R2_mean)) ## very slight negative skew
qqnorm(mean_r2_abund_mix$R2_mean) 
qqline(mean_r2_abund_mix$R2_mean) ## looks good

## run model
frug_abund2 <- lm(R2_mean ~ trait_number*prop_cont, data=mean_r2_abund_mix)
summary(frug_abund2) ## trait number and prop_cont significant positive relationship
## interaction non-significant
## for mix/categorical traits, as trait number & prop_cont increases, R2 increases
par(mfrow=c(2,2))
plot(frug_abund2)
## plot result
r2_predict <- predict(frug_abund2,interval="confidence")
newdata_abund <- cbind(data.frame(mean_r2_abund_mix), data.frame(r2_predict))

frug_abund_mix <- ggplot(newdata_abund, aes(x = prop_cont, y = R2_mean, color=as.factor(trait_number))) +
  geom_point(size=2.5) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .15, linetype = 0) +
  labs(x="Proportion of continuous traits", y =expression("Mean"~R^2))+
  labs(color='Number of traits') +
  geom_line(aes(y = fit), size = 1) +
  theme_classic()
frug_abund_mix
ggsave(file="../Graphs/Frugivore_meanR2_mix_traits_abund_det2.png", frug_abund_mix, height=5, width=7) 

## Model 3: mean R2 (all) ~ categorical/continuous variable type

## split into two groups - continuous OR mix/categorical
mean_r2_abund$prop_cont_group <- ifelse(mean_r2_abund$prop_cont==1, "continuous", "mix/categorical")
plot(hist(mean_r2_abund$R2_mean)) ## negative skew
mean_r2_abund$R2_mean_sqr <- (mean_r2_abund$R2_mean)^2
plot(hist(mean_r2_abund$R2_mean_sqr))
qqnorm(mean_r2_abund$R2_mean_sqr) 
qqline(mean_r2_abund$R2_mean_sqr) ## looks good

## log mean R2 is normally distributed - t test
t.test(R2_mean_sqr ~ prop_cont_group, data=mean_r2_abund) 
## significant differene in means 
## plot result (boxplot)
library(ggpubr)
trait_identity <- c("Continuous", "Mix/Categorical")
abund_prop_cont <- ggplot(mean_r2_abund, aes(x=prop_cont_group, y=R2_mean_sqr)) + 
  geom_boxplot() +
  labs(y=expression("Mean"~R^2), x="Variable type") + 
  scale_x_discrete(labels= trait_identity) +
  theme_classic()
abund_prop_cont
## mean R2 is higher for mix/categorical traits than continuous traits


###################################
##### presence/absence models ##### 
###################################

## split dataframe into those with continuous trait results and mix/categorical trait results
mean_r2_pa_cont <- mean_r2_PA[mean_r2_PA$prop_cont == "1", ]
mean_r2_pa_mix <- mean_r2_PA[!mean_r2_PA$prop_cont == "1", ]

## Model 1: R2 mean ~ no. continuous traits (polynomial model)
## check normality

## check distribution
plot(hist(mean_r2_pa_cont$R2_mean)) ## left skew
mean_r2_pa_cont$R2_mean_transf <- log(mean_r2_pa_cont$R2_mean)
plot(hist(mean_r2_pa_cont$R2_mean_transf)) ## looks good
qqnorm(mean_r2_pa_cont$R2_mean_transf)
qqline(mean_r2_pa_cont$R2_mean_transf)

frug_pa1 <- lm(R2_mean ~ trait_number, data=mean_r2_pa_cont)
summary(frug_pa1) ## non-significant
frug_pa2 <- lm(R2_mean ~ poly(trait_number,2), data=mean_r2_pa_cont)
summary(frug_pa2) ## significant
frug_pa3 <- lm(R2_mean ~ poly(trait_number,3), data=mean_r2_pa_cont)
summary(frug_pa3) ## significant
frug_pa4 <- lm(R2_mean ~ poly(trait_number,4), data=mean_r2_pa_cont)
summary(frug_pa4) ## one term becomes non-significant (use 2nd order - don't have biological reason for 3rd order)
AIC(frug_pa1,frug_pa2,frug_pa3,frug_pa4) ## 4 has lowest AIC

par(mfrow=c(2,2))
plot(frug_pa2) ## looks ok

## plot graph
mydf <- ggpredict(frug_pa2, terms = "trait_number")
mydf2 <- merge(mydf, mean_r2_pa_cont, by.x="x", by.y="trait_number")

frug_pa_cont_traits <- ggplot(mydf2, aes(x = x, y = predicted)) +
  geom_smooth(aes(x = x, y = predicted), color="black", stat="identity") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha=0.1) +
  geom_point(aes(x = x, y = R2_mean), size=2, colour="darkgrey") +
  geom_errorbar(aes(ymin=R2_mean-R2_sd, ymax=R2_mean+R2_sd), width=.2, colour="darkgrey") +
  labs(x="Number of traits", y =expression("Mean"~R^2))+
  theme_classic()
frug_pa_cont_traits

## Model 2: R2 ~ no. traits * prop_cont of mixed traits
## does the R2 depend on prop continuous and no traits?

## check distribution
plot(hist(mean_r2_pa_mix$R2_mean)) ## positive skew
mean_r2_pa_mix$R2_mean_transf <- log(mean_r2_pa_mix$R2_mean)
plot(hist(mean_r2_pa_mix$R2_mean_transf)) ## looks good
qqnorm(mean_r2_pa_mix$R2_mean_transf)
qqline(mean_r2_pa_mix$R2_mean_transf)

## run model
frug_pa2 <- lm(R2_mean_transf ~ trait_number*prop_cont, data=mean_r2_pa_mix)
summary(frug_pa2) ## trait_number negative, prop_cont positive 
## interaciton non-significant
## as trait number decreases and prop_cont increases, R2 increases
par(mfrow=c(2,2))
plot(frug_pa2)
## plot result
r2_predict <- predict(frug_pa2,interval="confidence")
newdata_pa <- cbind(data.frame(mean_r2_pa_mix), data.frame(r2_predict))

newdata_pa$fit <- exp(newdata_pa$fit)
newdata_pa$upr <- exp(newdata_pa$upr)
newdata_pa$lwr <- exp(newdata_pa$lwr)

frug_pa_mix <- ggplot(newdata_pa, aes(x = prop_cont, y = R2_mean, color=as.factor(trait_number))) +
  geom_point(size=2.5) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .15, linetype = 0) +
  labs(x="Proportion of continuous traits", y =expression("Mean"~R^2))+
  labs(color='Number of traits') +
  geom_line(aes(y = fit), size = 1) +
  theme_classic()
  #theme(legend.position = "none")
frug_pa_mix
ggsave(file="../Graphs/Frugivore_meanR2_interaction_pa_det2.png", frug_pa_mix, height=5, width=7) 

## Model 3: mean R2 (all) ~ categorical/continuous variable type

## split into two groups - continuous OR mix/categorical
mean_r2_PA$prop_cont_group <- ifelse(mean_r2_PA$prop_cont==1, "continuous", "mix/categorical")

## check distribution 
plot(hist(mean_r2_PA$R2_mean)) ## positive skew
mean_r2_PA$R2_mean_log <- log(mean_r2_PA$R2_mean)
plot(hist(mean_r2_PA$R2_mean_log))
qqnorm(mean_r2_PA$R2_mean_log)
qqline(mean_r2_PA$R2_mean_log)

## log mean R2 is normally distributed - t test
t.test(R2_mean_log ~ prop_cont_group, data=mean_r2_PA) 
## significant differene in means 
## plot result (boxplot)
trait_identity <- c("Continuous", "Mix/Categorical")
pa_prop_cont <- ggplot(mean_r2_PA, aes(x=prop_cont_group, y=R2_mean_log)) + 
  geom_boxplot() +
  labs(y=expression("(log) Mean"~R^2), x="Variable type") + 
  scale_x_discrete(labels= trait_identity) +
  theme_classic()
pa_prop_cont
## mean R2 is higher for continuous traits than mix/categorical traits

### put graphs together and save
library(cowplot)
results <- plot_grid(R2_abund, SD_abund, R2_PA, SD_PA, 
                     labels=c("(a)", "(b)", "(c)", "(d)"), ncol = 2, nrow = 2, hjust=0)
results
ggsave(file="../Graphs/Frugivores_heatmap_det2.png", results, height=7, width=10) 

# results2 <- plot_grid(slope_abund, slope_var_abund, slope_PA, slope_var_PA, 
#                      labels=c("(a)", "(b)", "(c)", "(d)"), ncol = 2, nrow = 2, hjust=0)
# results2
# ggsave(file="../Graphs/Frugivores_slope_results.png", results2, height=7, width=10) 



#############################################################################################################
## look at relationships between FDis and mean_abund for best and worst models (i.e. highest and lowest R2 values)
## Checked with Tom == relationships look linear so models are fine

FDis_abund2 <- read.csv("../Results/Frugivores/FDis_abund.csv", header=TRUE)
FDis_PA <- read.csv("../Results/Frugivores/FDis_PA_det2.csv", header=TRUE)
mean_abund <- read.csv("../Data/Seed_mean_abund_det2.csv", header=TRUE)

mean_abund$log_mean_abund <- sqrt(mean_abund$mean_abund) ## sqrt not log now
mean_abund2$log_mean_abund <- log(mean_abund2$mean_abund) 

## compare 5_93_1_1 between old and new

## new
FDis_abund$unique_code = paste(FDis_abund$number_comb, FDis_abund$mix_comb, FDis_abund$p_cont, sep="_")

new_r2 <- FDis_abund[FDis_abund$unique_code=="6_51_1_1",]
new_r2 <- new_r2[,c(1,2,7)] 
new_r2 <- merge(new_r2, mean_abund, by="GRIDREF")
new_r2_plot <- ggplot(new_r2, aes(x = FDis, y = log_mean_abund)) +
  geom_point(colour="black") +
  stat_smooth(method = "lm", col = "black") +
  labs(y="(log) Mean total community abundance", x=expression("FDis")) +
  ggtitle("New data") +
  theme_classic()
new_r2_plot ## no relationship

mod_new <- lm(log_mean_abund ~ FDis, new_r2)
summary(mod_new) ## significant, v low R2

## old
FDis_abund2$unique_code = paste(FDis_abund2$number_comb, FDis_abund2$mix_comb, FDis_abund2$p_cont, sep="_")

old_r2 <- FDis_abund2[FDis_abund2$unique_code=="6_51_1_1",]
old_r2 <- old_r2[,c(1,2,7)]
old_r2 <- merge(old_r2, mean_abund2, by="GRIDREF")
old_r2_plot <- ggplot(old_r2, aes(x = FDis, y = log_mean_abund)) +
  geom_point(colour="black") +
  stat_smooth(method = "lm", col = "black") +
  labs(y="(log) Mean total community abundance", x=expression("FDis")) +
  ggtitle("Old data") +
  theme_classic()
old_r2_plot ## positive relationship

mod_old <- lm(log_mean_abund ~ FDis, old_r2)
summary(mod_old) ## significant, R2 ~0.04

#### abundance 
## highest R2 = 0.10597867 for 4_26_2_0.75
## lowest R2 = 0.000000009491479 for 4_98_1_1
FDis_abund$unique_code = paste(FDis_abund$number_comb, FDis_abund$mix_comb, FDis_abund$p_cont, sep="_")

highest_r2 <- FDis_abund[FDis_abund$unique_code=="4_26_2_0.75",]
highest_r2 <- highest_r2[,c(1,2,7)]
highest_r2 <- merge(highest_r2, mean_abund, by="GRIDREF")
highest_r2_plot <- ggplot(highest_r2, aes(x = FDis, y = log_mean_abund)) +
  geom_point(colour="black") +
  stat_smooth(method = "lm", col = "black") +
  labs(y="(sqrt) Mean total community abundance", x=expression("FDis")) +
  ggtitle("Trait number = 4 \n Proportion of continuous traits = 0.75 \n R^2 = 0.12") +
  theme_classic()
highest_r2_plot ## negative relationship

mod <- lm(log_mean_abund ~ FDis, highest_r2)
summary(mod) ## significant

lowest_r2 <- FDis_abund[FDis_abund$unique_code=="4_98_1_1",]
lowest_r2 <- lowest_r2[,c(1,2,7)]
lowest_r2 <- merge(lowest_r2, mean_abund, by="GRIDREF")
lowest_r2_plot <- ggplot(lowest_r2, aes(x = FDis, y = log_mean_abund)) +
  geom_point(colour="black") +
  stat_smooth(method = "lm", col = "black") +
  labs(y="(sqrt) Mean total community abundance", x=expression("FDis")) +
  ggtitle("Trait number = 4 \n Proportion of continuous traits = 1 \n R^2 = 0.0000000095") +
  theme_classic()
lowest_r2_plot

mod2 <- lm(log_mean_abund ~ FDis, lowest_r2)
summary(mod2) ## non-significant

## check 9 continuous traits
nine_cont <- FDis_abund[FDis_abund$unique_code=="9_1_1_1",]
nine_cont <- nine_cont[,c(1,2,7)]
nine_cont <- merge(nine_cont, mean_abund, by="GRIDREF")

mod3 <- lm(log_mean_abund ~ FDis, nine_cont)
summary(mod3) ## non-significant
## plot back-transformed results
r2_predict <- predict(mod3,interval="confidence")
newdata <- cbind(data.frame(nine_cont), data.frame(r2_predict))

newdata$fit <- (newdata$fit)^2
newdata$upr <- (newdata$upr)^2
newdata$lwr <- (newdata$lwr)^2

nine_cont_plot <- ggplot(newdata, aes(x = FDis, y = mean_abund)) +
  geom_point(size=2.5) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .15, linetype = 0) +
  labs(x="FDis", y ="Mean abundance")+
  geom_line(aes(y = fit), size = 1) +
  theme_classic()
nine_cont_plot ## no relationship

## compare to 9 traits, 0.9 prop cont
nine_cont2 <- FDis_abund2[FDis_abund2$unique_code=="9_1_12_0.666666666666667",]
nine_cont2 <- nine_cont2[,c(1,2,7)]
nine_cont2 <- merge(nine_cont2, mean_abund2, by="GRIDREF")

mod4 <- lm(log_mean_abund ~ FDis, nine_cont2)
summary(mod4) ## non-significant
## plot back-transformed results
r2_predict <- predict(mod4,interval="confidence")
newdata2 <- cbind(data.frame(nine_cont2), data.frame(r2_predict))

newdata2$fit <- exp(newdata2$fit)
newdata2$upr <- exp(newdata2$upr)
newdata2$lwr <- exp(newdata2$lwr)

nine_cont_plot2 <- ggplot(newdata2, aes(x = FDis, y = mean_abund)) +
  geom_point(size=2.5) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .15, linetype = 0) +
  labs(x="FDis", y ="Mean abundance")+
  geom_line(aes(y = fit), size = 1) +
  theme_classic()
nine_cont_plot2 ## negative relationship

#### presence/absence 
## highest R2 = 0.15043289 for 2_14_4_0
## lowest R2 = 0.0000000005748083 for 6_72_59_0.333333333333333
FDis_PA$unique_code = paste(FDis_PA$number_comb, FDis_PA$mix_comb, FDis_PA$p_cont, sep="_")

highest_r2 <- FDis_PA[FDis_PA$unique_code=="5_81_1_1",]
highest_r2 <- highest_r2[,c(1,2,7)]
highest_r2 <- merge(highest_r2, mean_abund, by="GRIDREF")
highest_r2_plot <- ggplot(highest_r2, aes(x = FDis, y = log_mean_abund)) +
  geom_point(colour="black") +
  stat_smooth(method = "lm", col = "black") +
  labs(y="(sqrt) Mean total community abundance", x=expression("FDis")) +
  ggtitle("Trait number = 2 \n Proportion of continuous traits = 0 \n R^2 = 0.15") +
  theme_classic()
highest_r2_plot ## positive relationship

mod <- lm(log_mean_abund ~ FDis, highest_r2)
summary(mod) ## significant

lowest_r2 <- FDis_PA[FDis_PA$unique_code=="5_122_2_0.8",]
lowest_r2 <- lowest_r2[,c(1,2,7)]
lowest_r2 <- merge(lowest_r2, mean_abund, by="GRIDREF")
lowest_r2_plot <- ggplot(lowest_r2, aes(x = FDis, y = log_mean_abund)) +
  geom_point(colour="black") +
  stat_smooth(method = "lm", col = "black") +
  labs(y="(sqrt) Mean total community abundance", x=expression("FDis")) +
  ggtitle("Trait number = 6 \n Proportion of continuous traits = 0.33 \n R^2 = 0.00000000057") +
  theme_classic()
lowest_r2_plot ## no relationship

mod2 <- lm(log_mean_abund ~ FDis, lowest_r2)
summary(mod2) ## non-significant
