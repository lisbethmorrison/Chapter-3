###################################################################################
## Title: Plot model results for insectivores
## User: Lisbeth Hordley 
## email: l.hordley@pgr.reading.ac.uk
## Date: January 2020
###################################################################################

options(scipen=999)
library(ggplot2)
library(dplyr)
library(ggeffects)
library(cowplot)
rm(list = ls())

## read in data
model_results_abund <- read.csv("../Results/Insectivores/R2_abund_det2_gow_lmer.csv", header=TRUE)
model_results_PA <- read.csv("../Results/Insectivores/R2_PA_det2_lmer.csv", header=TRUE)

## create trait number and proportion continuous columns
model_results_abund$trait_number <- sub("\\_.*", "", model_results_abund$trait_combination) ## keep first values
model_results_abund$prop_cont <- sub(".*_|^[^_]*$", "", model_results_abund$trait_combination) ## keep last values
model_results_abund$trait_number <- as.numeric(model_results_abund$trait_number)
model_results_abund$prop_cont <- as.numeric(model_results_abund$prop_cont)
model_results_abund <- model_results_abund[,-2]

model_results_PA$trait_number <- sub("\\_.*", "", model_results_PA$trait_combination) ## keep first values
model_results_PA$prop_cont <- sub(".*_|^[^_]*$", "", model_results_PA$trait_combination) ## keep last values
model_results_PA$trait_number <- as.numeric(model_results_PA$trait_number)
model_results_PA$prop_cont <- as.numeric(model_results_PA$prop_cont)
model_results_PA <- model_results_PA[,-2]

####### PLOT HEATMAP OF R2 AND SD ##########
#### ABUNDANCE #####
## take average and SE of R2 per trait group
mean_r2_abund <- model_results_abund %>% group_by(trait_number, prop_cont) %>% 
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
colnames(mean_r2_abund)[3] <- "R2_mean"

## R2 values
R2_abund <- ggplot(mean_r2_abund, aes(prop_cont, trait_number, fill= R2_mean)) + 
  geom_tile(width=0.1) +
  scale_fill_viridis_c(name =expression(""~R^2)) +
  labs(y="Number of traits", x="Proportion of continuous traits") + 
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
  scale_y_continuous(breaks=seq(from=2,to=10,by=1)) + 
  theme_classic()
R2_abund
ggsave(file="../Graphs/Insectivores_R2_abund_det2.png", R2_abund, height=7, width=10) 

### heatmap of variation (shown as standard deviation)
mean_r2_abund2 <- na.omit(mean_r2_abund) ## remove NAs where variation cannot be calculated (only one combination e.g. all 9 traits)
SD_abund <- ggplot(mean_r2_abund2, aes(prop_cont, trait_number, fill= sd)) + 
  geom_tile(width=0.1) +
  scale_fill_viridis_c(name ="SD") +
  labs(y="Number of traits", x="Proportion of continuous traits") + 
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
  scale_y_continuous(breaks=seq(from=2,to=10,by=1)) +
  theme_classic()
SD_abund
ggsave(file="../Graphs/Insectivores_SD_abund_det2.png", SD_abund, height=7, width=10) 

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
# ggsave(file="../Graphs/Insectivores_slope_abund.png", slope_abund, height=7, width=10) 
# 
# slope_var_abund <- ggplot(mean_r2_abund2, aes(prop_cont, trait_number, fill= slope_sd)) + 
#   geom_tile(width=0.1) +
#   labs(y="Number of traits", x="Proportion of continuous traits") + 
#   scale_fill_viridis_c(name ="SD") +
#   scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
#   scale_y_continuous(breaks=seq(from=2,to=9,by=1)) + 
#   theme_classic()
# slope_var_abund
# ggsave(file="../Graphs/Insectivores_slope_var_abund.png", slope_var_abund, height=7, width=10) 

#### PRESENCE-ABSENCE #####

## take average and SE of R2 per trait group
mean_r2_PA <- model_results_PA %>% group_by(trait_number, prop_cont) %>% 
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
colnames(mean_r2_PA)[3] <- "R2_mean"

## take mean R2 values - shows the correct results
R2_PA <- ggplot(mean_r2_PA, aes(prop_cont, trait_number, fill= R2_mean)) + 
  geom_tile(width=0.1) +
  labs(y="Number of traits", x="Proportion of continuous traits") + 
  scale_fill_viridis_c(name =expression(""~R^2)) +
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
  scale_y_continuous(breaks=seq(from=2,to=9,by=1)) + 
  theme_classic()
R2_PA
ggsave(file="../Graphs/Insectivores_R2_PA_det2.png", R2_PA, height=7, width=10) 

### heatmap of variation (shown as standard deviation)
mean_r2_PA2 <- na.omit(mean_r2_PA)
SD_PA <- ggplot(mean_r2_PA2, aes(prop_cont, trait_number, fill= sd)) + 
  geom_tile(width=0.1) +
  scale_fill_viridis_c(name ="SD") +
  labs(y="Number of traits", x="Proportion of continuous traits") + 
  scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
  scale_y_continuous(breaks=seq(from=2,to=9,by=1)) +
  theme_classic()
SD_PA
ggsave(file="../Graphs/Insectivores_SD_PA_det2.png", SD_PA, height=7, width=10) 

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
# ggsave(file="../Graphs/Insectivores_slope_PA.png", slope_PA, height=7, width=10) 
# 
# slope_var_PA <- ggplot(mean_r2_PA2, aes(prop_cont, trait_number, fill= slope_sd)) + 
#   geom_tile(width=0.1) +
#   labs(y="Number of traits", x="Proportion of continuous traits") + 
#   scale_fill_viridis_c(name ="SD") +
#   scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) + 
#   scale_y_continuous(breaks=seq(from=2,to=9,by=1)) + 
#   theme_classic()
# slope_var_PA
# ggsave(file="../Graphs/Insectivores_slope_var_PA.png", slope_var_PA, height=7, width=10) 

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
plot(hist(mean_r2_abund_cont$R2_mean)) ## positive skew
## reflect the dataset
mean_r2_abund_cont$R2_mean_transf <- log(mean_r2_abund_cont$R2_mean)
plot(hist(mean_r2_abund_cont$R2_mean_transf)) ## positive skew
## can't transform - just leave it..

insect_abund1 <- lm(R2_mean ~ trait_number, data=mean_r2_abund_cont)
summary(insect_abund1) ## non-significant
insect_abund2 <- lm(R2_mean ~ poly(trait_number,2), data=mean_r2_abund_cont)
summary(insect_abund2) ## significant
insect_abund3 <- lm(R2_mean ~ poly(trait_number,3), data=mean_r2_abund_cont)
summary(insect_abund3) ## significant
insect_abund4 <- lm(R2_mean ~ poly(trait_number,4), data=mean_r2_abund_cont)
summary(insect_abund4) ## one term becomes non-significant (use 2nd order - no biological reason to choose 3rd order)
AIC(insect_abund1,insect_abund2,insect_abund3,insect_abund4) ## 3 has lowest AIC
## for noPCA - polynomial model not needed - linear fits best
## significant negative (but v small effect)

par(mfrow=c(2,2))
plot(insect_abund1) ## looks ok

## plot graph
mydf <- ggpredict(insect_abund2, terms = "trait_number")
mydf2 <- merge(mydf, mean_r2_abund_cont, by.x="x", by.y="trait_number")

insect_abund_cont_traits_PCA <- ggplot(mydf2, aes(x = x, y = predicted)) +
  geom_smooth(aes(x = x, y = predicted), color="black", stat="identity") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha=0.1) +
  geom_point(aes(x = x, y = R2_mean), size=2, colour="darkgrey") +
  geom_errorbar(aes(ymin=R2_mean-sd, ymax=R2_mean+sd), width=.2, colour="darkgrey") +
  labs(x="Number of traits", y =expression("Mean"~R^2))+
  theme_classic()
insect_abund_cont_traits_PCA

## Model 2: R2 ~ no. traits * prop_cont of mixed traits
## does the R2 depend on prop continuous and no traits?

## check distribution
plot(hist(mean_r2_abund_mix$R2_mean)) ## positive skew
mean_r2_abund_mix$R2_mean_transf <- sqrt(mean_r2_abund_mix$R2_mean)
plot(hist(mean_r2_abund_mix$R2_mean_transf)) ## looks good
qqnorm(mean_r2_abund_mix$R2_mean_transf)
qqline(mean_r2_abund_mix$R2_mean_transf)

insect_abund2 <- lm(R2_mean ~ trait_number*prop_cont, data=mean_r2_abund_mix)
summary(insect_abund2) ## trait_number + prop_cont significant negative relationship
## non-significant interaction (just, p=0.06)
## for mix/categorical traits, as trait number & prop_cont increases, R2 decreases
par(mfrow=c(2,2))
plot(insect_abund2)
hist(insect_abund2$residuals)
## plot result (trait_number)
r2_predict <- predict(insect_abund2,interval="confidence")
newdata_abund <- cbind(data.frame(mean_r2_abund_mix), data.frame(r2_predict))

newdata_abund$fit <- exp(newdata_abund$fit)
newdata_abund$upr <- exp(newdata_abund$upr)
newdata_abund$lwr <- exp(newdata_abund$lwr)

insect_abund_mix <- ggplot(newdata_abund, aes(x = prop_cont, y = R2_mean, color=as.factor(trait_number))) +
  geom_point(size=2.5) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .15, linetype = 0) +
  labs(x="Proportion of continuous traits", y =expression("Mean"~R^2))+
  labs(color='Number of traits') +
  geom_line(aes(y = fit), size = 1) +
  theme_classic()
insect_abund_mix
## save graph
ggsave(file="../Graphs/Insectivore_meanR2_interaction_abund_det2.png", insect_abund_mix, height=5, width=7) 

## Model 3: mean R2 (all) ~ categorical/continuous variable type

## split into two groups - continuous OR mix/categorical
mean_r2_abund$prop_cont_group <- ifelse(mean_r2_abund$prop_cont==1, "continuous", "mix/categorical")
qqnorm(mean_r2_abund$R2_mean)
qqline(mean_r2_abund$R2_mean)
## Mean R2 is not normally distributed
## use a Mann-Whitney test instead of t test
wilcox.test(R2_mean ~ prop_cont_group, data=mean_r2_abund) 
## significant differene in means 
## plot result (boxplot)
trait_identity <- c("Continuous", "Mix/Categorical")
abund_prop_cont_PCA <- ggplot(mean_r2_abund, aes(x=prop_cont_group, y=R2_mean)) + 
  geom_boxplot() +
  labs(y=expression("Mean"~R^2), x="Variable type") + 
  scale_x_discrete(labels= trait_identity) +
  theme_classic()
abund_prop_cont_PCA
## mean R2 is higher for continuous traits than mix/categorical traits

###################################
##### presence/absence results #####
###################################

## split dataframe into those with continuous trait results and mix/categorical trait results
mean_r2_pa_cont <- mean_r2_PA[mean_r2_PA$prop_cont == "1", ]
mean_r2_pa_mix <- mean_r2_PA[!mean_r2_PA$prop_cont == "1", ]

## Model 1: R2 mean ~ no. continuous traits (polynomial model)
## check normality
plot(hist(mean_r2_pa_cont$R2_mean)) ## negative skew
## can't transform

insect_pa1 <- lm(R2_mean ~ trait_number, data=mean_r2_pa_cont)
summary(insect_pa1) ## non-significant
insect_pa2 <- lm(R2_mean ~ poly(trait_number,2), data=mean_r2_pa_cont)
summary(insect_pa2) ## significant
insect_pa3 <- lm(R2_mean ~ poly(trait_number,3), data=mean_r2_pa_cont)
summary(insect_pa3) ## significant
insect_pa4 <- lm(R2_mean ~ poly(trait_number,4), data=mean_r2_pa_cont)
summary(insect_pa4) ## one term becomes non-significant (use 2nd order - no biological reason to use 3rd order)
AIC(insect_pa1,insect_pa2,insect_pa3,insect_pa4) ## 3 has lowest AIC

## for continuous traits, mean R2 increases as trait number increases
par(mfrow=c(2,2))
plot(insect_pa1) ## looks ok

## plot graph
mydf <- ggpredict(insect_pa2, terms = "trait_number")
mydf2 <- merge(mydf, mean_r2_pa_cont, by.x="x", by.y="trait_number")

insect_pa_cont_traits_PCA <- ggplot(mydf2, aes(x = x, y = predicted)) +
  geom_smooth(aes(x = x, y = predicted), color="black", stat="identity") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha=0.1) +
  geom_point(aes(x = x, y = R2_mean), size=2, colour="darkgrey") +
  geom_errorbar(aes(ymin=R2_mean-sd, ymax=R2_mean+sd), width=.2, colour="darkgrey") +
  labs(x="Number of traits", y =expression("Mean"~R^2))+
  theme_classic()
insect_pa_cont_traits_PCA

## Model 2: R2 ~ no. traits * prop_cont of mixed traits
## does the R2 depend on prop continuous and no traits?

## check distribution
plot(hist(mean_r2_pa_mix$R2_mean)) ## positive skew
mean_r2_pa_mix$R2_mean_transf <- log(mean_r2_pa_mix$R2_mean)
plot(hist(mean_r2_pa_mix$R2_mean_transf)) ## looks ok
qqnorm(mean_r2_pa_mix$R2_mean)
qqline(mean_r2_pa_mix$R2_mean)

## mix of traits interaction
insect_pa2 <- lm(R2_mean ~ trait_number*prop_cont, data=mean_r2_pa_mix)
summary(insect_pa2) ## trait_number non-significnt
## prop_cont significant positive
## interaction is significant
par(mfrow=c(2,2))
plot(insect_pa2)
hist(insect_pa2$residuals)
## plot result
r2_predict <- predict(insect_pa2,interval="confidence")
newdata_pa <- cbind(data.frame(mean_r2_pa_mix), data.frame(r2_predict))

insect_pa_mix <- ggplot(newdata_pa, aes(x = prop_cont, y = R2_mean, color=as.factor(trait_number))) +
  geom_point(size=2.5) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .15, linetype = 0) +
  labs(x="Proportion of continuous traits", y =expression("Mean"~R^2))+
  labs(color='Number of traits') +
  geom_line(aes(y = fit), size = 1) +
  theme_classic()
insect_pa_mix
ggsave(file="../Graphs/Insectivore_meanR2_interaction_pa_det2.png", insect_pa_mix, height=5, width=7) 

## Model 3: mean R2 (all) ~ categorical/continuous variable type

## split into two groups - continuous OR mix/categorical
mean_r2_PA$prop_cont_group <- ifelse(mean_r2_PA$prop_cont==1, "continuous", "mix/categorical")

## check distribution 
qqnorm(mean_r2_PA$R2_mean)
qqline(mean_r2_PA$R2_mean)

## Mean R2 is not normally distributed
## use a Mann-Whitney test instead of t test
wilcox.test(R2_mean ~ prop_cont_group, data=mean_r2_PA) 
## significant differene in means 
## plot result (boxplot)
trait_identity <- c("Continuous", "Mix/Categorical")
pa_prop_cont_PCA <- ggplot(mean_r2_PA, aes(x=prop_cont_group, y=R2_mean)) + 
  geom_boxplot() +
  labs(y=expression("Mean"~R^2), x="Variable type") + 
  scale_x_discrete(labels= trait_identity) +
  theme_classic()
pa_prop_cont_PCA
## mean R2 is higher for continuous traits than mix/categorical traits

########################################################################

### put heatmap graphs together
library(cowplot)
results <- plot_grid(R2_abund, SD_abund, R2_PA, SD_PA, 
                     labels=c("(a)", "(b)", "(c)", "(d)"), ncol = 2, nrow = 2, hjust=0)
results
ggsave(file="../Graphs/Final/Insectivores_FigureS2.png", results, height=7, width=10) 

## plot boxplots together
library(cowplot)
results2 <- plot_grid(insect_abund_cont_traits_noPCA, insect_pa_cont_traits_noPCA, 
                      insect_abund_cont_traits_PCA, insect_pa_cont_traits_PCA,
                      labels=c("(a)", "(b)", "(c)", "(d)"), ncol = 2, nrow = 2, hjust=0)
results2
ggsave(file="../Graphs/Final/Insectivores_cont_traits_FigureS6.png", results2, height=7, width=10)

results3 <- plot_grid(abund_prop_cont_noPCA, pa_prop_cont_noPCA, 
                      abund_prop_cont_PCA, pa_prop_cont_PCA,
                      labels=c("(a)", "(b)", "(c)", "(d)"), ncol = 2, nrow = 2, hjust=0)
results3
ggsave(file="../Graphs/Final/Insectivores_variabletype_FigureS4.png", results3, height=7, width=10)

## put graphs together (interaction results)
trait_number_mix <- plot_grid(frug_abund_mix, insect_pa_mix, 
                              ncol=2, nrow = 1, 
                              labels=c("(a)", "(b)"),hjust=0)
trait_number_mix
ggsave(file="../Graphs/Final/Figure3.png.png", trait_number_mix, height=4, width=10) 





