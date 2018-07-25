# Run Analysis:
# w/ clean dataset: prepared in line with Thewissen & Rueda
# w/ pure dataset: as few additional variables as possible

# Date: Feb 16, 2018
# Author: Jonas Markgraf
#########################

rm(list = ls())

# load libraries
library(mvtnorm)
library(car)
library(gtools)
library(foreign)  # to import dta file
library(readxl)  # to load data
library(repmis)  # to set wd
library(MASS)  # for 'polr' function
library(plyr)  # for list to df
library(dplyr)
library(tidyr)  # to turn data in long format
library(doBy)  # summary statistics
library(DataCombine)  # to calculate percentage change
library(zoo)  # to interpolate gini values
library(lme4)  # to estimate multilevel models
#library(sjPlot)  # to predict probabilities for models
library(effects)  # to plot MLM effects
library(arm)  # to use 'se.ranef' for extracting standard error from hierarchical model
library(stargazer)  # to export regression results
#library(xlsx)  # to export xlsx table
library(Amelia)
library (runjags)    # run Jags from R
library (MCMCpack)  # canned MCMC procedures
library (mcmcplots) # plot MC
library (parallel)
library (readstata13)
library (interplot)  # to plot interaction effects
library (mice)

# set wd
possibles <- c("~/Dropbox/CreditPreferences/")
set_valid_wd(possibles)
graphicsPath <- c("~/Dropbox/CreditPreferences/Data/plotsGraphs/")

# load dataset
ess <- read.dta13("Data/ESSdata/finalESS_ThewissenRueda.dta",  # to be refined by each of us.
  convert.factors = F)
ess.backup <- ess
ess <- ess.backup

dataKey <- data.frame(varLabels = attr(ess, "var.labels"))

# recode variables included in analysis
unique(ess$brwmny)
ess$brwmny <- ifelse(ess$brwmny > 5, NA, ess$brwmny)

unique (ess$rlgdgr)
ess$rlgdgr <- ifelse(ess$rlgdgr > 11, NA, ess$rlgdgr)

unique (ess$agea)
ess$agea <- ifelse(ess$agea > 120, NA, ess$agea)

# recode homeowner variable for robustness
ess$h.owner <- ifelse(ess$hhmodwl >2, NA, ess$hhmodwl)
ess$h.owner <- ifelse(ess$h.owner == 2, 0, 1)

unique(ess$gincdif2)
unique(ess$male)
hist(log(ess$perceqnatincdollar))
unique(ess$our)  # occupational unemployment rate
quantile(ess$our, c(.25,.75), na.rm = T)

# generate key variables
ess$cntry.yr <- paste0(ess$cou, ess$year)

ess$country.year.isco <- paste0 (ess$cou, ess$year, ess$iscoco2)

ess$oesch2a <- ifelse (ess$iscoco>=60000, NA, 
                       ifelse ((ess$iscoco >= 3430 & ess$iscoco <= 3431) |
                          (ess$iscoco >= 4000 & ess$iscoco <= 4195) |
                          (ess$iscoco >= 4210 & ess$iscoco <= 4215) |
                           ess$iscoco == 4223, 1, 0))

ess$oesch2b <- ifelse ((ess$iscoco >= 1 & ess$iscoco <= 100) |
                          (ess$iscoco >= 6100 & ess$iscoco <= 7113) |
                          (ess$iscoco >= 7200 & ess$iscoco <= 8290) |
                          (ess$iscoco >= 9000 & ess$iscoco <= 9001) |
                          (ess$iscoco >= 9150 & ess$iscoco <= 9151) |
                          (ess$iscoco >= 9153 & ess$iscoco <= 9161) |
                          (ess$iscoco >= 9200 & ess$iscoco <= 9311), 1, 
                       ifelse (ess$iscoco>=60000, NA, 0))

ess$Oeschroutine <- ifelse (ess$oesch2a==1 | ess$oesch2b==1, 1, 0)
ess$rti2 <- ifelse ( ess$iscoco2==61 | ess$iscoco2==62 | ess$iscoco2==92, 0.89,
                     ifelse (ess$iscoco2==11 | ess$iscoco2==23 | ess$iscoco2==33, -0.68,
                             ess$rti))
   
# Model to build latent risk variables 
complete.ess <- ess[!is.na(ess$iscoco2),]

# allRisk <- c("rti2", "our", "relskillspec", "Oeschroutine")
# offShoreRisk <- c("offshwalt","offsh")
# omitObs <- as.numeric (attr (na.omit (complete.ess[,c(allRisk, offShoreRisk)])
#                              , "na.action"))

# allRisk.PC <- princomp(complete.ess[-omitObs,allRisk])
# offShoreRisk.PC <- princomp(complete.ess[-omitObs,offShoreRisk])
# round (cor (cbind (allRisk.PC$scores[,1], offShoreRisk.PC$scores[,1], complete.ess[-omitObs,c(allRisk,offShoreRisk)])), 2)

# # Including all variables
# Risk.PC  <- princomp(complete.ess[-omitObs, c(allRisk,offShoreRisk)])
# round (cor (cbind (Risk.PC$scores[,1], Risk.PC$scores[,2], complete.ess[-omitObs,c(allRisk,offShoreRisk)])), 2)

# # Excluding offshwalt
# Risk.PC  <- princomp(complete.ess[-omitObs, c(allRisk,"offsh")])
# round (cor (cbind (Risk.PC$scores[,1], Risk.PC$scores[,2], complete.ess[-omitObs,c(allRisk,offShoreRisk)])), 2)

#### MCMC mixed factor analysis ###
# Turn ordinal variables into ordered factors
complete.ess$Oeschroutine.fac <- as.ordered (complete.ess$Oeschroutine)
complete.ess$offshwalt.fac <- as.ordered (cut(complete.ess$offshwalt, c(-5,20,40,60,80,100)))
levels (complete.ess$offshwalt.fac) <- c("1","2","3","4","5")
complete.ess$offshwalt.fac.num <- as.numeric (complete.ess$offshwalt.fac)

complete.ess <- complete.ess[is.na(complete.ess$wrongincome) | 
                                complete.ess$wrongincome==0,]


#####################################
#### Analysis by income category ####
#####################################

# T&R then use their relative measure of income and compare low income (50 percentile - sd)
quantile(complete.ess$perceqnatincdollar, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm =T) # p50 = 0.83
sd(complete.ess$perceqnatincdollar, na.rm =T) # sd = 0.74

# high income: perceqnatincdollar = 1.57; low income: perceqnatincdollar = 0.09
complete.ess$incomeTER.TR <- ifelse(complete.ess$perceqnatincdollar > 1.57, "High"
                                    , ifelse(complete.ess$perceqnatincdollar < 0.09, "Low", "Middle"))


# generate quintile groups for income by cntr.yr pair
complete.ess$incomeQNT.TR <- ntile(complete.ess$perceqnatincdollar, 5)
complete.ess$incomeQNT.TR <- relevel(as.factor(complete.ess$incomeQNT.TR), ref = 3)
complete.ess$incomeLOG <- log(complete.ess$perceqnatincdollar)

# generate quintile groups for income by cntr.yr pair
complete.ess <- ddply(complete.ess, .(cntry.yr), mutate, incomeQNT = ntile(income, 5))
complete.ess$incomeQNT <- relevel(as.factor(complete.ess$incomeQNT), ref = 3)

# generate tertile groups for income by cntry.yr pair
complete.ess <- ddply(complete.ess, .(cntry.yr), mutate, incomeTER = ntile(income, 3))
complete.ess$incomeTER <- relevel(as.factor(complete.ess$incomeTER), ref = 2)



#### FACTOR ANALYSIS ####
# Reduce dataset even further, to exclude observations with "wrong income"
# The main risk measure that we use is a factor that excludes OUR
AllFactor <- factanal (~rti2+relskillspec+offsh+Oeschroutine+offshwalt.fac.num 
                       , factors=1
                       , rotation="varimax"
                       , scores="regression"
                       , na.action=na.exclude
                       , data=complete.ess)

AllFactor$loadings
factor.scores <- AllFactor$scores

# Notation follows factanal help
# Useful math help from https://stats.stackexchange.com/questions/126885/methods-to-compute-factor-scores-and-what-is-the-score-coefficient-matrix-in
Lambda <- AllFactor$loadings
Phi    <- diag(AllFactor$uniquenesses)
Sigma <- Lambda %*% t(Lambda) + Phi
predictors <- complete.ess[,c("rti2","relskillspec","offsh","Oeschroutine","offshwalt.fac.num")]
std.predictors <- apply (predictors, 2, function(x) (x-mean(x, na.rm=T))/sd(x, na.rm=T))
corrMatrix <- AllFactor$correlation # similar to cor (predictors, use="p")

regression.scores <- std.predictors %*% Lambda  # These are easiest 
thompson.scores   <- std.predictors %*% solve (corrMatrix) %*% Lambda # This is what R does
bartlett.scores   <- std.predictors %*% t(solve (t(Lambda) %*% Phi %*% Lambda) %*% t(Lambda) %*% solve (Phi)) # correlate almost perfectly with thompson

# To fill in factor scores for all observations, even those with some missing values,
# we use as a predictor matrix a patched std.predictor matrix where NAs are substituted
# by column means (i.e., 0).

# First, we need to distinguish between units without any observed predictors,
# and units with some observed predictors
num.Missing <- apply (std.predictors, 1, function (x) sum (is.na(x)))
allMissing <- ifelse (num.Missing==5, 1, 0)

## Do not need to run multiple imputation again, if
## file in line 199 is re-loaded
# # Multiple imputation of observations that have at least one (1) "risk" measure
# mi.predictors <- mice (predictors[allMissing==0,], m=1
#                        , method=c("pmm","pmm","pmm","logreg","polr"))
# mi.predictors <- complete (mi.predictors)
# save (mi.predictors, file="~/Dropbox/CreditPreferences/Data/riskPCdata/riskPC.Rda")
load (file="~/Dropbox/CreditPreferences/Data/riskPCdata/riskPC.Rda")

# rownames (mi.predictors) <- rownames(std.predictors)[allMissing==0]

# Standardize multiply-imputed predictors (to construct factanal scores)
std.mi.predictors <- apply (mi.predictors, 2, function(x) (x-mean(x, na.rm=T))/sd(x, na.rm=T))

# Create Thompson scores (same as factanal "regression" option)
thompson.scores.all <- std.mi.predictors %*% solve (corrMatrix) %*% Lambda  
complete.ess$risk <- thompson.scores.all

pdf (paste0 (graphicsPath, "RiskHistogram.pdf"), h=5, w=7)
par (mar=c(4,4,1,1))
hist (complete.ess$risk, main="", xlab="Relative risk of losing job")
dev.off()

#### Create risk/income groups ####
risky.rich <- ifelse (complete.ess$risk > 0 
                      & complete.ess$incomeTER==3
                      & !is.na(complete.ess$incomeTER), 1, 0)
riskless.rich <- ifelse (complete.ess$risk < 0 
                      & complete.ess$incomeTER==3
                      & !is.na(complete.ess$incomeTER), 1, 0)
risky.poor <- ifelse (complete.ess$risk > 0 
                      & complete.ess$incomeTER==1
                      & !is.na(complete.ess$incomeTER), 1, 0)
riskless.poor <- ifelse (complete.ess$risk < 0 
                      & complete.ess$incomeTER==1
                      & !is.na(complete.ess$incomeTER), 1, 0)

#### Distribution of preference for redistribution among risk/income groups ####
par (mfrow=c(2,2), mar=c(3,3,1,1))
hist (complete.ess$gincdif2[risky.rich==1]
      , breaks=seq(0.5,5.5,by=1), ylim=c(0,14000)
      , main="Risky rich")
hist (complete.ess$gincdif2[risky.poor==1]
      , breaks=seq(0.5,5.5,by=1), ylim=c(0,14000)
      , main="Risky poor")
hist (complete.ess$gincdif2[riskless.rich==1]
      , breaks=seq(0.5,5.5,by=1), ylim=c(0,14000)
      , main="Riskless rich")
hist (complete.ess$gincdif2[riskless.poor==1]
      , breaks=seq(0.5,5.5,by=1), ylim=c(0,14000)
      , main="Riskless poor")



par (mfrow=c(2,2), mar=c(3,3,1,1))
hist (complete.ess$gincdif2bin[risky.rich==1]
      , ylim=c(0,30000)
      , main="Risky rich")
hist (complete.ess$gincdif2bin[risky.poor==1]
      , ylim=c(0,30000)
      , main="Risky poor")
hist (complete.ess$gincdif2bin[riskless.rich==1]
      , ylim=c(0,30000)
      , main="Riskless rich")
hist (complete.ess$gincdif2bin[riskless.poor==1]
      , ylim=c(0,30000)
      , main="Riskless poor")

print ("risky rich")
length(complete.ess$gincdif2bin[complete.ess$gincdif2bin==0 & risky.rich==1])/length(complete.ess$gincdif2bin[complete.ess$gincdif2bin==1 & risky.rich==1])
print ("riskless rich")
length(complete.ess$gincdif2bin[complete.ess$gincdif2bin==0 & riskless.rich==1])/length(complete.ess$gincdif2bin[complete.ess$gincdif2bin==1 & riskless.rich==1])
print ("risky poor")
length(complete.ess$gincdif2bin[complete.ess$gincdif2bin==0 & risky.poor==1])/length(complete.ess$gincdif2bin[complete.ess$gincdif2bin==1 & risky.poor==1])
print ("riskless poor")
length(complete.ess$gincdif2bin[complete.ess$gincdif2bin==0 & riskless.poor==1])/length(complete.ess$gincdif2bin[complete.ess$gincdif2bin==1 & riskless.poor==1])


#### Distribution of access to credit among risk/income groups ####
par (mfrow=c(2,2), mar=c(3,3,1,1))
hist (complete.ess$brwmny[risky.rich==1]
      , breaks=seq(0.5,5.5,by=1), ylim=c(0,14000)
      , main="Risky rich")
hist (complete.ess$brwmny[risky.poor==1]
      , breaks=seq(0.5,5.5,by=1), ylim=c(0,14000)
      , main="Risky poor")
hist (complete.ess$brwmny[riskless.rich==1]
      , breaks=seq(0.5,5.5,by=1), ylim=c(0,14000)
      , main="Riskless rich")
hist (complete.ess$brwmny[riskless.poor==1]
      , breaks=seq(0.5,5.5,by=1), ylim=c(0,14000)
      , main="Riskless poor")

#### Create dichotomous brwmny ####
complete.ess$brwmnybin <- ifelse (complete.ess$brwmny <= 3 & !is.na(complete.ess$brwmny), 0
                                  , ifelse (complete.ess$brwmny > 3 & !is.na(complete.ess$brwmny), 1, NA))

print ("risky rich")
mean(complete.ess$brwmnybin[risky.rich==1], na.rm=T)
print ("riskless rich")
mean(complete.ess$brwmnybin[riskless.rich==1], na.rm=T)
print ("risky poor")
mean(complete.ess$brwmnybin[risky.poor==1], na.rm=T)
print ("riskless poor")
mean(complete.ess$brwmnybin[riskless.poor==1], na.rm=T)





#### Means and standard deviations ####
# Risky rich
mn.risky.rich.credit <- mean (complete.ess$gincdif2bin[risky.rich==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)], na.rm=T)
mn.risky.rich.no.credit <- mean (complete.ess$gincdif2bin[risky.rich==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)], na.rm=T)
sd.risky.rich.credit <- sd (complete.ess$gincdif2bin[risky.rich==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)], na.rm=T)
sd.risky.rich.no.credit <- sd (complete.ess$gincdif2bin[risky.rich==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)], na.rm=T)
ln.risky.rich.credit <- length (complete.ess$gincdif2bin[!is.na(complete.ess$gincdif2bin) & risky.rich==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)])
ln.risky.rich.no.credit <- length (complete.ess$gincdif2bin[!is.na(complete.ess$gincdif2bin) & risky.rich==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)])

# Riskless rich
mn.riskless.rich.credit <- mean (complete.ess$gincdif2bin[riskless.rich==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)], na.rm=T)
mn.riskless.rich.no.credit <- mean (complete.ess$gincdif2bin[riskless.rich==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)], na.rm=T)
sd.riskless.rich.credit <- sd (complete.ess$gincdif2bin[riskless.rich==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)], na.rm=T)
sd.riskless.rich.no.credit <- sd (complete.ess$gincdif2bin[riskless.rich==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)], na.rm=T)
ln.riskless.rich.credit <- length (complete.ess$gincdif2bin[!is.na(complete.ess$gincdif2bin) & riskless.rich==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)])
ln.riskless.rich.no.credit <- length (complete.ess$gincdif2bin[!is.na(complete.ess$gincdif2bin) & riskless.rich==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)])

# Risky poor
mn.risky.poor.credit <- mean (complete.ess$gincdif2bin[risky.poor==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)], na.rm=T)
mn.risky.poor.no.credit <- mean (complete.ess$gincdif2bin[risky.poor==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)], na.rm=T)
sd.risky.poor.credit <- sd (complete.ess$gincdif2bin[risky.poor==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)], na.rm=T)
sd.risky.poor.no.credit <- sd (complete.ess$gincdif2bin[risky.poor==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)], na.rm=T)
ln.risky.poor.credit <- length (complete.ess$gincdif2bin[!is.na(complete.ess$gincdif2bin) & risky.poor==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)])
ln.risky.poor.no.credit <- length (complete.ess$gincdif2bin[!is.na(complete.ess$gincdif2bin) & risky.poor==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)])

# Riskless poor
mn.riskless.poor.credit <- mean (complete.ess$gincdif2bin[riskless.poor==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)], na.rm=T)
mn.riskless.poor.no.credit <- mean (complete.ess$gincdif2bin[riskless.poor==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)], na.rm=T)
sd.riskless.poor.credit <- sd (complete.ess$gincdif2bin[riskless.poor==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)], na.rm=T)
sd.riskless.poor.no.credit <- sd (complete.ess$gincdif2bin[riskless.poor==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)], na.rm=T)
ln.riskless.poor.credit <- length (complete.ess$gincdif2bin[!is.na(complete.ess$gincdif2bin) & riskless.poor==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)])
ln.riskless.poor.no.credit <- length (complete.ess$gincdif2bin[!is.na(complete.ess$gincdif2bin) & riskless.poor==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)])

means.cat <- c(mn.risky.rich.credit, mn.risky.rich.no.credit
  , mn.riskless.rich.credit, mn.riskless.rich.no.credit
  , mn.risky.poor.credit, mn.risky.poor.no.credit
  , mn.riskless.poor.credit, mn.riskless.poor.no.credit)

sd.cat <- c(sd.risky.rich.credit, sd.risky.rich.no.credit
               , sd.riskless.rich.credit, sd.riskless.rich.no.credit
               , sd.risky.poor.credit, sd.risky.poor.no.credit
               , sd.riskless.poor.credit, sd.riskless.poor.no.credit)

plot (c(0,5), c(0,1), type="n"
      , axes=F, main="")
mtext (side=1, line=0
       , at=c(1:4), text=c("risky","riskless","risky","riskless"))
mtext (side=1, line=1
       , at=c(1:4), text=c("rich","rich","poor","poor"))
axis (2)
points (y=means.cat
        , x=c(0.8,1.2,1.8,2.2,2.8,3.2,3.8,4.2)
        , col=rep(c("black","grey"),4)
        , pch=19)
segments (x0=c(0.8,1.2,1.8,2.2,2.8,3.2,3.8,4.2)
          , x1=c(0.8,1.2,1.8,2.2,2.8,3.2,3.8,4.2)
          , y0=means.cat+(sd.cat)^2
          , y1=means.cat-(sd.cat)^2
          , col=rep(c("black","grey"),4))

# Risky poor
mean (complete.ess$gincdif2bin[risky.poor==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)], na.rm=T)
mean (complete.ess$gincdif2bin[risky.poor==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)], na.rm=T)

# Risky rich
mean (complete.ess$gincdif2bin[riskless.poor==1 & complete.ess$brwmnybin==1 & !is.na(complete.ess$brwmnybin)], na.rm=T)
mean (complete.ess$gincdif2bin[riskless.poor==1 & complete.ess$brwmnybin==0 & !is.na(complete.ess$brwmnybin)], na.rm=T)

#### Plot differences across risk/income groups and credit/no credit groups ####



########################################
#### Bayesian mixed factor analysis #### 
########################################
# MCMCmixfactanal (~rti2+our #+relskillspec+Oeschroutine+offshwalt.fac+offsh
#                  , factors=2
#                  , lambda.constraints=list(rti2=list(2,"+"),
#                                            # our=list(2,"+")) #,
#                                            # relskillspec=list(2,"+"),
#                                            # Oeschroutine=list(2,"+"),
#                                            # offsh=list(2,"+"),
#                                            offshwalt=list(3,"+"))
#                  , data=complete.ess
#                  , burnin=10
#                  , mcmc=20
#                  , thin=10
#                  , verbose=1
#                  , seed=12345
#                  , lambda.start=NA
#                  , l0=0
#                  , L0=1
#                  , store.scores=F
#                  , store.lambda=F)
# 
# source ("BayesGRMmodel.R")
# Matrix with responses to risk items
# Y <- as.matrix (cbind(tmp$rti2, tmp$relskillspec, tmp$offsh, tmp$Oeschroutine, as.numeric(tmp$offshwalt.fac)))
# complete.missing <- apply (Y, 1, invalid)
# 
# Y <- Y[complete.missing==FALSE,]
# Y[,1:3] <- apply (Y[,1:3], 2, function (x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
# N <- nrow(Y)   # Number of observations
# J <- ncol(Y)   # Number of risk items
# L <- 3 # number of continuous variables
# intercept <- rep(1,N)
# 
# source ("~/Dropbox/CreditPreferences/Code/mixFactAnal.R")
# 
# # Priors
# jags.parameters <- c("lambda","deviance","factor.score","kappa","sigma","alpha")
# Data=list(Y=Y, J=J, N=N, L=L
#           , intercept=intercept)
# jags.data <- dump.format(Data) #
# jags.inits <- function()
# {
#    dump.format(
#       list(
#          lambda=runif(5,0,2)
#          , alpha=rnorm(4)  # To coincide with negative priors
#          , kappa.unsorted=rnorm(3)
#          #, beta=rnorm(1,0,2)
#       ) )
# }
# 
# # Run model (so far, failure to define node lambda)
# jags.model <- run.jags( model=mixFactAnal
#                         , monitor=jags.parameters
#                         , inits=jags.inits()#list(jags.inits(), jags.inits())
#                         , n.chains=1
#                         , method="parallel"
#                         , data=jags.data
#                         , adapt=1#000
#                         , thin=2#0
#                         , burnin=20#000
#                         , sample=50#0
#                         #, thin=1, burnin=10, sample=50
#                         , summarise=FALSE
#                         , plots=FALSE )
###########################################################################

# Some variables need to be rescaled
include.vars <- c("gincdif2bin","gincdif2","brwmny","male"
                  ,"agea","unemplindiv","perceqnatincdollar"
                  ,"eduyrs2","mbtru2","rlgdgr","socgdp","gdpc"
                  ,"gininet","epl","unempl"
                  ,"mainsample","incomeTER.TR","incomeQNT.TR"
                  ,"incomeLOG","incomeQNT","incomeTER"
                  ,"h.owner","socialorigin","essround"
                  ,"risk","income","cntry.yr","dweight","country.year.isco")
tmp <- complete.ess[,include.vars]
tmp$log.income <- log(tmp$income)
tmp$log.gdpc   <- log(tmp$gdpc)
colnames(tmp)[grep("incomeLOG", colnames(tmp))] <- "log.perceqnatincdollar"

vars2rescale <- c("agea","eduyrs2","socgdp"
                  ,"gininet")  #,"log.perceqnatincdollar"
                  #,"log.income","log.gdpc")

# The following function can be included inside lmer to rescale
rescale.var <- function (x) {
   res.x <- (x -mean(x, na.rm=T))/sd(x, na.rm=T)
   return (res.x)}

for (i in 1:length (vars2rescale)) {
   var <- vars2rescale[i]
   resc.var <- rescale.var(tmp[,grep (paste0("^",var,"$"), colnames(tmp))])
   tmp[,grep (paste0("^",var,"$"), colnames(tmp))] <- resc.var
}

#### MODELS ####
# Baseline model
fit.baseline <- lmer(gincdif2 ~ brwmny
  + log.income + male + agea + unemplindiv 
  + eduyrs2 + mbtru2 + rlgdgr 
  + socgdp + log.gdpc
  + (1 + brwmny | cntry.yr),
  weights = dweight, data=tmp)
summary(fit.baseline)

# Meltzer-Richard Fork
fit.gini <- lmer(gincdif2 ~ brwmny*gininet
                     + log.income + male + agea + unemplindiv 
                     + eduyrs2 + mbtru2 + rlgdgr 
                     + socgdp + log.gdpc
                     + (1 + brwmny | cntry.yr),
                     weights = dweight, data=tmp)
summary(fit.gini)
interplot(fit.gini, "brwmny", "gininet")



# Interaction with income
fit.income <- lmer(gincdif2 ~ brwmny + log.income
                   + brwmny:log.income
                   + male + agea + unemplindiv 
                   + eduyrs2 + mbtru2 + rlgdgr 
                   + socgdp + log.gdpc
                   + (1 + brwmny | cntry.yr),
                   weights = dweight, data=tmp)
summary(fit.income)
# interplot (fit.income, "brwmny", "log.income")

#### Graph fit.income ####
# Find the omitted data
omitted.obs <- as.numeric (attr (fit.income@frame, "na.action"))
# Every survey has at least one omitted variable. Nine surveys are
# completely omitted. These are:
noData <- unique(tmp$cntry.yr)[!is.element(unique(tmp$cntry.yr), unique (tmp[-omitted.obs,]$cntry.yr))]

incomeSims <- seq(-3,4,length=50)
interceptPred <- ranef (fit.income)$cntry.yr[,grep("^brwmny.std$", names (fixef (fit.income)))] +
   fixef(fit.income)[grep("^brwmny.std$", names (fixef (fit.income)))]
slopePred <- ranef (fit.income)$cntry.yr[,grep("^brwmny.std:log.income.std$", names (ranef (fit.income)$cntry.yr))] + fixef(fit.income)[grep("^brwmny.std:log.income.std$", names (fixef (fit.income)))]
which.brwmny <- grep("^brwmny.std$", names (fixef(fit.income)))
which.interac <- grep("^brwmny.std:log.income.std$", names (fixef(fit.income)))
se.pointPred <- se.ranef (fit.income)$cntry.yr[,2]
margEffect <- fixef (fit.income)[which.brwmny] + fixef (fit.income)[which.interac]*incomeSims
se.margEffect <- sqrt(vcov(fit.income)[which.brwmny,which.brwmny] +
                            vcov(fit.income)[which.interac,which.interac]*(incomeSims^2) +
                            vcov(fit.income)[which.brwmny,which.interac]*2*incomeSims)

calculateImpact <- fixef (fit.income)[which.interac]*quantile(tmp$log.income.std, prob=0.75, na.rm=T) -
                     fixef (fit.income)[which.interac]*quantile(tmp$log.income.std, prob=0.25, na.rm=T)
calculateImpact/0.044 # denominator is standard deviation of cross-survey random effects for credit

# pdf (paste0(graphicsPath, "creditEffectVIncome.pdf"), h=7, w=10)
par (mar=c(3,4,0.5,0.5), las=0)
plot (margEffect~incomeSims, bty="n", type="n"
      , xlab=""
      , ylab="", pch=19, axes=F
      , ylim=c(-0.3,0.1), xlim=c(-2,3), cex.axis=0.8)
for (i in 1:length(interceptPred)){
   abline (a=interceptPred[i], b=slopePred[i], col="lightgray")
}
polygon( y=c(margEffect+1.96*se.margEffect, rev(margEffect-1.96*se.margEffect))
         , x=c(incomeSims, rev(incomeSims))
         , col="darkgray"
         , border=NA)
points (xy.coords(incomeSims, margEffect), type="l")
axis (2)
axis (1)
mtext (side=1, line=2
       , text="Income (log)")
mtext (side=2, line=3
       , text="Marginal pooled effect of access to credit")
mtext (side=2, line=2
       , text="on redistributive preference")
abline (h=0, lty=2)
# dev.off ()



# # Interaction with Thewissen-Rueda tertiles
# fit.incomeTR <- lmer(gincdif2 ~ brwmny + incomeTER.TR
#                      + brwmny:incomeTER.TR
#                      + male + agea + unemplindiv 
#                      + eduyrs2 + mbtru2 + rlgdgr 
#                      + socgdp + log.gdpc
#                      + (1 + brwmny  | cntry.yr),
#                      weights = dweight, data=tmp)
# summary(fit.incomeTR)


# Interaction with  quintiles takes forever with actual random coefficients
# fit.incomeQNT <- lmer(gincdif2 ~ brwmny + incomeQNT
#                       + brwmny:incomeQNT
#                       + male + agea + unemplindiv 
#                       + eduyrs2 + mbtru2 + rlgdgr 
#                       + socgdp + log.gdpc
#                       + (1 + brwmny | cntry.yr),
#                       weights = dweight, data=tmp)
# summary(fit.incomeQNT)

# tmp$lo.quintile <- car::recode (as.character(tmp$incomeQNT), "'1'=1; NA=NA; else=0") #& !is.na(tmp$incomeQNT), 1, ifelse (is.na(tmp$incomeQNT), NA, 0))
# tmp$hi.quintile <- car::recode (as.character(tmp$incomeQNT), "'5'=1; NA=NA; else=0") #& !is.na(tmp$incomeQNT), 1, ifelse (is.na(tmp$incomeQNT), NA, 0))
# tmp$lo.quintile <- (tmp$lo.quintile-mean(tmp$lo.quintile, na.rm=T))/sd(tmp$lo.quintile, na.rm=T)
# tmp$hi.quintile <- (tmp$hi.quintile-mean(tmp$hi.quintile, na.rm=T))/sd(tmp$hi.quintile, na.rm=T)
# tmp$interac.lo <- (tmp$brwmny.std*tmp$lo.quintile-mean(tmp$brwmny.std*tmp$lo.quintile,na.rm=T))/sd(tmp$brwmny.std*tmp$lo.quintile,na.rm=T)
# tmp$interac.hi <- (tmp$brwmny.std*tmp$hi.quintile-mean(tmp$brwmny.std*tmp$hi.quintile,na.rm=T))/sd(tmp$brwmny.std*tmp$hi.quintile,na.rm=T)
# 
# fit.incomeQNT.alt <- lmer(gincdif2 ~ brwmny.std + lo.quintile
#                       + interac.lo + hi.quintile
#                       + interac.hi
#                       + male + agea + unemplindiv 
#                       + eduyrs2 + mbtru2 + rlgdgr 
#                       + socgdp + log.gdpc
#                       + (1 + brwmny + lo.quintile
#                          + interac.lo + hi.quintile
#                          + interac.hi | cntry.yr),
#                       weights = dweight, data=tmp)
# summary(fit.incomeQNT.alt)

#### Graph fit.incomeQNT.alt ####
# Find the omitted data
# omitted.obs <- as.numeric (attr (fit.incomeQNT.alt@frame, "na.action"))
# # Every survey has at least one omitted variable. Nine surveys are
# # completely omitted. These are:
# noData <- unique(tmp$cntry.yr)[!is.element(unique(tmp$cntry.yr), unique (tmp[-omitted.obs,]$cntry.yr))]
# 
# pointPred.lo <- ranef (fit.incomeQNT.alt)$cntry.yr[,grep("^brwmny.std$", names (fixef (fit.incomeQNT.alt)))] +
#             fixef(fit.incomeQNT.alt)[grep("^brwmny.std$", names (fixef (fit.incomeQNT.alt)))] +
#             max(unique(tmp$lo.quintile), na.rm=T)* (ranef (fit.incomeQNT.alt)$cntry.yr[,grep("^interac.lo$", names (fixef (fit.incomeQNT.alt)))] + fixef(fit.incomeQNT.alt)[grep("^interac.lo$", names (fixef (fit.incomeQNT.alt)))])
# which.brwmny <- grep("^brwmny.std$", names (fixef(fit.incomeQNT.alt)))
# which.interac <- grep("interac.lo", names (fixef(fit.incomeQNT.alt)))
# se.pointPred.lo <- sqrt(vcov(fit.incomeQNT.alt)[which.brwmny,which.brwmny] +
#                         vcov(fit.incomeQNT.alt)[which.interac,which.interac]*(max(unique(tmp$lo.quintile), na.rm=T)^2) +
#                         vcov(fit.incomeQNT.alt)[which.brwmny,which.interac]*2*max(unique(tmp$lo.quintile), na.rm=T))
# margEffect.lo <- fixef (fit.incomeQNT.alt)[which.brwmny] + fixef (fit.incomeQNT.alt)[which.interac]*max(unique(tmp$lo.quintile), na.rm=T)
# se.margEffect.lo <- sd (pointPred.lo)
# 
# correctOrder <- order (pointPred.lo)
# 
# pointPred.hi <- ranef (fit.incomeQNT.alt)$cntry.yr[,grep("^brwmny.std$", names (fixef (fit.incomeQNT.alt)))] +
#             fixef(fit.incomeQNT.alt)[grep("^brwmny.std$", names (fixef (fit.incomeQNT.alt)))] +
#             (ranef (fit.incomeQNT.alt)$cntry.yr[,grep("^interac.hi$", names (fixef (fit.incomeQNT.alt)))] + fixef(fit.incomeQNT.alt)[grep("^interac.hi$", names (fixef (fit.incomeQNT.alt)))])*max(unique(tmp$hi.quintile), na.rm=T)
# which.interac <- grep("interac.hi", names (fixef(fit.incomeQNT.alt)))
# se.pointPred.hi <- sqrt(vcov(fit.incomeQNT.alt)[which.brwmny,which.brwmny] +
#                         vcov(fit.incomeQNT.alt)[which.interac,which.interac]*(max(unique(tmp$hi.quintile), na.rm=T)^2) +
#                         vcov(fit.incomeQNT.alt)[which.brwmny,which.interac]*2*max(unique(tmp$hi.quintile), na.rm=T))
# margEffect.hi <- fixef (fit.incomeQNT.alt)[which.brwmny] + fixef (fit.incomeQNT.alt)[which.interac]*max(unique(tmp$hi.quintile), na.rm=T)
# se.margEffect.hi <- sd (pointPred.hi)
# 
# 
# 
# # pdf (paste0(graphicsPath, "incomeQuintile.pdf"), h=7, w=10)
# par (mar=c(3,4,0.5,0.5), las=0)
# plot (pointPred.lo[correctOrder], bty="n"
#       , xlab=""
#       , ylab="", pch=19, axes=F
#       , ylim=c(-0.3,0.1), cex.axis=0.8)
# points (pointPred.hi[correctOrder] 
#         , pch=19, col="grey")
# axis (2)
# axis (1, at=c(1,length(pointPred)), labels=NA)
# mtext (side=1, line=1
#        , text="Country-years")
# mtext (side=2, line=3
#        , text="Marginal effect of access to credit")
# mtext (side=2, line=2
#        , text="on redistributive preference")
# abline (h=margEffect.lo)
# abline (h=margEffect.hi, col="gray")
# # polygon( y=c(margEffect.lo+1.96*se.margEffect.lo
# #              , margEffect.lo+1.96*se.margEffect.lo
# #              , margEffect.lo-1.96*se.margEffect.lo
# #              , margEffect.lo-1.96*se.margEffect.lo)
# #          , x=c(1, length(pointPred.lo), length(pointPred.lo), 1)
# #          , col="gray"
# #          , border=NA)
# # polygon( y=c(margEffect.hi+1.96*se.margEffect.hi
# #              , margEffect.hi+1.96*se.margEffect.hi
# #              , margEffect.hi-1.96*se.margEffect.hi
# #              , margEffect.hi-1.96*se.margEffect.hi)
# #          , x=c(1, length(pointPred.hi), length(pointPred.hi), 1)
# #          , col="gray"
# #          , border=NA)
# # points (xy.coords(giniSim, margEffect), type="l")
# # points (xy.coords(giniPoints, pointPred), pch=19)
# segments (x0=1:length(pointPred.lo), x1=1:length(pointPred.lo)
#           , y0=pointPred.lo[correctOrder] + 1.96*se.pointPred.lo
#           , y1=pointPred.lo[correctOrder] - 1.96*se.pointPred.lo)
# segments (x0=1:length(pointPred.hi), x1=1:length(pointPred.hi)
#           , y0=pointPred.hi[correctOrder] + 1.96*se.pointPred.hi
#           , y1=pointPred.hi[correctOrder] - 1.96*se.pointPred.hi, col="grey")
# abline (h=0, lty=2)
# legend ("bottomright", bty="n", legend=c("Poorest","Richest"), pch=19, col=c("black","grey"))
# # dev.off ()






#### Insurance-as-risk fork ####
fit.risk.macro <- lmer(gincdif2 ~ brwmny*epl
                 + log.income + male + agea + unemplindiv 
                 + eduyrs2 + mbtru2 + rlgdgr 
                 + socgdp + log.gdpc
                 + (1 + brwmny | cntry.yr),
                 weights = dweight, na.action="na.omit", data=tmp)
summary(fit.risk.macro)
interplot(fit.risk.macro, "brwmny", "epl")

#### Graph fit.risk.macro ####
# Find the omitted data
omitted.obs <- as.numeric (attr (fit.risk.macro@frame, "na.action"))
# Every survey has at least one omitted variable. Nine surveys are
# completely omitted. These are:
noData <- unique(tmp$cntry.yr)[!is.element(unique(tmp$cntry.yr), unique (tmp[-omitted.obs,]$cntry.yr))]

pointPred <- ranef (fit.risk.macro)$cntry.yr[,2] + fixef(fit.risk.macro)[grep("^brwmny$", names (fixef (fit.risk.macro)))]
se.pointPred <- se.ranef (fit.risk.macro)$cntry.yr[,2]
eplPoints <- as.numeric (unlist (by (tmp$epl, tmp$cntry.yr, unique)))
# Eliminate those countries that are not in the estimation
eplPoints <- eplPoints[is.element(unique(tmp$cntry.yr), unique (tmp[-omitted.obs,]$cntry.yr))]
correctOrder <- order (pointPred)

# marginal effect of access-to-credit conditional on gini: dY/dbrwmny|gininet
eplSim <- seq (0,5, length=50)
which.brwmny <- grep("^brwmny$", names (fixef(fit.risk.macro)))
which.interac <- grep("brwmny:epl", names (fixef(fit.risk.macro)))
margEffect <- fixef (fit.risk.macro)[which.brwmny] + fixef (fit.risk.macro)[which.interac]*eplSim
se.margEffect <- sqrt(vcov(fit.risk.macro)[which.brwmny,which.brwmny] 
                      + vcov(fit.risk.macro)[which.interac,which.interac]*(eplSim^2)
                      + vcov(fit.risk.macro)[which.brwmny,which.interac]*2*eplSim)

# pdf (paste0(graphicsPath, "eplEffect.pdf"), h=7, w=10)
par (mar=c(3,4,0.5,0.5), las=0)
plot (pointPred~eplPoints, type="n", bty="n"
      , xlab=""
      , ylab=""
      , ylim=c(-0.30,0.1), cex.axis=0.8)
mtext (side=1, line=2
       , text="Employment protection laws")
mtext (side=2, line=3
       , text="Marginal effect of access to credit")
mtext (side=2, line=2
       , text="on redistributive preference")
polygon( y=c(margEffect+1.96*se.margEffect, rev(margEffect-1.96*se.margEffect))
         , x=c(eplSim, rev(eplSim))
         , col="gray"
         , border=NA)
points (xy.coords(eplSim, margEffect), type="l")
points (xy.coords(eplPoints, pointPred), pch=19)
segments (x0=eplPoints, x1=eplPoints
          , y0=pointPred + 1.96*se.pointPred
          , y1=pointPred - 1.96*se.pointPred)
abline (h=0, lty=2)
# dev.off ()




# Change optimizer: doesn't work
# lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

# Rescale variables
tmp$risk.std <- (tmp$risk-mean(tmp$risk, na.rm=T))/sd(tmp$risk, na.rm=T)

fit.risk <- lmer(gincdif2 ~ brwmny.std + risk.std + brwmny.std:risk.std
                     + male + agea + log.income.std
                     + unemplindiv
                     + eduyrs2 + mbtru2+ rlgdgr 
  + socgdp + log.gdpc
  + (1 + brwmny | cntry.yr),
  weights = dweight, na.action="na.omit", data=tmp)

summary(fit.risk)

# interplot (fit.risk, "brwmny", "risk")

#### Graph fit.risk.macro ####
# Find the omitted data
omitted.obs <- as.numeric (attr (fit.risk@frame, "na.action"))
# Every survey has at least one omitted variable. Nine surveys are
# completely omitted. These are:
noData <- unique(tmp$cntry.yr)[!is.element(unique(tmp$cntry.yr), unique (tmp[-omitted.obs,]$cntry.yr))]

#### Graph fit.risk ####
riskSims <- seq (-3, 3, length=50)
interceptPred <- ranef (fit.risk)$cntry.yr[,grep("^brwmny.std$", names (fixef (fit.risk)))] +
   fixef(fit.risk)[grep("^brwmny.std$", names (fixef (fit.risk)))]
slopePred <- ranef (fit.risk)$cntry.yr[,grep("^brwmny.std:risk.std$", names (ranef (fit.risk)$cntry.yr))] + fixef(fit.risk)[grep("^brwmny.std:risk.std$", names (fixef (fit.risk)))]
which.brwmny <- grep("^brwmny.std$", names (fixef(fit.risk)))
which.interac <- grep("^brwmny.std:risk.std$", names (fixef(fit.risk)))
margEffect <- fixef (fit.risk)[which.brwmny] + fixef (fit.risk)[which.interac]*riskSims
se.margEffect <- sqrt(vcov(fit.risk)[which.brwmny,which.brwmny] +
                         vcov(fit.risk)[which.interac,which.interac]*(riskSims^2) +
                         vcov(fit.risk)[which.brwmny,which.interac]*2*riskSims)

calculateImpact <- fixef (fit.risk)[which.brwmny] + fixef (fit.risk)[which.interac]*quantile(tmp$risk.std, prob=0.75) -
   fixef (fit.risk)[which.brwmny] + fixef (fit.risk)[which.interac]*quantile(tmp$risk.std, prob=0.25)
calculateImpact/0.044 # denominator is standard deviation of cross-survey random effects for credit

# pdf (paste0(graphicsPath, "creditEffectVrisk.pdf"), h=7, w=10)
par (mar=c(3,4,0.5,0.5), las=0)
plot (margEffect~riskSims, bty="n", type="n"
      , xlab=""
      , ylab="", pch=19, axes=F
      , ylim=c(-0.3,0.1), xlim=c(-2,2), cex.axis=0.8)
for (i in 1:length(interceptPred)){
   abline (a=interceptPred[i], b=slopePred[i], col="lightgray")
}
polygon( y=c(margEffect+1.96*se.margEffect, rev(margEffect-1.96*se.margEffect))
         , x=c(riskSims, rev(riskSims))
         , col="darkgray"
         , border=NA)
points (xy.coords(riskSims, margEffect), type="l")
axis (2)
axis (1)
mtext (side=1, line=2
       , text="Risk of job elimination index")
mtext (side=2, line=3
       , text="Marginal pooled effect of access to credit")
mtext (side=2, line=2
       , text="on redistributive preference")
abline (h=0, lty=2)
# dev.off ()


# Final interaction 
fit.risk.rich <- lmer(gincdif2 ~ brwmny + risk + brwmny:risk
                 + male + agea + unemplindiv 
                 + eduyrs2 + mbtru2 + rlgdgr 
                 + socgdp + log.gdpc
                 + (1 + brwmny | cntry.yr),
                 weights = dweight, data=tmp, subset=incomeQNT==4 | incomeQNT==5)
summary(fit.risk.rich)


fit.risk.poor <- lmer(gincdif2 ~ brwmny + risk + brwmny:risk
                      + male + agea + unemplindiv 
                      + eduyrs2 + mbtru2 + rlgdgr 
                      + socgdp + log(gdpc)
                      + (1 + brwmny  | cntry.yr),
                      weights = dweight, data=tmp, subset=incomeQNT==1 | incomeQNT==2)
summary(fit.risk.poor)


#### Gather models in a small numer of tables ####


stargazer (fit.baseline, fit.gini, fit.risk)
stargazer (fit.income, fit.incomeTR, fit.incomeQNT, fit.risk.rich, fit.risk.poor)

# Dichotomous risk variable
tmp$dich.risk <- as.factor (ifelse (tmp$risk<0, 0, 1))
with (tmp, table(brwmny, incomeQNT, dich.risk))

fit.risk.full <- lmer(gincdif2 ~ brwmny*dich.risk*incomeQNT
                      + male + agea + unemplindiv 
                      + eduyrs2 + mbtru2 + rlgdgr 
                      + socgdp + log.gdpc
                      + (1 | cntry.yr),
                 weights = dweight, data=tmp)
summary(fit.risk.full)

sims <- rmvnorm (100, mean=fixef(fit.risk.full)
                 , sigma=as.matrix(vcov(fit.risk.full)))

x.male <- max (tmp$male, na.rm=T)
x.agea <- median (rescale.var(tmp$agea), na.rm=T)
x.unemplindiv <- max (tmp$unemplindiv, na.rm=T)
x.eduyrs2 <- median (rescale.var(tmp$eduyrs2))
x.mbtru2  <- max (tmp$mbtru2, na.rm=T)
x.rlgdgr  <- median (rescale.var(tmp$rlgdgr), na.rm=T)
x.socgdp  <- mean (rescale.var(tmp$socgdp), na.rm=T)
x.loggdp  <- mean (log(tmp$gdpc))
constant <- c(x.male, x.agea, x.unemplindiv, x.eduyrs2, x.mbtru2
              , x.rlgdgr, x.socgdp, x.loggdp)

brwmny <- c(1,5)
risk   <- c(0,1)
income.1 <- c(1,0,0,0,0)
income.2 <- c(0,1,0,0,0)
income.4 <- c(0,0,0,1,0)
income.5 <- c(0,0,0,0,1)

hat.y <- array (NA, dim=c(2,2,5))
for (i in 1:2){
   for (j in 1:2){
      for (k in 1:5){
         obs.values <- c(1, brwmny[i], risk[j], income.1[k], income.2[k]
                         , income.4[k], income.5[k]
                         , constant
                         , brwmny[i]*risk[j], brwmny[i]*income.1[k]
                         , brwmny[i]*income.2[k], brwmny[i]*income.4[k]
                         , brwmny[i]*income.5[k], risk[j]*income.1[k]
                         , risk[j]*income.2[k], risk[j]*income.4[k]
                         , risk[j]*income.5[k], brwmny[i]*risk[j]*income.1[k]
                         , brwmny[i]*risk[j]*income.2[k]
                         , brwmny[i]*risk[j]*income.4[k]
                         , brwmny[i]*risk[j]*income.5[k])
         hat.y[i,j,k] <-  mean (sims %*% obs.values)
      }
   }
}
round (hat.y,2)

###########################################################

# ROBUSTNESS

# 1) social network effects
fit.socialorigin <- lmer(gincdif2 ~ brwmny.std + socialorigin
  + log.income.std + male + agea + unemplindiv 
  + eduyrs2 + mbtru2 + rlgdgr 
  + socgdp + log.gdpc
  + (1 + brwmny | cntry.yr),
  weights = dweight, data=tmp)
summary(fit.socialorigin)

# 1) wealth effects
fit.wealth <- lmer(gincdif2 ~ brwmny + h.owner
  + log.income + male + agea + unemplindiv 
  + eduyrs2 + mbtru2 + rlgdgr 
  + socgdp + log.gdpc
  + (1 + brwmny | cntry.yr),
  weights = dweight, data=tmp, subset = essround == 2)
summary(fit.wealth)

fit.robustFull <- lmer(gincdif2 ~ brwmny + h.owner + socialorigin
  + log.income + male + agea + unemplindiv 
  + eduyrs2 + mbtru2 + rlgdgr 
  + socgdp + log.gdpc
  + (1 + brwmny + h.owner + socialorigin | cntry.yr),
  weights = dweight, data=tmp, subset = essround == 2)
summary(fit.robustFull)

# export findings
stargazer(fit.socialorigin, fit.wealth, fit.robustFull)

###########################################################
#### RUN REGRESSIONS BY "CNTRY.YR" and STORE ESTIMATES ####
###########################################################

# generate variable indicating missing income information
tmp <- tmp %>%
  group_by(cntry.yr) %>%
  mutate(n_unique = n_distinct(incomeQNT))
tmp$incomeNA <- ifelse(tmp$n_unique > 2, 0,1)

# with "incomeQNT"
incomeQNT_results <- dlply(tmp[tmp$incomeNA == 0,], "cntry.yr", function(df)
  lm(gincdif2 ~ brwmny + incomeQNT + brwmny : incomeQNT
    + male + agea + unemplindiv 
    + eduyrs2 + mbtru2 + rlgdgr 
    + socgdp + log.gdpc, data = df) 
)

# with "incomeTER"
incomeTER_results <- dlply(tmp[tmp$incomeNA == 0,], "cntry.yr", function(df)
  lm(gincdif2 ~ brwmny + incomeTER + brwmny : incomeTER
    + male + agea + unemplindiv 
    + eduyrs2 + mbtru2 + rlgdgr 
    + socgdp + log.gdpc, data = df) 
)

# with "log.income"
incomeLog_results <- dlply(tmp[tmp$incomeNA == 0,], "cntry.yr", function(df)
  lm(gincdif2 ~ brwmny + log.income
  + brwmny:log.income
  + male + agea + unemplindiv 
  + eduyrs2 + mbtru2 + rlgdgr 
  + socgdp + log.gdpc, data = df)
  )

# collect all results in single dataframe

incomeQNT <- lapply(incomeQNT_results, function(x) coef(x)[grep ("brwmny", names(coef(x)))])
incomeLog <- lapply(incomeLog_results, function(x) coef(x)[grep ("brwmny", names(coef(x)))])
incomeTER <- lapply(incomeTER_results, function(x) coef(x)[grep ("brwmny", names(coef(x)))])

incomeQNT.vcov <- lapply(incomeQNT_results, function(x) vcov(x)[grep ("brwmny", colnames(vcov(x))), grep ("brwmny", colnames(vcov(x)))])
incomeLog.vcov <- lapply(incomeLOG_results, function(x) vcov(x)[grep ("brwmny", colnames(vcov(x))), grep ("brwmny", colnames(vcov(x)))])
incomeTER.vcov <- lapply(incomeTER_results, function(x) vcov(x)[grep ("brwmny", colnames(vcov(x))), grep ("brwmny", colnames(vcov(x)))])



















# # With dichotomous "support for welfare state" outcome variable
# fit.risk.dich <- glmer(gincdif2bin ~ brwmny*dich.risk*incomeQNT.TR
#                       + male + agea + unemplindiv 
#                       + eduyrs2 + mbtru2 + rlgdgr 
#                       + socgdp + log(gdpc)
#                       + (1 | cntry.yr)
#                       , family=binomial
#                       , weights = dweight, data=tmp)
# summary(fit.risk.dich)
# # Model does not converge after a long while

###########################################################################
###########################################################################

# INCOME ANALYSIS: redistribution story, micro-level step





























###########################################################################
###########################################################################
# Comment these lines out if need to run model on full sample
# Otherwise, the following lines keep a smaller sample
set.seed (2018)
fullEss <- ess
smallEss <- ess %>%
  group_by(cntry.yr) %>%
  sample_n(50)
# ess <- smallEss  # Uncomment for Bayes
###########################################################################
###########################################################################

load ("~/Dropbox/CreditPreferences/Data/InequalityData/fullInequalityData.RData")
datInequal$cntry.yr <- paste0(datInequal$ISO3, datInequal$year)
datInequal <- datInequal[order(datInequal$country, datInequal$year),]

# Countries with inequality data but no survey data
extraCountries <- unique(datInequal$ISO3[!is.element(datInequal$ISO3, ess$cou)])
# Countries with survey data but no inequality data
unique(ess$cou[!is.element(ess$cou, datInequal$ISO3)]) # No inequality data for SVK and for LUX

datInequal <- datInequal[!is.element(datInequal$ISO3,extraCountries),]

# # These observations have duplicates # NO LONGER ANY DUPLICATES
# # DEU2004 FIN2002 FIN2004 GBR2002 GBR2004 IRL2004 SWE2004
# datInequal$cntry.yr[duplicated (datInequal[,c(1,2)])]
# # We need to check this; the problem is with the oecd data,
# # which reports two different sets of percentile ratios for
# # these observations. Which one is correct? Who knows. Here,
# # I only use the first observation

# datInequal <- datInequal[!duplicated (datInequal[,c(1,2)]) ,]



# Build a function to extend past values forward
extendFwd <- function (x) {
  length.x <- length (x)
  new.x <- c()
  new.x[1] <- x[1]
  for (i in 2:length(x)){
    new.x[i] <- ifelse (!is.na(x[i]), x[i], new.x[i-1])
  }
  return (new.x)
}

countries <- unique (datInequal$ISO3)
lis.90.10 <- lis.90.50 <- wiid.90.50 <- wiid.50.10 <- c()
for (i in 1:length(countries)){
  ctr <- countries[i]
  tmp1 <- extendFwd(datInequal$p90p10_lis[datInequal$ISO3==ctr])
  tmp2 <- extendFwd(datInequal$p90p50_lis[datInequal$ISO3==ctr])
  tmp3 <- extendFwd(datInequal$p90p50_wiid[datInequal$ISO3==ctr])
  tmp4 <- extendFwd(datInequal$p50p10_wiid[datInequal$ISO3==ctr])
  lis.90.10 <- c(lis.90.10, tmp1)
  lis.90.50 <- c(lis.90.50, tmp2)
  wiid.90.50 <- c(wiid.90.50, tmp3)
  wiid.50.10 <- c(wiid.50.10, tmp4)
}
skew.90.50.10 <- wiid.90.50/wiid.50.10

attach (datInequal)
Inequality <- data.frame (cntry.yr=cntry.yr, year=year, gini.net.solt=gini_net_solt
                          , gini.mkt.solt=gini_market_solt
                          , lis.90.10=lis.90.10, lis.90.50=lis.90.50
                          , wiid.50.10=wiid.50.10, wiid.90.50=wiid.90.50
                          , skew.90.50.10=skew.90.50.10)
detach (datInequal)

Inequality <- Inequality[is.element(Inequality$year, unique (ess$year)),]

# Merge survey-level data into individual-level data
ess <- merge (ess, Inequality, by="cntry.yr", all.x=T)

# Check Solt measures (we have two, one from Thewissen-Rueda, one from Jonas)
with (ess,
      plot (gininet~gini.net.solt, pch=19))
with (ess,
      plot (ginimarket~gini.mkt.solt, pch=19))
# don't quite line up, but there are no major deviations either




# Build group counter (grouping variable is municipality)
group.counter.individual <- as.numeric (as.factor (ess$cntry.yr))
by (group.counter.individual, ess$cntry.yr, unique) # check that string order is preserved
ess$group.counter.individual <- group.counter.individual

# Some questions:
# what is mbtru2?
# what is rlgdgr, and why does it have values 77,88,99?
covariates <- c("gincdif2","eduyrs2","brwmny","perceqnatincdollar","mbtru2"
                ,"rlgdgr","unemplindiv","male"
                ,"agea","country","year.x")
smallData <- ess[,is.element(colnames(ess), covariates)]

# country.covariates (keep separate)
country.covariates <- c("gininet","ginimarket","gdpc"
                        ,"unempl","wiid.50.10","wiid.90.50")
countryData <- ess[,is.element(colnames(ess)
                                     , country.covariates)]

a.out <- amelia(smallData
                , m=5
                , id=c("country","year.x")
                , noms=c("male","mbtru2","unemplindiv")
                , ords=c("gincdif2","brwmny","rlgdgr")
                , bounds=matrix(data=c(c(3,1,11) # agea
                                       , c(8,0,25) # eduyrs: capped at 25
                                       , c(11,0,20)) # income (may need to be done categorically)
                                , byrow=T, nrow=3)
                , p2s=1
                , parallel="multicore")

summary (a.out$imputations[[3]])

a.out.list <- list(na.omit(a.out$imputations[[1]]),
                   na.omit(a.out$imputations[[2]]),
                   na.omit(a.out$imputations[[3]]),
                   na.omit(a.out$imputations[[4]]),
                   na.omit(a.out$imputations[[5]]))

short.ess <- a.out.list[[3]]
# JM: don't we need to use all 5 imputations for all models following this?

# Outcome variable
# recode outcome, since there are very few 1s
short.ess$gincdif2 <- car::recode (short.ess$gincdif2, "2=1; 3=2; 4=3; 5=4")  # this is against the literature standard that always uses 5-point scale: what's the advantage of combining 2 and 1?
y <- as.numeric (short.ess$gincdif2)
hist (y, breaks=c(0.5,1.5,2.5,3.5,4.5))

# Control variables: Predictors with completely-pooled effects
pooled.controls <- c("perceqnatincdollar", "male", "agea", "unemplindiv")
X <- apply (short.ess[,is.element(colnames(short.ess), pooled.controls)], 2, function(x) (x-mean(x))/sd(x)  )
X.names <- c("income","male","age","unemployed")


# Control variables: Predictors with unmodeled random coefficients across country years
controls <- c("eduyrs2", "mbtru2", "rlgdgr")
W <- apply (short.ess[,is.element(colnames(short.ess), controls)], 2, function(x) (x-mean(x, na.rm=T))/sd(x, na.rm=T)) # standardize all predictors
W.names <- c("education","union","religiosity") # keep names, to plot later on

# Main predictor, with modeled random coefficient
V <- (short.ess$brwmny-mean(short.ess$brwmny, na.rm=T)) / sd (short.ess$brwmny, na.rm=T)
V.name <- "borrow money"

# To revise later: I submit that we don't need survey-level averages of individual-level
# variables for which we estimate a single completely-pooled coefficient
# Add country-year-level means for all individual-level predictors (these are controls)
m.eduyrs <- as.numeric ( by (short.ess$eduyrs2
                             , group.counter.individual, mean, na.rm=T))
# m.income <- as.numeric ( by (short.ess$perceqnatincdollar
#                              , group.counter.individual, mean, na.rm=T))
m.church <- as.numeric ( by (short.ess$rlgdgr
                             , group.counter.individual, mean, na.rm=T))
# m.male   <- as.numeric ( by (short.ess$male
#                              , group.counter.individual, mean, na.rm=T))
# m.age    <- as.numeric ( by (short.ess$agea
#                              , group.counter.individual, mean, na.rm=T))
# m.unemp  <- as.numeric ( by (short.ess$unemplindiv
#                              , group.counter.individual, mean, na.rm=T))
m.mbtru  <- as.numeric ( by (short.ess$mbtru2
                             , group.counter.individual, mean, na.rm=T))


#########################################################################
#### LMER RUNS ####
#########################################################################
# MODEL ESTIMATION #############

# baseline model
# completely-pooled effects (W): "perceqnatincdollar", "male", "agea", "unemplindiv"
# unmodeled random effects (X): "eduyrs2", "mbtru2", "rlgdgr"
# level-2 predictors (ess): ginimarket (gini.net.solt, wiid.50.10, wiid.90.50)
#                     unempl, gdpc, socgdp
# Outcome: y
# Brwmny:  V

standardize <- function (x) {(x -mean(x))/sd(x)}

modelData <- data.frame (y=y, brwmny=short.ess$brwmny
                         , ineqmkt=I(ess$ginimarket/100)
                         , ineqnet=I(ess$gininet/100)
                         , hi.end=ess$wiid.90.50
                         , lo.end=ess$wiid.50.10
                         , skew=ess$skew.90.50.10
                         , welfare=standardize(ess$socgdp)
                         , develop=standardize(ess$gdpc)
            #              , ineqmkt=I(standardize(short.ess$ginimarket))
            # , ineqnet=I(standardize(short.ess$gininet))
            # , skew=I(standardize(short.ess$skew.90.50.10))
            # , welfare=I(standardize(short.ess$socgdp))
            # , develop=I(standardize(short.ess$gdpc))
            , income=X[,4]
            , male=X[,2], age=X[,1], unemploy=X[,3]
            , educ=W[,2], union=W[,3], religion=W[,1]
            , cntry.yr=ess$cntry.yr, dweight=ess$dweight)

ineqmkt <- as.numeric (by(modelData$ineqmkt, modelData$cntry.yr, unique))
ineqnet <- as.numeric (by(modelData$ineqnet, modelData$cntry.yr, unique))
skew    <- as.numeric (by(modelData$skew, modelData$cntry.yr, unique))
hi.end  <- as.numeric (by(modelData$hi.end, modelData$cntry.yr, unique))
lo.end  <- as.numeric (by(modelData$lo.end, modelData$cntry.yr, unique))

fit.baseline <- lmer(y ~ brwmny
                + income + male + age + unemploy 
                + educ + union + religion 
                + welfare + develop
                + (1 + brwmny | cntry.yr), 
                # + (1 + brwmny + educ + union + religion | cntry.yr), 
                weights = dweight, data=modelData)

fit.ineqmkt <- lmer(y ~ brwmny*ineqmkt
                 + brwmny*welfare
                 + brwmny*develop
                 + income + male + age + unemploy 
                 + educ + union + religion 
                 + welfare + develop
                 + (1 + brwmny | cntry.yr), 
                 # + (1 + brwmny + educ + union + religion | cntry.yr), 
                 weights = dweight, data=modelData)

fit.ineqnet <- lmer(y ~ brwmny*ineqnet
                    + brwmny*welfare
                    + brwmny*develop
                    + income + male + age + unemploy 
                    + educ + union + religion 
                    + welfare + develop
                    + (1 + brwmny | cntry.yr), 
                    # + (1 + brwmny + educ + union + religion | cntry.yr), 
                    weights = dweight, data=modelData)

fit.skew <- lmer(y ~ brwmny*skew
                    + brwmny*welfare
                    + brwmny*develop
                    + income + male + age + unemploy 
                    + educ + union + religion 
                 + welfare + develop
                 + (1 + brwmny | cntry.yr), 
                    # + (1 + brwmny + educ + union + religion | cntry.yr), 
                    weights = dweight, data=modelData)

fit.hiend <- lmer(y ~ brwmny*hi.end
                 + brwmny*welfare
                 + brwmny*develop
                 + income + male + age + unemploy 
                 + educ + union + religion 
                 + welfare + develop
                 + (1 + brwmny | cntry.yr), 
                 # + (1 + brwmny + educ + union + religion | cntry.yr), 
                 weights = dweight, data=modelData)

fit.loend <- lmer(y ~ brwmny*lo.end
                  + brwmny*welfare
                  + brwmny*develop
                  + income + male + age + unemploy 
                  + educ + union + religion 
                  + welfare + develop
                  + (1 + brwmny | cntry.yr), 
                  # + (1 + brwmny + educ + union + religion | cntry.yr), 
                  weights = dweight, data=modelData)

fit.bothend <- lmer(y ~ brwmny*lo.end
                    + brwmny*hi.end
                  + brwmny*welfare
                  + brwmny*develop
                  + income + male + age + unemploy 
                  + educ + union + religion 
                  + welfare + develop
                  + (1 + brwmny | cntry.yr), 
                  # + (1 + brwmny + educ + union + religion | cntry.yr), 
                  weights = dweight, data=modelData)

fit.fullTest <- lmer(y ~ brwmny*lo.end
  + brwmny*hi.end
  + brwmny*ineqnet
  + brwmny*welfare
  + brwmny*develop
  + income + male + age + unemploy 
  + educ + union + religion 
  + welfare + develop
  + (1 + brwmny | cntry.yr), 
  # + (1 + brwmny + educ + union + religion | cntry.yr), 
  weights = dweight, data=modelData)


summary(fit.baseline)
summary(fit.ineqmkt)
summary(fit.ineqnet)
summary(fit.skew)
summary(fit.hiend)
summary(fit.loend)
summary(fit.bothend)
summary(fit.fullTest)

covar.names <- c("Borrowing Ability", "Gini", "90-50 ratio", "50-10 ratio", "Income"
  , "Male", "Age", "Unemployed", "Years Schooling", "Union Membership"
  , "Religiousness", "Borrowing Ability * Gini", "Borrowing Ability * 90-50 ratio"
  , "Borrowing Ability * 50-10 ratio", "Borrowing Ability * Social Expenditure"
  , "Borrowing Ability * per capita GDP", "Social Expenditure", "per capita GDP")

stargazer (fit.baseline, fit.ineqnet, fit.hiend, fit.loend, fit.bothend
  , label = "T:resultsMain"
  , title = "Multilevel linear model of the effect of \emph{access to credit} on \emph{preference for redistribution}. Statistics are restricted maximum likelihood estimates. Individuals are nested with 83 level-2 groups, corresponding to different country-year pairs."
  , font.size = "tiny"
  , dep.var.caption = "Redistribution Demand"
  , covariate.labels = covar.names
  , align = T)



#### Graphs ####
graphPath <- paste0(getwd (),"/Draft/draftMPSA/")

# Plot country average support for redistribution
tmp$cou <- substr(tmp$cntry.yr, start=1,stop=3)
tmp <- merge(tmp, unique(ess[,c("cou", "coun")]), by = "cou")
redist.cntry <- prop.table(table(tmp$coun, tmp$gincdif2), margin = 1)
colnames(redist.cntry) <- c("strongly disagree", "disagree", "neither", "agree", "strongly agree")

redist.cntry <- redist.cntry[order(redist.cntry[,5], redist.cntry[,4]),]

pdf(paste0(graphPath, "redistByCntry.pdf"), h=5, w=7)
par(mar = c(4, 4, 0, 0))
barplot(as.matrix(t(redist.cntry)), horiz = T, cex.names = 0.5, axes = F, las = 2)
par(fig = c(0, 1, 0, 1), oma = c(0,0,6,0), mar = c(0, 0, 0, 0), new = TRUE)
legend("bottom", colnames(redist.cntry), xpd=T, horiz = F
  , fill = gray.colors(5)
  , ncol = 5, cex = 0.7, bty ="n")
dev.off()

# Credit Access-Fin Depth (Appendix Figure)
crdt2gdp <- read_excel("~/Dropbox/aftermathCrisis/Data/FinStructure_2016.xlsx", sheet = 3)[,c("cncode","year","pcrdbgdp")]
tmp <- ddply(tmp,.(cntry.yr),transform,m.brwmny = mean(brwmny, na.rm = T))
tmp$year <- ifelse(tmp$essround ==1, 2002
  , ifelse(tmp$essround== 2, 2004
    , ifelse(tmp$essround==3, 2006
      ,ifelse(tmp$essround==4,2008, 2010))))

tmp <- merge(tmp, crdt2gdp, by.x=c("cou","year"), by.y= c("cncode", "year"))

data <- unique(tmp[,c("pcrdbgdp", "m.brwmny", "cntry.yr")])
reg <- lm(m.brwmny ~ pcrdbgdp, data = data)
newx = seq(min(data$pcrdbgdp, na.rm= T),max(data$pcrdbgdp, na.rm= T),by = 0.05)
conf_interval <- predict(reg, newdata=data.frame(pcrdbgdp=newx), interval="confidence",
  level = 0.95)

pdf(paste0(graphPath, "CorrCreditFindepth.pdf"), h=5, w=7)
plot(data$pcrdbgdp, data$m.brwmny, type = "p", axes = F,
  ylab = "Mean Reported Access to Credit", xlab = "Bank Credit to GDP (%)", 
  ylim = c(0,5), xlim = c(0,220), pch="+")
axis(2, seq(0,5,0.5))
axis(1, seq(0,200, 50), pos =0)
abline(reg)
matlines(newx, conf_interval[,2:3], col = "blue", lty=2)
dev.off()

# Credit Access by Year (Appendix Table)
tmp$crdt <- ifelse(tmp$brwmny == 4 | tmp$brwmny == 5, 1
  , ifelse(tmp$brwmny >= 1 & tmp$brwmny <= 3, 0, NA))
prop.table(table(tmp$essround, tmp$crdt), margin = 1)

# Credit Access by Cnty-Yr (Appendix Table)
# average crdtAcs & Redistribution
crdt.CntryYr <- ddply(tmp, .(cntry.yr), summarize, m.crdt=mean(crdt, na.rm=T))
crdt.CntryYr$cntry <- substr(crdt.CntryYr$cntry.yr, start = 1, stop = 3)
crdt.CntryYr$year <- substr(crdt.CntryYr$cntry.yr, start = 4, stop = 7)
crdt.CntryYr$cntry.yr <- NULL

# long2wide
crdt.CntryYr <- reshape(crdt.CntryYr, idvar="cntry",timevar="year",direction="wide")

# plot crdtAcs by cntry.yr
pdf(paste0(graphPath, "crdtAcs_cntryYr.pdf"), h=8, w=14)
barplot(as.matrix(t(crdt.CntryYr[, c(2:ncol(crdt.CntryYr))])), beside = T, col=grey.colors(5),
  ylab="Ability to Borrow", names.arg = crdt.CntryYr$cntry, ylim=c(0,1))
box(bty="l")
legend("topleft", legend = c("2002", "2004", "2006", "2008", "2010"), fill = grey.colors(5)
  , cex = 1, bty = "n")
dev.off()



# Plot the random coefficients for brwmny, along with standard errors
# BAsed on model fit.baseline
ranef.brwmny <- coef(fit.baseline)$cntry.yr[,2]
se.ranef.brwmny <- se.ranef(fit.baseline)$cntry.yr[,2]

increasing.order <- order (ranef.brwmny)

pdf(paste0(graphPath, "randomCoefBaselineLMER.pdf"), h=7, w=11)
plot (c(1,length(ranef.brwmny))
      , c(min(ranef.brwmny-2*se.ranef.brwmny),max(ranef.brwmny+2*se.ranef.brwmny))
      , type="n", xlab="Survey (country-year)", ylab="Random coefficient of access to credit", axes=F)
segments (x0=1:length(ranef.brwmny), x1=1:length(ranef.brwmny)
          , y0=ranef.brwmny[increasing.order]-1.96*se.ranef.brwmny[increasing.order]
          , y1=ranef.brwmny[increasing.order]+1.96*se.ranef.brwmny[increasing.order])
points (xy.coords(1:length(ranef.brwmny), ranef.brwmny[increasing.order]), pch=19)
axis (2)
mtext (tolower(names (se.ranef.brwmny))[increasing.order], side=1, line=0, las=2, at=1:length(ranef.brwmny), cex=0.8)
abline (h=0, lty=3)
dev.off()

# Preliminary graphs
plot (ineqmkt, ranef.brwmny, pch=19); abline(lm(ranef.brwmny~ineqmkt), col="red")
plot (ineqnet, ranef.brwmny, pch=19); abline(lm(ranef.brwmny~ineqnet), col="red")
plot (skew, ranef.brwmny, pch=19); abline(lm(ranef.brwmny~skew), col="red")
plot (hi.end, ranef.brwmny, pch=19); abline(lm(ranef.brwmny~hi.end), col="red")
plot (lo.end, ranef.brwmny, pch=19); abline(lm(ranef.brwmny~lo.end), col="red")



# Plot the random coefficients for brwmny, along with standard errors
# based on model fit.ineqnet
M4 <- fit.ineqnet
b.hat.M4 <- fixef(M4)[2] + fixef(M4)[13]*ineqnet + ranef(M4)$cntry.yr[,2]
b.se.M4 <- se.ranef(M4)$cntry.yr[,2]

# plot on Figure 13.2(a)
# a.hat.M4 <- fixef(M4)[1] + fixef(M4)[3]*ineqnet + ranef(M4)$cntry.yr[,1]
# a.se.M4 <- se.ranef(M4)$cntry.yr[,1]
# lower <- a.hat.M4 - a.se.M4
# upper <- a.hat.M4 + a.se.M4
# par (mar=c(5,5,1,1)+.1)
# plot (ineqnet, a.hat.M4, cex.lab=1, cex.axis=1.1, ylim=range(lower,upper), 
#       xlab="Gini (net)", ylab="Random ''access to credit'' coefficients", 
#       pch=20)
# axis (2, c(0,1,1.5,2))
# curve (fixef(M4)[1] + fixef(M4)[3]*x, lwd=1, col="black", add=TRUE)
# segments (ineqnet, lower, ineqnet, upper, lwd=.5, col="gray10")
# mtext ("Intercepts", line=1)

# plot on Figure 13.2(b)
par (mar=c(5,5,1,1)+.1)
lower <- b.hat.M4 - b.se.M4
upper <- b.hat.M4 + b.se.M4
pdf(paste0(graphPath, "randomCoefIneqnetLMER.pdf"), h=7, w=11)
par (mar=c(5,5,4,2)+.1)
plot (ineqnet, b.hat.M4, cex.lab=1.2, cex.axis=1.1, ylim=range(lower,upper),
      xlab="Gini (net)", ylab="Random ''access to credit'' coefficients", pch=20, col="gray")
curve (fixef(M4)[2] + fixef(M4)[13]*x, lwd=3, col="black", add=TRUE)
segments (ineqnet, lower, ineqnet, upper, lwd=.5, col="gray")
abline (h=0, lty=3)
# mtext ("Slopes", line=1)
# add confidence interval for level-2 regression
# curve (fixef(M4)[2] + fixef(M4)[11]*x
#        + sqrt (vcov(M4)[2,2] + vcov(M4)[11,11]*x + 2*vcov(M4)[2,11]*x)
#        , lty=3, col="black", add=TRUE)
# curve (fixef(M4)[2] + fixef(M4)[11]*x
#        - sqrt (vcov(M4)[2,2] + vcov(M4)[11,11]*x + 2*vcov(M4)[2,11]*x)
#        , lty=3, col="black", add=TRUE)
dev.off()


# based on model fit.loend
par (mar=c(5,5,1,1)+.1)
M4 <- fit.loend
b.hat.M4 <- fixef(M4)[2] + fixef(M4)[13]*lo.end + ranef(M4)$cntry.yr[,2]
b.se.M4 <- se.ranef(M4)$cntry.yr[,2]
lower <- b.hat.M4 - b.se.M4
upper <- b.hat.M4 + b.se.M4
par (mar=c(5,5,4,2)+.1)
pdf(paste0(graphPath, "randomCoefLoendLMER.pdf"), h=7, w=11)
plot (lo.end, b.hat.M4, cex.lab=1.2, cex.axis=1.1, ylim=range(lower,upper),
      xlab="50/10 ratio", ylab="Random ''access to credit'' coefficients", pch=20, col="gray")
curve (fixef(M4)[2] + fixef(M4)[13]*x, lwd=3, col="black", add=TRUE)
segments (lo.end, lower, lo.end, upper, lwd=.5, col="gray")
abline (h=0, lty=3)
dev.off()

# based on model fit.hiend
par (mar=c(5,5,1,1)+.1)
M4 <- fit.hiend
b.hat.M4 <- fixef(M4)[2] + fixef(M4)[13]*hi.end + ranef(M4)$cntry.yr[,2]
b.se.M4 <- se.ranef(M4)$cntry.yr[,2]
lower <- b.hat.M4 - b.se.M4
upper <- b.hat.M4 + b.se.M4
par (mar=c(5,5,4,2)+.1)
pdf(paste0(graphPath, "randomCoefHiendLMER.pdf"), h=7, w=11)
plot (hi.end, b.hat.M4, cex.lab=1.2, cex.axis=1.1, ylim=range(lower,upper),
      xlab="90/50 ratio", ylab="Random ''access to credit'' coefficients", pch=20, col="gray")
curve (fixef(M4)[2] + fixef(M4)[13]*x, lwd=3, col="black", add=TRUE)
segments (hi.end, lower, hi.end, upper, lwd=.5, col="gray")
abline (h=0, lty=3)
dev.off()


# based on model fit.hiend
par (mar=c(5,5,1,1)+.1)
M4 <- fit.bothend
b.hat.M4 <- fixef(M4)[2] + fixef(M4)[15]*hi.end + ranef(M4)$cntry.yr[,2]
b.se.M4 <- se.ranef(M4)$cntry.yr[,2]
lower <- b.hat.M4 - b.se.M4
upper <- b.hat.M4 + b.se.M4
par (mar=c(5,5,4,2)+.1)
pdf(paste0(graphPath, "randomCoefBothendLMER.pdf"), h=7, w=11)
plot (hi.end, b.hat.M4, cex.lab=1.2, cex.axis=1.1, ylim=range(lower,upper),
      xlab="90/50 ratio", ylab="Random ''access to credit'' coefficients", pch=20, col="gray")
curve (fixef(M4)[2] + fixef(M4)[15]*x, lwd=3, col="black", add=TRUE)
segments (hi.end, lower, hi.end, upper, lwd=.5, col="gray")
abline (h=0, lty=3)
dev.off()

#########################################################################
#########################################################################
#########################################################################



#### Bayesian fit ####
# Variables at level-2


# Add other country-year variables
gini      <- as.numeric ( by (ess$ginimarket
                              , group.counter.individual, unique))
pcGDP     <- as.numeric ( by (ess$gdpc
                              , group.counter.individual, unique))

# country/year-level predictors, including all country/year-level means
D <- data.frame (m.eduyrs=m.eduyrs
                 # , m.income=m.income
                 , m.church=m.church
                 # , m.male=m.male
                 # , m.age=m.age
                 # , m.unemp=m.unemp
                 , m.mbtru=m.mbtru
                 , gini=gini
                 , pcGDP=log(pcGDP))  # Note: logged values
#, soc2gdp=soc2gdp
#, econ.conf=econ.conf  # lots of missing values
#, crdt2gdp=log(crdt2gdp)  # Note: logged values 
#, finInst=finInst)

D <- apply (D, 2, function(x) (x-mean(x))/sd(x)  )
D <- cbind (rep(1,nrow(D)), D)
D.names <- colnames (D); D.names[1] <- "Intercept"



# Gather other objects for jags run
n.cut <- max (y, na.rm=T) - 1
N <- length (y)
R <- max (group.counter.individual)
B <- ncol(X)
Tau.theta <- Tau.delta <- diag(ncol(D))*1
mu.theta  <- mu.delta  <- rep(0,ncol(D))

# Set parameter monitors
jags.parameters <- c("theta", "beta"
  , "delta"
  , "deviance"
  , "gamma"
  , "alpha"#
  , "kappa"
  , "tau","mu.gamma"
  , "sigmaAlpha", "sigmaKappa", "sigmaGamma"
  , "mu.gamma")   #, "Y.star")

# Collect data into dump format
jags.data <- dump.format(list(y=y, WCol=ncol(W), DCol=ncol(D)
                              , n.cut=n.cut, V=V, X=X
                              , W=W, N=N, R=R, D=D, B=B
                              , group=group.counter.individual
                              , mu.delta=mu.delta, Tau.delta=Tau.delta
                              , mu.theta=mu.theta, Tau.theta=Tau.theta)) 

# Function to create initial values
jags.inits <- function() {
  dump.format(
    list(
       #%#%#%tau.unsorted=sort(runif(n.cut,-2,2))
      #%#%#%, alpha=rnorm(1)
      tau.unsorted=cbind(runif(R,-2,2), runif(R,-2,2), runif(R,-2,2))
      , alpha=rnorm(R,0,1)
      , kappa=rnorm(R,0,1)
      , delta=rnorm (ncol(D),0,1)
      , theta=rnorm (ncol(D),0,1)
      , mu.gamma=rnorm(ncol(W),0,1)
      , sigmaAlpha=runif(1,2,5)
      , sigmaKappa=runif(1,2,5)
      , sigmaGamma=runif(ncol(W),2,5)
      ,.RNG.name="base::Wichmann-Hill"
      ,.RNG.seed=1971)
  )
}

# Read model
source ("~/Dropbox/CreditPreferences/Code/HLMprobitSelfReportedCreditSomeRandom.txt")

#######################################################
# Try one small run with one chain to test jags
#######################################################
# trialRun <- run.jags (
# 	model=HLMprobit,
# 	monitor=jags.parameters, n.chains=1,
# 	data=jags.data, inits=jags.inits(),
# 	adapt=10, thin=1, burnin=50, sample=100,
# 	summarise=TRUE, plots=FALSE )
#######################################################

# 
# Estimate model using two cores for two chains
system.time(
  out <- mclapply(1:2, function(x) {
    model.jags.re <- try(run.jags( model=HLMprobit
                                   , monitor=jags.parameters
                                   , inits=jags.inits()
                                   , n.chains=1
                                   , monitor.deviance=TRUE
                                   , data=jags.data
                                   #, thin=1, burnin=100, sample=100
                                   , thin=10, burnin=500, sample=250
                                   #, thin=240, burnin=60000, sample=120000
                                   , plots=FALSE
    ))
    if(inherits(model.jags.re,"try-error")) {return()}
    return(model.jags.re)
  }, mc.cores=2 )
)

# user   system  elapsed 
# 1558.406    5.050 1606.547 

save (out, file=paste0(getwd(), "/Code/Files4Server/giniMarket.RData"))


