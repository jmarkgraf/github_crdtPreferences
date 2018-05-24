# Run Analysis:
# w/ clean dataset: prepared in line with Thewissen & Rueda
# w/ pure dataset: as few additional variables as possible

# Date: Feb 16, 2018
# Author: Jonas Markgraf
#########################

rm(list = ls())

# load libraries
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
library(interplot)  # to plot interaction effects


# set wd
possibles <- c("~/Dropbox/CreditPreferences/")
set_valid_wd(possibles)

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

allRisk <- c("rti2", "our", "relskillspec", "Oeschroutine")
offShoreRisk <- c("offshwalt","offsh")
omitObs <- as.numeric (attr (na.omit (complete.ess[,c(allRisk, offShoreRisk)])
                             , "na.action"))

allRisk.PC <- princomp(complete.ess[-omitObs,allRisk])
offShoreRisk.PC <- princomp(complete.ess[-omitObs,offShoreRisk])
round (cor (cbind (allRisk.PC$scores[,1], offShoreRisk.PC$scores[,1], complete.ess[-omitObs,c(allRisk,offShoreRisk)])), 2)

# Including all variables
Risk.PC  <- princomp(complete.ess[-omitObs, c(allRisk,offShoreRisk)])
round (cor (cbind (Risk.PC$scores[,1], Risk.PC$scores[,2], complete.ess[-omitObs,c(allRisk,offShoreRisk)])), 2)

# Excluding offshwalt
Risk.PC  <- princomp(complete.ess[-omitObs, c(allRisk,"offsh")])
round (cor (cbind (Risk.PC$scores[,1], Risk.PC$scores[,2], complete.ess[-omitObs,c(allRisk,offShoreRisk)])), 2)

#### MCMC mixed factor analysis ###
# Turn ordinal variables into ordered factors
complete.ess$Oeschroutine.fac <- as.ordered (complete.ess$Oeschroutine)
complete.ess$offshwalt.fac <- as.ordered (cut(complete.ess$offshwalt, c(-5,20,40,60,80,100)))
levels (complete.ess$offshwalt.fac) <- c("1","2","3","4","5")
complete.ess$offshwalt.fac.num <- as.numeric (complete.ess$offshwalt.fac)


#### FACTOR ANALYSIS ####
# Reduce dataset even further, to exclude observations with "wrong income"
# The main risk measure that we use is a factor that excludes OUR
tmp <- complete.ess[complete.ess$wrongincome != 1,]
AllFactor <- factanal (~rti2+relskillspec+offsh+Oeschroutine+offshwalt.fac.num 
                       , factors=1
                       , rotation="varimax"
                       , scores="regression"
                       , na.action=na.exclude
                       , data=tmp)


AllFactor$loadings
factor.scores <- AllFactor$scores

# Notation follows factanal help
# Useful math help from https://stats.stackexchange.com/questions/126885/methods-to-compute-factor-scores-and-what-is-the-score-coefficient-matrix-in
Lambda <- AllFactor$loadings
Phi    <- diag(AllFactor$uniquenesses)
Sigma <- Lambda %*% t(Lambda) + Phi
predictors <- tmp[,c("rti2","relskillspec","offsh","Oeschroutine","offshwalt.fac.num")]
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
patched.std.predictors <- apply (std.predictors, c(1,2), function (x) ifelse (is.na(x), 0, x))
patched.std.predictors.valid <- matrix (NA, ncol=ncol(patched.std.predictors), nrow=nrow(patched.std.predictors))
for (i in 1:nrow(patched.std.predictors)){
   patched.std.predictors.valid[i,] <- ifelse (allMissing[i]==1, rep(NA, ncol(patched.std.predictors)), patched.std.predictors[i,])   
}

thompson.scores.all <- patched.std.predictors.valid %*% solve (corrMatrix) %*% Lambda  

tmp$risk <- thompson.scores.all

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
Y <- as.matrix (cbind(tmp$rti2, tmp$relskillspec, tmp$offsh, tmp$Oeschroutine, as.numeric(tmp$offshwalt.fac)))
complete.missing <- apply (Y, 1, invalid)

Y <- Y[complete.missing==FALSE,]
Y[,1:3] <- apply (Y[,1:3], 2, function (x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
N <- nrow(Y)   # Number of observations
J <- ncol(Y)   # Number of risk items
L <- 3 # number of continuous variables
intercept <- rep(1,N)

source ("~/Dropbox/CreditPreferences/Code/mixFactAnal.R")

# Priors
jags.parameters <- c("lambda","deviance","factor.score","kappa","sigma","alpha")
Data=list(Y=Y, J=J, N=N, L=L
          , intercept=intercept)
jags.data <- dump.format(Data) #
jags.inits <- function()
{
   dump.format(
      list(
         lambda=runif(5,0,2)
         , alpha=rnorm(4)  # To coincide with negative priors
         , kappa.unsorted=rnorm(3)
         #, beta=rnorm(1,0,2)
      ) )
}

# Run model (so far, failure to define node lambda)
jags.model <- run.jags( model=mixFactAnal
                        , monitor=jags.parameters
                        , inits=jags.inits()#list(jags.inits(), jags.inits())
                        , n.chains=1
                        , method="parallel"
                        , data=jags.data
                        , adapt=1#000
                        , thin=2#0
                        , burnin=20#000
                        , sample=50#0
                        #, thin=1, burnin=10, sample=50
                        , summarise=FALSE
                        , plots=FALSE )
###########################################################################


fit.baseline <- lmer(gincdif2 ~ brwmny
  + income + male + agea + unemplindiv 
  + eduyrs2 + mbtru2 + rlgdgr 
  + socgdp + gdpc
  + (1 + brwmny | cntry.yr),
  weights = dweight, data=complete.ess)
summary(fit.baseline)



fit.risk <- lmer(gincdif2 ~ brwmny*risk
                     + log(income) + male + agea + unemplindiv 
                     + eduyrs2 + mbtru2 + rlgdgr 
                     + socgdp + log(gdpc)
                     + (1 + brwmny*risk | cntry.yr),
                     weights = dweight, data=tmp)
summary(fit.risk)

interplot (fit.risk, "brwmny", "risk")

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


