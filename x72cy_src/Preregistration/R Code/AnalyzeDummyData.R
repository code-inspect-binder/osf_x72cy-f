###########################################################
##               Analyze Dummy Data                      ##
##                                                       ##
## 0. Set Working Directory And Load Relevant R Packages ##
###########################################################


rm(list = ls())
# setwd("...")

library(polspline)                 # to fit posterior density
source("samplingFunctions.R")      # to estimate Bayesian Spearman's rho & Omega

dat <- read.csv('Data/DummyData_Full.csv', header = TRUE)
# create variable 'guessed correctly'; 1 == 'correct', 0 = 'incorrect'
dat$guessed.correctly <- ifelse(dat$replication.belief == dat$replication.outcomes, 1, 0)
nstudies              <- length(unique(dat$study))
# create variable 'confidence.scale.wide' for the spearman analysis range [-1, 1]
# if participants think a study replicates, have positive confidence.ratings
# if participants think a study does not replicate, have negative confidence.ratings

#########################################
## 1. Quality Check of Data And Design ##
#########################################


# Research question: Do lay people perform worse than guessing level? 
# He: Omega is not smaller than .5
# Hr: Omega is smaller than .5
# The quality check fails, if we collect strong evidence (BF >= 10) in favor for Hr
#
# Prior distribution: Omega ~ Beta(1, 1); Kappa ~ Gamma(0.01, 0.01)
alpha.omega <- 1
beta.omega  <- 1

# generate MCMC chain 
mcmc.chain <-  genMCMC(data          = dat, 
                       sName         = "id",                # subject id
                       yName         = "guessed.correctly", # responses '0' or '1
                       aOmega        = alpha.omega,
                       bOmega        = beta.omega,
                       numSavedSteps = 20000 , 
                       thinSteps     = 10 )
# check convergence of chains
posterior.omega.check <- cbind(mcmc.chain[[1]][,'omega'], mcmc.chain[[2]][,'omega'])
nsamples              <- nrow(posterior.omega.check)
plot( 1:nsamples, posterior.omega.check[,1], type = 'l', las = 1, bty = 'n')
lines(1:nsamples, posterior.omega.check[,2], col = 'red')
# extract posterior samples of omega 
posterior.omega <- c(posterior.omega.check[, 1], posterior.omega.check[, 2])

# parameter estimation
median.and.credible.interval <- quantile(posterior.omega, c(0.025, 0.5, 0.975))
median.and.credible.interval
hist(posterior.omega, las = 1)

# hypothesis testing: encompassing prior approach (Klugkist et al., 2005)
Iprior <- 0.5 # Iprior = 50% since prior distribution is centered around 0.5
Ipost  <- sum(posterior.omega < 0.5)/length(posterior.omega)
BFre   <- Ipost/Iprior; BFre

rm(posterior.omega.check, posterior.omega)

########################################
## 2.1 Test Differences in Conditions ##
########################################

# This analysis will be conducted in JASP using a Bayesian independent samples t-test
# with grouping variable 'condition' and dependent variable 'log.brier.score'
# 
# Prior destribution: Default cauchy prior (width 0.707) on the effect size delta

################################
## 2.2 Estimate Accuracy Rate ##
################################


# Research question: Is the accuracy rate higher than chance?
# H0: Omega is exactly 0.5
# Hr: Omega is higher than .5
# If the date suggest strong evidence (BF > 10) in favor for the hypothesis that the Brier scores differ across 
# conditions we will analyze the accuracy rate for each condition separately; otherwise, we use the full dataset
#
# Prior distribution: Omega ~ Beta(10, 10); Kappa ~ Gamma(0.01, 0.01)
alpha.omega <- 10
beta.omega  <- 10

## Condition 'DescriptionOnly' ##

dat.description.only         <- dat[dat$condition == 'DescriptionOnly', ]
nsubjects                    <- length(unique(dat.description.only$id))
dat.description.only$subject <- rep(1:nsubjects, each = nstudies)
# generate MCMC chain 
mcmc.description.only <-  genMCMC(data  =  dat.description.only, 
                                  sName = "subject" ,           # subject id
                                  yName = "guessed.correctly" , # responses '0' or '1
                                  aOmega        = alpha.omega,
                                  bOmega        = beta.omega,
                                  numSavedSteps = 20000 , 
                                  thinSteps     = 10 )
# convergence check
posterior.omega.check <- cbind(mcmc.description.only[[1]][,'omega'], mcmc.description.only[[2]][,'omega'])
nsamples              <- nrow(posterior.omega.check)
plot( 1:nsamples, posterior.omega.check[,1], type = 'l', las = 1, bty = 'n')
lines(1:nsamples, posterior.omega.check[,2], col = 'red')
# extract posterior samples of omega 
posterior.omega.description.only <- c(posterior.omega.check[, 1], posterior.omega.check[, 2])
# hypothesis testing: Savage-Dickey density ratio
posterior.density           <- logspline(posterior.omega.description.only)
density.zero.point          <- dlogspline(    0.5, posterior.density)
correction.factor.posterior <- 1 - plogspline(0.5, posterior.density)
prior.density.zero.point    <- dbeta(0.5, alpha.omega, beta.omega)
correction.factor.prior     <- pbeta(0.5, alpha.omega, beta.omega, lower.tail = FALSE)
# Bayes factor
BFr0.description.only       <- (prior.density.zero.point / correction.factor.prior)/
                               (density.zero.point / correction.factor.posterior) 
# parameter estimation
median.and.credible.interval <- quantile(posterior.omega.description.only, c(0.025, 0.5, 0.975))
median.and.credible.interval
# plot prior and posterior
plotPosterior(posterior.omega.description.only, 
              posteriorDensity = posterior.density, 
              xlabName = 'Omega', xLow = 0, xHigh = 1, yHigh = 25)
# add prior distribution
curve(dbeta(x, alpha.omega, beta.omega), 0, 1, lwd = 2, lty = 3, add = TRUE)
# remove variables
rm(posterior.omega.check, posterior.omega.description.only,
   posterior.density, density.zero.point, correction.factor.posterior, correction.factor.prior)

## Condition 'DescriptionPlusStatistics' ##

dat.description.plus.stats         <- dat[dat$condition == 'DescriptionPlusStatistics', ]
nsubjects                          <- length(unique(dat.description.plus.stats$id))
dat.description.plus.stats$subject <- rep(1:nsubjects, each = nstudies)
# generate MCMC chain 
mcmc.description.plus.stats <-  genMCMC(data  =  dat.description.plus.stats, 
                                        sName = "subject" ,           # subject id
                                        yName = "guessed.correctly" , # responses '0' or '1
                                        aOmega        = alpha.omega,
                                        bOmega        = beta.omega,
                                        numSavedSteps = 20000 , 
                                        thinSteps     = 10 )
# convergence check
posterior.omega.check <- cbind(mcmc.description.plus.stats[[1]][,'omega'], mcmc.description.plus.stats[[2]][,'omega'])
nsamples              <- nrow(posterior.omega.check)
plot( 1:nsamples, posterior.omega.check[,1], type = 'l', las = 1, bty = 'n')
lines(1:nsamples, posterior.omega.check[,2], col = 'red')
# extract posterior samples of omega 
posterior.omega.description.plus.stats <- c(posterior.omega.check[, 1], posterior.omega.check[, 2])
# hypothesis testing: Savage-Dickey density ratio
posterior.density           <- logspline(posterior.omega.description.plus.stats)
density.zero.point          <- dlogspline(   0.5, posterior.density)
correction.factor.posterior <- 1 - plogspline(0.5, posterior.density)
prior.density.zero.point    <- dbeta(0.5, alpha.omega, beta.omega)
correction.factor.prior     <- pbeta(0.5, alpha.omega, beta.omega, lower.tail = FALSE)
# Bayes factor
BFr0.description.only       <- (prior.density.zero.point / correction.factor.prior)/
                               (density.zero.point / correction.factor.posterior) 
# parameter estimation
median.and.credible.interval <- quantile(posterior.omega.description.plus.stats, c(0.025, 0.5, 0.975))
median.and.credible.interval
# plot prior and posterior
plotPosterior(posterior.omega.description.plus.stats, 
              posteriorDensity = posterior.density, 
              xlabName = 'Omega', xLow = 0, xHigh = 1, yHigh = 25)
# add prior distribution
curve(dbeta(x, alpha.omega, beta.omega), 0, 1, lwd = 2, lty = 3, add = TRUE)

#########################################
## 3. Estimate Bayesian Spearman's Rho ##
#########################################
  
# Research question: Can lay people capture the evidential value of the study?
# H0: Rho is 0 
# Hr: Rho is bigger than 0 
# If the data suggest strong evidence (BF > 10) in favor for the hypothesis that the Brier scores differ across 
# conditions we will analyze the correlation between (averaged) confidence ratings and replication effect sizes  
# for each condition separately; otherwise, we use the full dataset
#
# Prior distribution: Rho ~ Uniform(0, 1) for testing; Rho ~ Uniform(-1, 1) for estimation
  
confidence.rating       <- aggregate(confidence.rating ~ study, data = dat, FUN = mean)$confidence.rating # (range: 0,1)
replication.effectsizes <- aggregate(replication.effectsizes ~ study, data = dat, FUN = mean)$replication.effectsizes # (range: -1,1)
# generate MCMC chain (assuming wide uniform Beta distribution)
output.gibbs.sampler <- spearmanGibbsSampler(confidence.rating, replication.effectsizes, progBar = TRUE) 
posterior.rho        <- output.gibbs.sampler$rhoSamples
# hypothesis testing: Savage-Dickey density ratio
posterior.density           <- logspline(posterior.rho)
density.zero.point          <- dlogspline(   0, posterior.density)
correction.factor.posterior <- 1 -plogspline(0, posterior.density)
prior.density.zero.point    <- dbeta(0, 1, 1)
# Bayes factor
BFr0.rho                    <- prior.density.zero.point/(density.zero.point / correction.factor.posterior)
# parameter estimation
median.and.credible.interval <- quantile(posterior.rho, c(0.025, 0.5, 0.975))
median.and.credible.interval
# plot prior and posterior
plotPosterior(posterior.rho, 
              posteriorDensity = posterior.density, 
              xlabName = 'Spearmans Rho', xLow = -1, xHigh = 1, yHigh = 4)
# add prior distribution
curve(dunif(x, -1, 1), -1, 1, lwd = 2, lty = 3, add = TRUE)
# plot data
plot(confidence.rating, replication.effectsizes, xlim = c(0, 1), ylim = c(-0.2, 1), pch = 21, bg = 'grey', cex = 2)


