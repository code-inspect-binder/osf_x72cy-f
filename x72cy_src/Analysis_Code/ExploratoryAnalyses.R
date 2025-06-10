#############################################################
##                Exploratory Analyses                     ##
##                                                         ##
##  0. Set Working Directory And Load Relevant R Packages  ##
#############################################################

  rm(list = ls())
  # setwd('...') # set working directory to 'Analysis_Code' folder
  source('ConfirmatoryAnalyses.R') # some results from confirmatory analyses are required (takes a while)

##########################
## 0.1. Global Settings ##
##########################
  
  saveFigures <- FALSE # save figures in separate files

#################################################################################################
## 1. Accuracy Rate of Experts from SSRP (Camerer et al., 2018) and ML2 (Forsell et al., 2018) ##
#################################################################################################
  alpha.omega <- 10
  beta.omega  <- 10
  
# For Camerer et al. (2018) SSRP
  dat0       <- read.csv('../Data/Expert_Data/SSRP_ExpertPredictions.csv'  , header = TRUE)
  dat1       <- read.csv('../Data/Expert_Data/SSRP_ReplicationOutcomes.csv', header = TRUE)
  
  dat0       <- dat0[dat0$treatment == 'm3', ]           # only include data from treatment 2 (all data)
  
  # extract relevant information
  dat2                    <- dat0[, c('study', 'pid', 'm3_b3')]
  colnames(dat2)[3]       <- 'beliefs' 
  dat2$beliefs            <- (1 - dat2$beliefs)
  # recode participant id
  dat2$pid                <- dat2$pid[, drop = TRUE]
  dat2$pid                <- factor(dat2$pid, levels = levels(dat2$pid), labels = 1:length(levels(dat2$pid)))
  dat2$pid                <- as.numeric(as.character(dat2$pid))
  
  # add replication belief and replication outcome
  dat2$replication.belief  <- ifelse(dat2$beliefs > 0.5, 1, 0)
  replication.outcome     <- dat1$rep_sr_rp
  
  dat2$replication.outcome <- NA
  for(i in 1:21){
    dat2$replication.outcome[dat2$study == i] <- replication.outcome[i]
  }
  # add variable that coded whether participant guessed correctly
  dat2$guessed.correctly <- ifelse(dat2$replication.outcome == dat2$replication.belief, 1, 0)
  
  # calculate group accuracy
  group.accuracy    <- aggregate(guessed.correctly ~ study, data = dat2, FUN = mean) 
  ssrp.prop.accuracy <- length(which(group.accuracy$guessed.correctly > .5))/21
  
  # Estimate hierarchical Omega For SSRP study
  # Prior distribution:   Omega ~ Beta(1, 1); Kappa ~ Gamma(0.01, 0.01)
  
  set.seed(4491)
  # generate MCMC chain 
  mcmc.chain <-  genMCMC(data          = dat2, 
                         sName         = 'pid',               # subject id
                         yName         = 'guessed.correctly', # responses '0' or '1
                         aOmega        = alpha.omega,
                         bOmega        = beta.omega,
                         numSavedSteps = 20000 , 
                         thinSteps     = 10,
                         setSeed       = 4491)
  # check convergence of chains
  posterior.omega.check <- cbind(mcmc.chain[[1]][,'omega'], mcmc.chain[[2]][,'omega'])
  nsamples              <- nrow(posterior.omega.check)
  plot( 1:nsamples, posterior.omega.check[,1], type = 'l', las = 1, bty = 'n')
  lines(1:nsamples, posterior.omega.check[,2], col = 'red')
  # extract posterior samples of omega 
  posterior.ssrp <- c(posterior.omega.check[, 1], posterior.omega.check[, 2])
  # save output
  # write.table(as.data.frame(posterior.ssrp), file = '../R_Output/Samples/posterior_omega_ssrp.txt', quote = FALSE, row.names = FALSE)
  
  # parameter estimation
  median.and.credible.interval <- quantile(posterior.ssrp, c(0.025, 0.5, 0.975))
  median.and.credible.interval
  
  rm(dat0, dat1, dat2, mcmc.chain, posterior.omega.check, median.and.credible.interval)
  
# For Forsell (2018) ML2
# In this dataset, each row is an expert and predictions are represented in columns. 
# explanation: Q-studynumber-A is prediction of replication on 0-100 scale. 
# Only 24 studies are included. The outcomes are recorded in the 'ML2_ReplicationOutcomes.csv' file. 

  dat0 <- read.csv('../Data/Expert_Data/ML2_ExpertPredictions.csv', header = TRUE, check.names = FALSE)
  dat1 <- read.csv2('../Data/Expert_Data/ML2_ReplicationOutcomes.csv', header = TRUE)
  
  colnames(dat0) <- substr(colnames(dat0), 1,4)
  
  dat2          <- dat0[grep('A', names(dat0))]              # select columns that contain 'A'
  dat2$IPAd     <- NULL                                      # remove this random column
  dat2$pid      <- 1:nrow(dat2)                               # change partipant ID to integer
  dat2          <- dat2 %>% gather(study, beliefs, Q1A.:Q28A) # convert to long format
  dat2$study.id <- rep(1:28, each = 91)                      # change study numbers
  exclude      <- c(1, 7, 19, 20)                           # specify studies that were excluded in Forsell paper (1, 7, 19, 20)
  dat2          <- dat2[which(!dat2$study.id %in% exclude),]   # remove selected studies
  dat2$study    <- rep(1:24, each = 91)                      # renumber studies
  
  # add replication belief and replication outcome
  dat2$replication.belief  <- ifelse(dat2$beliefs > 50, 1, 0)
  replication.outcome     <- dat1$Replication_Outcome
  
  dat2$replication.outcome <- NA
  for(i in 1:24){
    dat2$replication.outcome[dat2$study == i] <- replication.outcome[i]
  }
  # add variable that coded whether participant guessed correctly
  dat2$guessed.correctly <- ifelse(dat2$replication.outcome == dat2$replication.belief, 1, 0)
  
  # calcuta group accuracy
  group.accuracy     <- aggregate(guessed.correctly ~ study, data = dat2, FUN = mean) 
  ml2.prop.accuracy <- length(which(group.accuracy$guessed.correctly > .5))/24
  
  # Estimate hierarchical Omega For Forsell study
  # Prior distribution:   Omega ~ Beta(1, 1); Kappa ~ Gamma(0.01, 0.01)
  
  set.seed(4491)
  # generate MCMC chain 
  mcmc.chain <-  genMCMC(data          = dat2, 
                         sName         = 'pid',               # subject id
                         yName         = 'guessed.correctly', # responses '0' or '1
                         aOmega        = alpha.omega,
                         bOmega        = beta.omega,
                         numSavedSteps = 20000 , 
                         thinSteps     = 10,
                         setSeed       = 4491)
  # check convergence of chains
  posterior.omega.check <- cbind(mcmc.chain[[1]][,'omega'], mcmc.chain[[2]][,'omega'])
  nsamples              <- nrow(posterior.omega.check)
  plot( 1:nsamples, posterior.omega.check[,1], type = 'l', las = 1, bty = 'n')
  lines(1:nsamples, posterior.omega.check[,2], col = 'red')
  # extract posterior samples of omega 
  posterior.ml2 <- c(posterior.omega.check[, 1], posterior.omega.check[, 2])
  # save output
  # write.table(as.data.frame(posterior.ml2), file = '../R_Output/Samples/posterior_omega_ml2.txt', quote = FALSE, row.names = FALSE)
  
  # parameter estimation
  median.and.credible.interval <- quantile(posterior.ml2, c(0.025, 0.5, 0.975))
  median.and.credible.interval
  
  rm(dat0, dat1, dat2, mcmc.chain, posterior.omega.check, median.and.credible.interval)
  
## Figure: Accuracy Rates of Experts and Laypeople
  # dataframe with posterior samples of two expert and two non-expert groups
  df                   <- data.frame(posterior.omega.description.only, 
                                     posterior.omega.description.plus.evidence,
                                     posterior.ssrp, 
                                     posterior.ml2)
  df                   <- df %>% gather(condition, omega, 1:4)
  df$condition         <- as.factor(df$condition)
  levels(df$condition) <- c('Experts ML2', 'Laypeople (Description)',
                            'Laypeople (Evidence)', 'Experts SSRP')
  
  # set colors 
  figureColors[['cols']]             <- c('grey',figureColors[['col']][2], figureColors[['col']][8], 'grey') # order in figure: Experts ML2, Laypeople Description, Laypeople Evidence, Experts SSRP
  figureColors[['cols.transparent']] <- add.alpha(figureColors[['cols']], 0.7)

  # Save and export image as .eps (8.18 x 4.88 inch)
  if(saveFigures) cairo_ps(file = '../R_Output/Images/AccuracyPerProject_Revision.eps', onefile = TRUE, fallback_resolution = 600, 
           width = 8.18, height = 4.88)
  
  par(cex.main = 2.2, mar = c(5, 5.5, 2, 2) + 0.1, mgp = c(3, 1, 0), 
      cex.lab = 2, font.lab = 1, cex.axis = 1.6, bty = 'n', las = 1)
  plot(x = c(.5, .8), y = c(0, 50), type = 'n', 
       xlab = expression(paste("Accuracy Rate ", omega, " (in %)")), ylab = 'Density', 
       main = '', axes = FALSE) 
  axis(2, lwd.tick = 0, labels = FALSE)
  axis(1, at = seq(.5, .8, by = .05), labels = seq(50, 80, by = 5))
  d <- lapply(split(df$omega, df$condition), density)
  Map(function(dens, col) polygon(dens, col = col), 
      dens = d, col = figureColors[['cols.transparent']])
  text(.64,44.5,'Laypeople\n(Evidence)', cex = 1.3)
  text(.55,44.5,'Laypeople\n(Description)', cex = 1.3)
  text(.621,25,'Experts\nML2',  cex = 1.3)
  text(.75,25, 'Experts\nSSRP', cex = 1.3)
  
  graphics.off()
  # Image is saved as specified in cairo_ps
  
  par(mfrow=c(1,1)) # reset layout of plots to default
  
  rm(df, d)
  
###################################################################
## 2. Accuracy Rate Of The Bayes Factor Without Human Evaluation ##
###################################################################

  # This analysis will be conducted in JASP using a Bayesian binomial test
  #
  # corresponding JASP file: 'PredictionAccuracy_BayesFactorWithCutoff10.jasp'
  # datafile: 'DescriptionsAndStatisticsPerStudy.csv'
  #
  # The dependent variable is 'Bayes Factor Cutoff 10 Prediction Outcome' which encodes correct and false predictions
  # of the Bayes factor with a cutoff of 10.
  # That is, if the original study has a BF < 10, the prediction is 'will not replicate' 
  #          if the original study has a BF > 10, the prediction is 'will replicate' 
  # 
  # Prior destribution (parameter estimation): Beta(1, 1) distribution on the accuracy rate theta
  # Prior destribution (hypothesis testing)  : truncated Beta(1, 1)I[0.5, 1] distribution on the accuracy rate theta
  
###########################################
##  3. Signal Detection Theory Analysis  ##
###########################################
  
  # load data:
  dat        <- out.full 
  
  # specify hits, false alarms, misses and correct rejections 
  dat$h      <- ifelse(dat$replication.outcome == 1 & dat$replication.belief == 1, 1, 0) #hits
  dat$f      <- ifelse(dat$replication.outcome == 0 & dat$replication.belief == 1, 1, 0) #false alarms
  dat$MI     <- ifelse(dat$replication.outcome == 1 & dat$replication.belief == 0, 1, 0) #misses
  dat$CR     <- ifelse(dat$replication.outcome == 0 & dat$replication.belief == 0, 1, 0) #correct rejections 
  
  # create aggregated data set 
  dat.agg    <- aggregate(h  ~ sub.number + condition, data = dat, FUN = sum)            # subject number, condition, number of hits
  dat.agg$f  <- aggregate(f  ~ sub.number + condition, data = dat, FUN = sum)[,3]        # add false alarms 
  dat.agg$MI <- aggregate(MI ~ sub.number + condition, data = dat, FUN = sum)[,3]        # add misses
  dat.agg$CR <- aggregate(CR ~ sub.number + condition, data = dat, FUN = sum)[,3]        # add correct rejections 
  dat.agg$sr <- aggregate(replication.outcome ~ sub.number + condition, data = dat, FUN = sum)[,3]                    # total number of successful replications 
  dat.agg$fr <- aggregate(replication.outcome ~ sub.number + condition, data = dat, FUN = function(x) sum(x==0))[,3]  # total number of failed replications 
  
  # set settings for sampling 
  niter   <- 10000
  nburnin <- 1000
  
  set.seed(2019)
  for (dataset in 1:2){ #analyze both conditions
    if (dataset == 1)
      # the Description Only condition data
      data <- dat.agg[dat.agg$condition == 'DescriptionOnly',] 
    if (dataset == 2)
      # the Description Plus Bayes Factor condition data
      data <- dat.agg[dat.agg$condition == 'DescriptionPlusEvidence',] 
    
    h  <- data[,'h']
    f  <- data[,'f']
    MI <- data[,'MI']
    CR <- data[,'CR']
    sr <- data[,'sr']
    fr <- data[,'fr']
    k  <- nrow(data)	
    
    data <- list('h', 'f', 'sr', 'fr', 'k') # to be passed on to JAGS
    myinits <- list(
      list(deltac = rep(0, k), deltad = rep(0, k), xic = 0.5, xid = 0.5, muc = 0, mud = 0, lambdac = 1, lambdad = 1,
      # set seed
      .RNG.name = "base::Wichmann-Hill", 
      .RNG.seed = 4491))  
    
    # parameters to be monitored:	
    parameters <- c('muc', 'mud', 'sigmac', 'sigmad')
    
    if (dataset == 1) # Description Only
    {
      # The following command calls JAGS with specific options.
      # For a detailed description see the R2jags documentation.
      description.samples <- jags(data, inits=myinits, parameters,
                                  model.file ='SDT_JagsModel.txt',
                                  n.chains=1, n.iter=niter, n.burnin=nburnin, n.thin=1)				
    }
    
    if (dataset == 2) # Description Plus Evidence
    {
      # The following command calls JAGS with specific options.
      # For a detailed description see the R2jags documentation.
      evidence.samples <- jags(data, inits=myinits, parameters,
                               model.file ='SDT_JagsModel.txt',
                               n.chains=1, n.iter=niter, n.burnin=nburnin, n.thin=1)				
    }
  }		
  # Now the values for the monitored parameters are in the 'description.samples' and 
  # 'evidence.samples' objects, ready for inspection.
  
  # extract relevant parameters for sensitivity and bias per condition 
  description.mud   <- description.samples$BUGSoutput$sims.array[,,'mud']
  description.muc   <- description.samples$BUGSoutput$sims.array[,,'muc']
  evidence.mud      <- evidence.samples$BUGSoutput$sims.array[,,'mud']
  evidence.muc      <- evidence.samples$BUGSoutput$sims.array[,,'muc']
  
  # save output
  # write.table(as.data.frame(description.mud), file = '../R_Output/SDT output/sensitivity_descriptiononly.txt'    , quote = FALSE, row.names = FALSE)
  # write.table(as.data.frame(description.muc), file = '../R_Output/SDT output/bias_descriptiononly.txt'           , quote = FALSE, row.names = FALSE)
  # write.table(as.data.frame(evidence.mud)   , file = '../R_Output/SDT output/sensitivity_descriptionEvidence.txt', quote = FALSE, row.names = FALSE)
  # write.table(as.data.frame(evidence.muc)   , file = '../R_Output/SDT output/bias_descriptionEvidence.txt'       , quote = FALSE, row.names = FALSE)
  
  # median and 95% credible interval
  quantile(description.mud, c(0.5, 0.025, 0.975))
  quantile(evidence.mud   , c(0.5, 0.025, 0.975))
  quantile(description.muc, c(0.5, 0.025, 0.975))
  quantile(evidence.muc   , c(0.5, 0.025, 0.975))
  
  # function to calculate AUC (based on Wickens (2001), p.68, Equation 4.6)
  groupAUC.description <- pnorm(description.mud / sqrt(2)) # signal distribution has sigma = 1
  quantile(groupAUC.description, c(0.025, 0.5, 0.975))
  groupAUC.evidence <- pnorm(evidence.mud / sqrt(2)) # signal distribution has sigma = 1
  quantile(groupAUC.evidence, c(0.025, 0.5, 0.975))
  
  ## Figure 8 SDT Plot with Sensitivity and Bias per Condition ##
  
  #  R Code adopted from:
  #  Lee, M. D., & Wagenmakers, E.-J. (2013).Bayesian Cognitive Modeling: A Practical Course.Cambridge University Press
  #  Figure 11.5 (p. 163)
  
  # set indices of samples we want to plot 
  set.seed(2019)
  keepi <- 1000
  niter <- 10000
  keep  <- sample(niter, keepi)
  
  # densities for bias
  dens.description.muc <- density(description.muc)
  dens.evidence.muc    <- density(evidence.muc)
  # densities for sensitivity
  dens.description.mud <- density(description.mud)
  dens.evidence.mud    <- density(evidence.mud)
  
  # Save and export image as .pdf (8.18 x 6 inch)
  if(saveFigures) pdf(file = '../R_Output/Images/SDTposteriors_Revision.pdf', 
             width = 8.14, height = 6)
  
  layout(matrix(c(1,2,3,0), 2, 2, byrow=T), width=c(2/3, 1/3), heights=c(2/3,1/3))
  
  # scatter plot of samples
  par(mar=c(2,1,1,0))
  
  plot(  description.mud, description.muc, pch = '.', cex = 2, col = figureColors[['col.des.transparent']][1],  xlab = '', ylab='', axes = F, xlim = c(0,1.25), ylim = c(-.6,0))
  points(evidence.mud, evidence.muc      , pch = '.', cex = 2, col = figureColors[['col.evi.transparent']][1])
  box(lty = 1, col = 'black')
  
  # density for bias parameter
  par(mar=c(2,2,1,6))
  plot(dens.description.muc$y, dens.description.muc$x, xlim = rev(c(0,12)), type = 'n',
       axes = F, xlab = '', ylab = '', ylim = c(-.6,0))
  axis(4, at = c(-.6, -.4, -.2, 0), cex.axis = 1.5)
  mtext('Bias', side = 4, line = 3.5, cex = 1.5)
  # shade area under the curve
  polygon(dens.description.muc$y, dens.description.muc$x, col = figureColors[['col.des.transparent']][1])
  polygon(dens.evidence.muc$y   , dens.evidence.muc$x   , col = figureColors[['col.evi.transparent']][1])
  
  # density for sensitivity parameter
  par(mar=c(6,1,0,0))
  plot(density(description.mud), zero.line = F, main = '', ylab = '', xlab = '', cex.lab = 1.3,
       axes = F, xlim = c(0,1.25), ylim = c(0,8), lwd = 3, type = 'n')
  axis(1, at = c(0, .25, .5, .75, 1, 1.25), cex.axis = 1.5)
  mtext('Discriminability', side = 1, line = 3.5, cex = 1.5)
  # shade area under the curve
  polygon(dens.description.mud$x, dens.description.mud$y, col = figureColors[['col.des.transparent']])
  polygon(dens.evidence.mud$x   , dens.evidence.mud$y, col = figureColors[['col.evi.transparent']])
  
  graphics.off()
  # Image is saved as specified in cairo_ps
  
  par(mfrow=c(1,1)) # reset layout of plots to default
  
  ## Figure 9a : Signal and Noise distribution per condition
  
  # set condition (select one at a time)
  condition <- 'DescriptionOnly'
  # condition <- 'DescriptionPlusEvidence'
  
  if (condition == 'DescriptionOnly') {
    mud <- description.mud
    muc <- description.muc
    col <- figureColors[['col.des.transparent']]
  }
  if (condition == 'DescriptionPlusEvidence') {
    mud <- evidence.mud
    muc <- evidence.muc
    col <- figureColors[['col.evi.transparent']]
  }
  
  # get relevant information 
  # compute mode of signal distribution
  estimate_mode <- function(x) {
    d <- density(x)
    d$x[which.max(d$y)]
  }
  
  modeD           <- estimate_mode(mud)
  medianC         <- median(muc)
  empiricalLambda <- (0.5*modeD) + medianC
  
  # get height of noise distribution at lambda
  heightLambda           <- dnorm(empiricalLambda, 0,1)
  # get height of noise distribution at optimal criterion d/2
  heightOptimalCriterion <- dnorm(0.5*modeD, 0,1)
  
  densityNoise  <- density(rnorm(5e6))
  densitySignal <- density(rnorm(5e6, modeD, 1))
  
  # Save and export image as .eps (6 x 6 inch)
  if(saveFigures) cairo_ps(file = paste0('../R_Output/Images/SignalNoise_',condition,'_Revision.eps'), onefile = TRUE, fallback_resolution = 600,
           width = 6, height = 4)
  par(mfrow=c(1,1), mar=c(2.5, 0.5, 0, 0.5), mgp=c(2.4, 0.8, 0))
  plot(densityNoise, zero.line = FALSE, main = '', ylab = '', xlab = '', cex.lab = 1.3,
       axes = FALSE, xlim = c(-3.8, 4.5), ylim = c(0, 0.5), lwd = 3, type = 'n')
  mtext('Feeling of Replicability', side = 1, line = 1, cex = 1.5)
  # shade area under the curve
  polygon(densitySignal$x, densitySignal$y, col = col, lwd = 1.5)
  polygon(densityNoise$x , densityNoise$y, lwd = 1.5)
  # include optimal criterion
  lines(x = c(0.5 * modeD, 0.5 * modeD)        , y = c(0, heightOptimalCriterion), lty = 3, lwd = 1.5)
  lines(x = c(empiricalLambda, empiricalLambda), y = c(0, heightLambda)          , lty = 2, lwd = 1.5)
  # add text criteria
  arrows(x0 = -1, x1 = empiricalLambda-0.05, y0 = 0.45, y1 = (heightLambda-0.03), length = .1)
  text(x = -2.2, y = 0.47, 'adopted criterion', cex = 1.2)
  arrows(x0 = (0.5 * modeD), x1 = (0.5 * modeD), y0 = 0.45, y1 = (heightOptimalCriterion+.01), length = .1)
  text(x = 1.45,    y = 0.47, 'optimal criterion', cex = 1.2)
  # bias
  arrows(x0 = empiricalLambda, x1 = (0.5*modeD), y0 = .1, y1 = .1, length = .05, code = 3)
  text(x = 0.5*modeD+.4, y = .1, 'bias', cex = 1.2)
  
  graphics.off()
  # Image is saved as specified in cairo_ps
  
  
  ## Figure 9b: ROC plots per condition ##
  
  #  R Code adopted from:
  #   Selker, R., van den Bergh, D., Criss, A. H., & Wagenmakers, E. J. (2019). Parsimonious estimation of signal 
  #   detection models from confidence ratings. Behavior Research Methods, 1-15. Retrieved from https://osf.io/v3b76/
  
  # set condition (select one at a time)
  # condition <- 'DescriptionOnly'
  condition <- 'DescriptionPlusEvidence'
  
  if (condition == 'DescriptionOnly') {
    mud <- description.mud
    muc <- description.muc
    col <- figureColors[['col.des']]
    col.trans <- figureColors[['col.des.transparent']]
  }
  if (condition == 'DescriptionPlusEvidence') {
    mud <- evidence.mud
    muc <- evidence.muc
    col <- figureColors[['col.evi']]
    col.trans <- figureColors[['col.evi.transparent']]
  }
  
  # get relevant information and set points for ROC curve
  medianD  <- quantile(mud, 0.5)
  nPoints  <- 1e2
  lambdas  <- seq(-3.8, 3.8, length.out = nPoints)
  c        <- lambdas - (0.5 * medianD)
  
  # proportion of false alarms for each sample (based on median mu_d)
  pf <- pnorm((-0.5*medianD) - c)
  # proportion of hits for each sample (based on median mu_d)
  ph <- pnorm((0.5*medianD) - c)
  
  # calculates the quantiles directly (without storing the ROC intermediately)
  nShades  <- 50
  probs    <- c(seq(0.025, .5, length.out = nShades + 1), seq(.5, 0.975, length.out = nShades + 1)[-1])
  midPoint <- nShades+1
  
  qph <- matrix(NA, nrow = nPoints, ncol = length(probs))
  # compute CIs for each value of lambda
  probD      <- quantile(mud, probs = probs)
  for (i in seq_len(nPoints)){
    c        <- lambdas[i] - (0.5 * probD)
    qph[i, ] <- pnorm((0.5 * probD) - c)
  }
  
  # dataframe for ROC line
  df2plot <- data.frame(
    x = pf,
    y = qph[, midPoint]
  )
  
  # dataframe for shading
  df2shade <- data.frame(
    ylower = c(qph[, 1:nShades]),
    yupper = c(qph[, nShades:1 + midPoint]),
    x      = pf,
    shade  = factor(rep(1:nShades, each = length(pf)))
  )
  
  # make a color gradient
  rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
  intensity <- dbeta(seq(0, .5, length.out = nShades), 2.5, 2.5)
  intensity <- intensity / max(intensity)
  ourCols <- apply(colorRamp(c(col[2], col[1]))(intensity), 1, rgb2hex)
  names(ourCols) <- levels(df2shade$shade)
  
  fontsize <- 24
  ROCgraph <-
    ggplot(data = df2plot, aes(x = x, y = y)) +
    geom_ribbon(inherit.aes = FALSE, data = df2shade, show.legend = FALSE,
                aes(x = x, ymin = ylower, ymax = yupper, group = shade, fill = shade),
                alpha = .7) +
    scale_fill_manual(values = ourCols) +
    geom_line(size = 2, lineend = 'round', alpha = .5, color = 'black') +
    geom_abline(linetype = 'dashed') +
    scale_x_continuous('False Alarm Rate', breaks = 0:1) +
    scale_y_continuous('Hit Rate', breaks = 0:1) +
    geom_segment(inherit.aes = FALSE, x = 0, y = -Inf, xend = 1, yend = -Inf) +
    geom_segment(inherit.aes = FALSE, x = -Inf, y = 0, xend = -Inf, yend = 1) +
    theme(
      panel.grid       = element_blank(),
      panel.background = element_rect(colour = 'white', fill = 'white'),
      axis.text        = element_text(size = fontsize),
      axis.title       = element_text(size = fontsize),
      axis.line        = element_blank()
    )
  
  # function to calculate AUC for one sample ( based on Wickens (2001), p.68 (equation 4.6))
  grpAUCsamps <- pnorm(mud / sqrt(2)) # signal distribution has sigma = 1
  dd          <- density(grpAUCsamps)
  dfAUC       <- data.frame(x = dd$x, y = dd$y)
  
  # density plot of AUC samples
  fontsize <- 18
  AUCgraph <- ggplot(data = dfAUC, aes(x = x, y = y)) +
    geom_line(size = 1) +
    geom_polygon(fill = col.trans[1]) +
    scale_y_continuous(name = '', expand=c(0,0), limits = c(0, max(dd$y) + 1)) +
    scale_x_continuous('AUC', breaks = seq(.5, 1, .1), limits = c(.5, 1)) +
    theme_void() +
    theme(
      axis.title   = element_text(size = fontsize),
      axis.line.x  = element_line(),
      axis.ticks.x = element_line(),
      axis.text.x  = element_text(size = fontsize, angle = 0, debug = FALSE)
    )
  
  # overlay plots
  fullGraph <- ggdraw(ROCgraph) +
    draw_plot(AUCgraph, x = .2, y = -0.15, scale = .4)
  
  # Save and export image as .eps (6 x 6 inch)
  if(saveFigures) cairo_ps(file = paste0('../R_Output/Images/ROC_',condition,'_Revision.eps'), onefile = TRUE, fallback_resolution = 600,
           width = 6, height = 6)
  fullGraph
  graphics.off() # Image is saved as specified in cairo_ps
  
##################################################
##  4. Logistic Regression With Random Effects  ##
##################################################
  
  ## Generate models 
  # Priors in the model (not for testing):
  # weakly informative prior --> student_t(3,0,2.5) on the intercept and coefficient for condition 
   Prior <- c(set_prior('student_t(3,0,2.5)', class = 'b'),
             set_prior('student_t(3,0,10)',  class = 'sd'))
  
  # logistic regression model with random intercepts for subject and study and random effects of condition for study
  # note: takes quite some time! 
  m.studies <- brm(formula = guessed.correctly ~ 0 + intercept + condition +
                     (0 + intercept + condition | study) + (0 + intercept | subject), 
                   data = dat, family = bernoulli(link = "logit"),
                   prior = Prior,
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   seed = 2020, control = list(adapt_delta = 0.9))
  
  plot(m.studies)
  
  ## posteriors of coefficients
  quantile(posterior_samples(m.studies, 'b_intercept')[,1], c(0.025,0.5,0.975))
  quantile(posterior_samples(m.studies, 'b_conditionDescriptionPlusEvidence')[,1], c(0.025,0.5,0.975))
  
  ## get the posterior distributions of the intercepts, transformed to rate scale
  post.description <- inv_logit_scaled(posterior_samples(m.studies, 'b_intercept'))[,1]
  post.evidence    <- inv_logit_scaled(posterior_samples(m.studies, 'b_intercept') + posterior_samples(m.studies, 'b_conditionDescriptionPlusEvidence'))[,1]
  
  # parameter estimation
  median.and.credible.interval <- quantile(post.description, c(0.025, 0.5, 0.975))
  median.and.credible.interval
  
  # parameter estimation
  median.and.credible.interval <- quantile(post.evidence, c(0.025, 0.5, 0.975))
  median.and.credible.interval
  
  ## Figure S1: Accuracy Per Condition For Logistic Regression Model ##
  
  ## get the posterior distributions of the intercepts, transformed to rate scale
  post.description <- inv_logit_scaled(posterior_samples(m.studies, 'b_intercept'))[,1]
  post.evidence    <- inv_logit_scaled(posterior_samples(m.studies, 'b_intercept') + posterior_samples(m.studies, 'b_conditionDescriptionPlusEvidence'))[,1]
  
  # plot prior and posterior
  # Save and export image as .eps (8.5 x 5 inch)
  if(saveFigures) cairo_ps(file = '../R_Output/Images/AccuracyRatesRandomEffects_Revision.eps', onefile = TRUE, fallback_resolution = 600, 
           width = 8.5, height = 5)
  par(cex.main = 2.2, mar = c(5, 5.5, 2, 2) + 0.1, mgp = c(3, 1, 0), 
      cex.lab = 2, font.lab = 1, cex.axis = 1.6, bty = 'n', las = 1)
  plot(density(post.description), main = '', xlab = expression(paste("Accuracy Rate (in %)")), 
       ylab = 'Density', bty = 'n', las = 1,
       ylim = c(0,9.5), xlim = c(0,1), lwd = 2, axes = F)
  polygon(density(post.description), col = figureColors[['col.des.transparent']][1])
  lines(density(post.evidence), lwd = 2)
  polygon(density(post.evidence), col = figureColors[['col.evi.transparent']][1])
  axis(1)
  axis(2, labels = FALSE, lwd.ticks = 0, at = c(0,9.5))
  graphics.off()
  
  ## Figure S2: Random Effects Per Study and Condition ##  

  # prepare estimates for figure with random effects per condition
  Replications <- replication.results$Study
  newdata      <- data.frame(study = rep(1:27,2), condition = rep(c('DescriptionOnly','DescriptionPlusEvidence'), each=27))
  
  # get fitted accurary rate per study 
  fits <- fitted(m.studies, newdata = newdata, re_formula = ~ 0 + intercept + condition + (0 + intercept + condition | study)) %>%
    as_tibble() %>%
    bind_cols(newdata) %>% #add 'condition' and 'study'
    mutate(study = rep(Replications, 2)) #change to study labels 
  # add the average accuracy rate per condition 
  fit <- fitted(m.studies, newdata = data.frame(condition = c('DescriptionOnly','DescriptionPlusEvidence')), re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(study = rep('Average',2), condition = as.factor(c('DescriptionOnly','DescriptionPlusEvidence'))) %>%
    bind_rows(fits)
  # order based on average accuracy rate (over conditions)
  fit$order <- reorder(as.factor(fit$study), fit$Estimate, mean) 
  fit$order <- relevel(fit$order, ref = 'Average')
  
  # colours
  cols <- c('grey38',figureColors[['col.des']][1], figureColors[['col.evi']][1])
  # add the observed average accuracy per study and per condition 
  fit$obs.accuracy <- c(aggregate(dat$guessed.correctly, list(condition = dat$condition), mean)[,'x'],
                        aggregate(dat$guessed.correctly, list(study = dat$study, condition = dat$condition), mean)[,'x']) 
  
  # Save and export image as .eps (8.5 x 5 inch)
  if(saveFigures) cairo_ps(file = '../R_Output/Images/RandomEffect_Revision.eps', onefile = TRUE, fallback_resolution = 600, 
           width = 8, height = 10)
  fit %>%
    ggplot(aes(y = order, x = Estimate, xmin = Q2.5, xmax = Q97.5, colour = condition)) +
    geom_pointrangeh(position = position_dodgev(height = .5), size = .8, alpha = .7) +
    geom_point(data = subset(fit, study == 'Average'), stat = 'identity', size = 6, position = position_dodgev(height = .5)) +
    geom_point(aes(y = order, x = obs.accuracy), stat = 'identity', shape = 4, size = 3, position = position_dodge2v(height = .5), colour = 'black') +
    theme_classic(base_size = 19) + 
    theme(legend.position = 'none') +
    ylab('') +
    xlab('Accuracy Rate') + 
    scale_color_manual(values = cols) +
    geom_segment(aes(x=0,xend=1,y=-Inf,yend=-Inf, colour = 'black')) +
    geom_segment(aes(y=1, yend=28, x=-Inf, xend=-Inf, colour = 'black')) +
    theme(axis.line=element_blank()) + 
    theme(plot.margin = unit(c(0.3,0.3,0.3,-0.3), 'cm'))
  graphics.off()
  