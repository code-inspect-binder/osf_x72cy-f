#############################################################
##               Confirmatory Analyses                     ##
##                                                         ##
##  0. Set Working Directory And Load Relevant R Packages  ##
#############################################################

  rm(list = ls())
  # setwd('...') # set working directory to 'Analysis_Code' folder
  source('loadPackages.R')
  source('PreprocessingQualtricsDataAnonymized.R') # preprocess the data

  # load data and create all relevant datasets for analyses
  dat                     <- out.full
  dat$confidence.rating   <- (dat$confidence.rating -.5) * 2 # transform scale of the confidence ratings from 0,1 to -1,1 
  
##########################
## 0.1. Global Settings ##
##########################
  
  saveFigures <- FALSE # save figures in separate files
  
  # load sampling functions
  source('samplingFunctions.R') # to estimate Bayesian Spearman's rho & Omega
  
  # load replication info
  replication.results       <- na.omit(read.csv(file = '../Data/Processed_Data/DescriptionAndStatisticsPerStudy.csv',stringsAsFactors = FALSE))[1:27,] # get the results for the 27 studies 
  replication.effectsizes   <- aggregate(replication.effectsize ~ study, data = dat, FUN = mean)$replication.effectsize # (range: -1,1)
  replication.outcomes      <- as.factor(replication.results$RS.outcome) # 1 = replicated, 0 = did not replicate
  replication.results$Study <- as.character(replication.results$Study)
  # Add study name to full data 
  for (i in 1:27){
    dat$study.label[dat$study == i] <- replication.results$Study[replication.results$Study.number == i]
  }
  dat$study.confidence_order <- reorder(dat$study.label, dat$confidence.rating, mean)
  dat$rep.outcome.label      <- ifelse(dat$replication.outcome == 0, 'did not replicate','replicated')
  
  # Select data for each condition 
  dat.description       <- dat[dat$condition == 'DescriptionOnly', ]
  dat.evidence          <- dat[dat$condition == 'DescriptionPlusEvidence', ]
  # order according to average confidence in the Description Only condition, and also use this as order in Evidence condition
  dat.description$study.confidence_order <- reorder(dat.description$study.label, dat.description$confidence.rating, mean)
  dat.evidence$study.confidence_order    <- reorder(dat.evidence$study.label, dat.evidence$confidence.rating, mean)
  
  # Colors for Figures
  figureColors <- list()
  figureColors[["col"]]     <- brewer.pal(n = 8, name = 'Paired')
  figureColors[["col.des"]] <- c(figureColors[["col"]][2], figureColors[["col"]][1]) # Descriptive Only:          dark blue and light blue
  figureColors[["col.evi"]] <- c(figureColors[["col"]][8], figureColors[["col"]][7]) # Descriptive Plus Evidence: dark orange and light orange 
  
  ## Add an alpha value to a colour
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  figureColors[["col.des.transparent"]] <- add.alpha(figureColors[["col.des"]], 0.8)
  figureColors[["col.evi.transparent"]] <- add.alpha(figureColors[["col.evi"]], 0.8)
  
#######################################################
## 0.2. Descriptives of Confidence Ratings per Study ##
#######################################################
  
  # Save and export image as .eps (8.5 x 11 inch)
  if(saveFigures) cairo_ps(file = '../R_Output/Images/ConfidenceDensitiesPerStudy.eps', onefile = TRUE, fallback_resolution = 600,
                           width = 8.5, height = 11)
  
  # Density plot for confidence ratings per study
  ggplot(dat, aes(x = 100*confidence.rating, y = study.confidence_order)) +
    geom_density_ridges(
      aes(fill = rep.outcome.label), 
      alpha = .8, from = -100, to = 100)  +
    geom_vline(xintercept = 0, linetype = 'dotted') + 
    theme_classic(base_size = 19) +   
    theme(legend.position='none') +
    coord_cartesian(xlim=c(-100.01, 100.01), clip = 'off') +
    ylab('') +
    xlab('Confidence Rating') + 
    scale_x_continuous(breaks=seq(-100,100, 50)) +
    scale_fill_cyclical(values = c('grey36','grey')) +
    geom_segment(aes(x=-100,xend=100,y=-Inf,yend=-Inf)) +
    geom_segment(aes(y=1, yend=27, x=-Inf, xend=-Inf)) +
    theme(axis.line=element_blank()) + 
    theme(plot.margin = unit(c(0.3,0.3,0.3,-0.3), 'cm'))
  graphics.off()
  # Image is saved as specified in cairo_ps
  
  ## Descriptives of Confidence Ratings per Study, Split per Condition ##  
  
  # Save and export image as .eps (8.5 x 11 inch)
  if(saveFigures) cairo_ps(file = '../R_Output/Images/ConfidenceDensitiesPerStudy_Description_Revision.eps', onefile = TRUE, fallback_resolution = 600,
           width = 8.5, height = 11)
  
  # Description Only density plot for confidence ratings per study
  ggplot(dat.description, aes(x = 100*confidence.rating, y = study.confidence_order)) +
    geom_density_ridges(
      aes(fill = rep.outcome.label), 
      alpha = .8, from = -100, to = 100)  +
    geom_vline(xintercept = 0, linetype = 'dotted') + 
    theme_classic(base_size = 19) +   
    theme(legend.position='none') +
    coord_cartesian(xlim=c(-100.01, 100.01), clip = 'off') +
    ylab('') +
    xlab('Confidence Rating') + 
    scale_x_continuous(breaks=seq(-100,100, 50)) +
    scale_fill_cyclical(values = figureColors[['col.des.transparent']]) +
    geom_segment(aes(x=-100,xend=100,y=-Inf,yend=-Inf)) +
    geom_segment(aes(y=1, yend=27, x=-Inf, xend=-Inf)) +
    theme(axis.line=element_blank()) + 
    theme(plot.margin = unit(c(0.3,0.3,0.3,-0.3), 'cm'))
  graphics.off()
  # Image is saved as specified in cairo_ps
  
  # Save and export image as .eps (8.5 x 11 inch)
  if(saveFigures) cairo_ps(file = '../R_Output/Images/ConfidenceDensitiesPerStudy_Evidence_Revision.eps', onefile = TRUE, fallback_resolution = 600,
           width = 8.5, height = 11)
  
  # Density plot for confidence ratings per study
  ggplot(dat.evidence, aes(x = 100*confidence.rating, y = study.confidence_order)) +
    geom_density_ridges(
      aes(fill = rep.outcome.label), 
      alpha = .8, from = -100, to = 100)  +
    geom_vline(xintercept = 0, linetype = 'dotted') + 
    theme_classic(base_size = 19) +   
    theme(legend.position='none') +
    coord_cartesian(xlim=c(-100.01, 100.01), clip = 'off') +
    ylab('') +
    xlab('Confidence Rating') + 
    scale_x_continuous(breaks=seq(-100,100, 50)) +
    scale_fill_cyclical(values = figureColors[['col.evi.transparent']]) +
    geom_segment(aes(x=-100,xend=100,y=-Inf,yend=-Inf)) +
    geom_segment(aes(y=1, yend=27, x=-Inf, xend=-Inf)) +
    theme(axis.line=element_blank()) + 
    #theme(axis.text.y = element_blank()) + 
    #theme(axis.ticks.y = element_blank()) +
    theme(plot.margin = unit(c(0.3,0.3,0.3,-0.3), 'cm'))
  graphics.off()
  
  par(mfrow=c(1,1)) # reset layout of plots to default
  
  ## Detailed Descriptives for Three Selected Studies ##  
  
  # Select three studies to illustrate
  Gneezy  <- dat[dat$study.label=='Gneezy et al. (2014)',]
  Gneezy$confidence.rating <- 100* Gneezy$confidence.rating
  Tversky <- dat[dat$study.label=='Tversky & Gati (1978)',]
  Tversky$confidence.rating <- 100* Tversky$confidence.rating
  Shah    <- dat[dat$study.label=='Shah et al. (2012)',]
  Shah$confidence.rating <- 100* Shah$confidence.rating
  
  # Save and export image as .eps (8.5 x 11 inch)
  if(saveFigures) cairo_ps(file = '../R_Output/Images/ThreeStudiesDescriptives.eps', onefile = TRUE, fallback_resolution = 600,
           width = 8.5, height = 11)
  
  # Three different studies plotted for details 
  par(mfrow = c(3, 1))
  
  # Gneezy et al. (2014)
  par(cex.main = 2.2, mar = c(5, 5.5, 5, 2) + 0.1, mgp = c(3, 1, 0), 
      cex.lab = 1.5, font.lab = 2, cex.axis = 1.8, bty = 'n', las = 1)
  hist(Gneezy$confidence.rating, freq = F, main = Gneezy$study.label[1], xlab = '', ylab = ' ', 
       xlim = c(-100, 100), axes = FALSE, breaks = 20, 
       ylim = c(0, 0.02), 
       yaxt = 'n', 
       col = 'grey')
  rug(jitter(Gneezy$confidence.rating), ticksize = -0.05)
  axis(1, at = seq(-100, 100, 50), lwd = 2, lwd.ticks = 2, line = 0.8)
  axis(2, labels = FALSE, lwd.ticks = 0)
  mtext('Confidence Rating', side = 1, line = 3.5, cex = 1.5, font = 1, adj = 0.5)
  mtext('Density', side = 2, line = 1.5, cex = 1.5, font = 1, las = 0)
  lines(density(Gneezy$confidence.rating, from = -100.01, to = 100.01), lwd = 3)
  lines(rep(mean(Gneezy$confidence.rating), 2), c(0, 2), lty = 2, col = 'black', 
        lwd = 1)
  
  # Tversky & Gati (1978) 
  par(cex.main = 2.2, mar = c(5, 5.5, 5, 2) + 0.1, mgp = c(3, 1, 0), 
      cex.lab = 2, font.lab = 2, cex.axis = 1.8, bty = 'n', las = 1)
  hist(Tversky$confidence.rating, freq = F, main = Tversky$study.label[1], xlab = '', ylab = ' ', 
       xlim = c(-100, 100), axes = FALSE, breaks = 20, ylim = c(0, 0.02), yaxt = 'n', 
       col = 'grey46')
  rug(jitter(Tversky$confidence.rating), ticksize = -0.05)
  axis(1, at = seq(-100, 100, 50), lwd = 2, lwd.ticks = 2, line = 0.8)
  axis(2, labels = FALSE, lwd.ticks = 0)
  mtext('Confidence Rating', side = 1, line = 3.5, cex = 1.5, font = 1, adj = 0.5)
  mtext('Density', side = 2, line = 1.5, cex = 1.5, font = 1, las = 0)
  lines(density(Tversky$confidence.rating, from = -100.01, to = 100.01), lwd = 3)
  lines(rep(mean(Tversky$confidence.rating), 2), c(0, 2), lty = 2, col = 'black', 
        lwd = 1)
  
  # Shah et al. (2012)
  par(cex.main = 2.2, mar = c(5, 5.5, 5, 2) + 0.1, mgp = c(3, 1, 0), 
      cex.lab = 1.5, font.lab = 2, cex.axis = 1.8, bty = 'n', las = 1)
  hist(Shah$confidence.rating, freq = F, main = Shah$study.label[1], xlab = '', ylab = ' ', 
       xlim = c(-100, 100), axes = FALSE, breaks = 20, ylim = c(0, 0.02), yaxt = 'n', 
       col = 'grey46')
  rug(jitter(Shah$confidence.rating), ticksize = -0.05)
  axis(1, at = seq(-100, 100, 50), lwd = 2, lwd.ticks = 2, line = 0.8)
  axis(2, labels = FALSE, lwd.ticks = 0)
  mtext('Confidence Rating', side = 1, line = 3.5 , cex = 1.5, font = 1, adj = 0.5)
  mtext('Density', side = 2, line = 1.5, cex = 1.5, font = 1, las = 0)
  lines(density(Shah$confidence.rating, from = -100.01, to = 100.01), lwd = 3)
  lines(rep(mean(Shah$confidence.rating), 2), c(0, 2), lty = 2, col = 'black', 
        lwd = 1)
  
  graphics.off()
  # Image is saved as specified in cairo_ps
  
  par(mfrow=c(1,1)) # reset layout of plots to default
  
  rm(Gneezy, Shah, Tversky)

##################################################################################################################
## 0.3. Preparatory Analysis: Relationship Between Bayes Factor of Original Studies and Replication Effect Size ##
##################################################################################################################
  
  original.bf.log             <- log(as.numeric(replication.results$OS.Bayes.factor))

  set.seed(4491)
  output.gibbs.sampler        <- spearmanGibbsSampler(original.bf.log, replication.effectsizes, progBar = TRUE)
  posterior.rho.bf            <- output.gibbs.sampler$rhoSamples

  # parameter estimation for Spearmans Rho
  median.and.credible.interval <- quantile(posterior.rho.bf, c(0.025, 0.5, 0.975))
  median.and.credible.interval
  # hypothesis testing: Savage-Dickey density ratio
  posterior.density           <- logspline(posterior.rho.bf)
  density.zero.point          <- dlogspline(   0, posterior.density)
  correction.factor.posterior <- 1 -plogspline(0, posterior.density)
  prior.density.zero.point    <- dbeta(0, 1, 1)
  # Bayes factor
  BFr0.rho.original.bf        <- prior.density.zero.point/(density.zero.point / correction.factor.posterior)
  BFr0.rho.original.bf
  
  ## Figure: Scatterplot of Log Bayes Factos and Replication Effect Size ##
  original.bf             <- as.numeric(replication.results$OS.Bayes.factor)
  
  # function to create tick marks on a log-scale
  log10Tck <- function(side, type){
    lim <- switch(side,
                  x = par('usr')[1:2],
                  y = par('usr')[3:4],
                  stop('side argument must be "x" or "y"'))
    at <- floor(lim[1]) : ceiling(lim[2])
    at <- 0:8
    return(switch(type,
                  minor = outer(1:9, 10^(min(at):max(at)))[,1:8],
                  major = 10^at,
                  stop('type argument must be "major" or "minor"')
    ))
  }
  
  if(saveFigures) cairo_ps(file = '../R_Output/Images/OriginalBFScatter.eps', onefile = TRUE, fallback_resolution = 600,
           width = 7.3, height = 4.16)
  
  # plot settings
  par(mar = c(5, 5, .1, 1) + 0.1, mgp = c(3.5, 1, 0))
  
  # Plot: original Bayes factor and replication effect size
  plot(original.bf, replication.effectsizes,
       xlim    = c(1, 10^8),
       ylim    = c(-0.2, 1),
       axes    = FALSE,
       log     = 'x',
       cex     = 2,
       cex.lab = 1.6,
       pch     = 21,
       bg      = c('grey36', 'grey')[replication.outcomes],
       xlab    = 'Bayes Factor Original Study',
       ylab    = 'Replication Effect Size (r)')
  # add axes and vertical line
  labels <- sapply((0:8),function(i)
    as.expression(bquote(10^ .(i))))
  #lines(c(10, 10), c(-.2,1), lty = 2, lwd = 1)
  axis(1, at=log10Tck('x','major'), tcl= -0.5, lwd = 2, lwd.ticks = 2, cex.axis = 1.4, las = 1, labels = labels)
  axis(1, at=log10Tck('x','minor'), tcl= -0.3, labels=NA, lwd = 2, lwd.ticks = 2, cex.axis = 1.4, las = 1)
  axis(2, lwd = 2, lwd.ticks = 2, cex.axis = 1.4, las = 1)
  
  graphics.off()
  # Image is saved as specified in cairo_ps
  
  par(mfrow=c(1,1)) # reset layout of plots to default

  rm(BFr0.rho.original.bf, correction.factor.posterior, density.zero.point,
     labels, median.and.credible.interval, prior.density.zero.point, posterior.rho.bf)

#########################################
## 1. Quality Check of Data And Design ##
#########################################

  # Research question: Do non-experts perform worse than guessing level? 
  # He: Omega is not smaller than .5
  # Hr: Omega is smaller than .5
  # The quality check fails, if we collect strong evidence (BF >= 10) in favor for Hr
  #
  # Prior distribution: Omega ~ Beta(1, 1); Kappa ~ Gamma(0.01, 0.01)
  alpha.omega <- 1
  beta.omega  <- 1

  # generate MCMC chain 
  set.seed(4491)
  mcmc.chain <-  genMCMC(data          = dat, 
                         sName         = 'sub.number',        # subject id
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
  posterior.nonexperts <- c(posterior.omega.check[, 1], posterior.omega.check[, 2])
  # save output
  # write.table(as.data.frame(posterior.nonexperts), file = '../R_Output/Samples/posterior_omega_nonexperts.txt', quote = FALSE, row.names = FALSE)

  # parameter estimation
  median.and.credible.interval <- quantile(posterior.nonexperts, c(0.025, 0.5, 0.975))
  median.and.credible.interval
  hist(posterior.nonexperts, las = 1)

  # hypothesis testing: encompassing prior approach (Klugkist et al., 2005)
  Iprior <- 0.5 # Iprior = 50% since prior distribution is centered around 0.5
  Ipost  <- sum(posterior.nonexperts < 0.5)/length(posterior.nonexperts)
  BFre   <- Ipost/Iprior; BFre
  
  rm(alpha.omega, beta.omega, Iprior, Ipost, BFre, median.and.credible.interval)

########################################
## 2.1 Test Differences in Conditions ##
########################################
  
  # This analysis will be conducted in JASP using a Bayesian independent samples t-test
  # with grouping variable 'condition' and dependent variable 'log.brier.score'
  #
  # corresponding JASP file: 'FinalData_BrierScores.jasp'
  # datafile: 'FinalData_BrierScores.csv'
  # 
  # Prior destribution: Default cauchy prior (width 0.707) on the effect size delta
  
  # Conclusion: The conditions differ, so we split the dataset in two parts
  nsubjects <- length(unique(dat.description$subject))
  subjects  <- unique(dat.description$subject)
  for (i in 1:nsubjects){
    dat.description$id[dat.description$subject==subjects[i]] <- i
  }
  
  nsubjects  <- length(unique(dat.evidence$subject))
  subjects   <- unique(dat.evidence$subject)
  for (i in 1:nsubjects){
    dat.evidence$id[dat.evidence$subject==subjects[i]] <- i
  }
  
  ## Figure: Boxplot per condition ## 
  if(saveFigures) pdf(file = '../R_Output/Images/Boxplots_Revision.pdf',
           width = 7.13, height = 6.12)
  
  # plot settings
  par(mar = c(5, 5, .1, 1) + 0.1, mgp = c(3.5, 1, 0), las = 1, cex.lab = 1.8, cex.axis = 1.5)
  # boxplot
  boxplot(data = out.brierscores, log.brier.score~condition, 
          col = c(figureColors[['col.des.transparent']][1],figureColors[['col.evi.transparent']][1]),
          outline = F, lwd = 2, lty = 1, cex = 2,
          ylim = c(-2.2,-0.6), axes = FALSE,
          boxwex = 0.5, xlab = 'Condition', ylab = 'Log (Brier Score)')
  axis(2, at = seq(-2.2, -0.6, by = 0.2), font =1)
  axis(1, at = c(1,2), labels = c('Description Only','Description Plus Evidence'))
  # Add data points
  mylevels <- levels(out.brierscores$condition)
  levelProportions <- summary(out.brierscores$condition)/nrow(out.brierscores)
  for(i in 1:length(mylevels)){
    thislevel <- mylevels[i]
    thisvalues <- out.brierscores[out.brierscores$condition==thislevel, "log.brier.score"]
    # take the x-axis indices and add a jitter, proportional to the N in each level
    myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/6)
    points(myjitter, thisvalues, pch=1, col=rgb(0,0,0,.9), cex = 1.5, lwd = 1.5) 
  }
  graphics.off()
  
  rm(myjitter, thislevel, thisvalues, mylevels, levelProportions, nsubjects, subjects)

#######################################
## 2.2 Estimate Accuracy Rate: Omega ##
#######################################
  
  # Research question: Is the accuracy rate higher than chance?
  # H0: Omega is exactly 0.5
  # Hr: Omega is higher than .5
  # If the data suggest strong evidence (BF > 10) in favor for the hypothesis that the Brier scores differ across 
  # conditions we will analyze the accuracy rate for each condition separately; otherwise, we use the full dataset
  #
  # Prior distribution: Omega ~ Beta(10, 10); Kappa ~ Gamma(0.01, 0.01)
  alpha.omega <- 10
  beta.omega  <- 10

## 2.2.1 Condition 'DescriptionOnly' ##

  # generate MCMC chain 
  set.seed(4491)
  mcmc.description.only <-  genMCMC(data  =  dat.description, 
                                    sName = 'id' ,           # subject id
                                    yName = 'guessed.correctly' , # responses '0' or '1
                                    aOmega        = alpha.omega,
                                    bOmega        = beta.omega,
                                    numSavedSteps = 20000 , 
                                    thinSteps     = 10,
                                    setSeed       = 4491)
  # convergence check
  posterior.omega.check <- cbind(mcmc.description.only[[1]][,'omega'], mcmc.description.only[[2]][,'omega'])
  nsamples              <- nrow(posterior.omega.check)
  plot( 1:nsamples, posterior.omega.check[, 1], type = 'l', las = 1, bty = 'n')
  lines(1:nsamples, posterior.omega.check[, 2], col = 'red')
  # extract posterior samples of omega 
  posterior.omega.description.only <- c(posterior.omega.check[, 1], posterior.omega.check[, 2])
  # save output
  # write.table(as.data.frame(posterior.omega.description.only), file = '../R_Output/Samples/posterior_omega_descriptiononly.txt', quote = FALSE, row.names = FALSE)
  
  # hypothesis testing: Savage-Dickey density ratio
  posterior.density           <- logspline(posterior.omega.description.only)
  density.zero.point          <- dlogspline(    0.5, posterior.density)
  correction.factor.posterior <- 1 - plogspline(0.5, posterior.density)
  prior.density.zero.point    <- dbeta(0.5, alpha.omega, beta.omega)
  correction.factor.prior     <- pbeta(0.5, alpha.omega, beta.omega, lower.tail = FALSE)
  # Bayes factor
  BFr0.description.only       <- (prior.density.zero.point / correction.factor.prior)/
                                 (density.zero.point / correction.factor.posterior) 
  BFr0.description.only
  # parameter estimation
  median.and.credible.interval <- quantile(posterior.omega.description.only, c(0.025, 0.5, 0.975))
  median.and.credible.interval
  # plot prior and posterior
  plot(density(posterior.omega.description.only), main = '', xlab = expression(omega), ylab = 'Density', bty = 'n', las = 1)
  # add prior distribution
  curve(dbeta(x, alpha.omega, beta.omega), 0, 1, lwd = 2, lty = 3, add = TRUE)
  # remove variables
  rm(posterior.omega.check, density.zero.point, correction.factor.posterior, correction.factor.prior,
     BFr0.description.only, median.and.credible.interval, mcmc.description.only)

## 2.2.2 Condition 'DescriptionPlusEvidence' ##

  # generate MCMC chain 
  set.seed(4491)
  mcmc.description.plus.evidence <-  genMCMC(data   =  dat.evidence, 
                                             sName  = 'id' ,                # subject id
                                             yName  = 'guessed.correctly' , # responses '0' or '1
                                             aOmega = alpha.omega,
                                             bOmega = beta.omega,
                                             numSavedSteps = 20000 , 
                                             thinSteps     = 10,
                                             setSeed       = 4491)
  # convergence check
  posterior.omega.check <- cbind(mcmc.description.plus.evidence[[1]][,'omega'], mcmc.description.plus.evidence[[2]][,'omega'])
  nsamples              <- nrow(posterior.omega.check)
  plot( 1:nsamples, posterior.omega.check[,1], type = 'l', las = 1, bty = 'n')
  lines(1:nsamples, posterior.omega.check[,2], col = 'red')
  # extract posterior samples of omega 
  posterior.omega.description.plus.evidence <- c(posterior.omega.check[, 1], posterior.omega.check[, 2])
  # save output
  # write.table(as.data.frame(posterior.omega.description.plus.evidence), file = '../R_Output/Samples/posterior_omega_descriptionEvidence.txt', quote = FALSE, row.names = FALSE)

  # hypothesis testing: Savage-Dickey density ratio
  posterior.density           <- logspline(posterior.omega.description.plus.evidence)
  density.zero.point          <- dlogspline(   0.5, posterior.density)
  correction.factor.posterior <- 1 - plogspline(0.5, posterior.density)
  prior.density.zero.point    <- dbeta(0.5, alpha.omega, beta.omega)
  correction.factor.prior     <- pbeta(0.5, alpha.omega, beta.omega, lower.tail = FALSE)
  # Bayes factor
  BFr0.description.plus.evidence <- (prior.density.zero.point / correction.factor.prior)/
                                    (density.zero.point / correction.factor.posterior) 
  BFr0.description.plus.evidence
  # parameter estimation
  median.and.credible.interval <- quantile(posterior.omega.description.plus.evidence, c(0.025, 0.5, 0.975))
  median.and.credible.interval
  # plot prior and posterior
  plot(density(posterior.omega.description.plus.evidence), main = '', xlab = expression(omega), ylab = 'Density', bty = 'n', las = 1)
  # add prior distribution
  curve(dbeta(x, alpha.omega, beta.omega), 0, 1, lwd = 2, lty = 3, add = TRUE)
  # remove variables
  rm(posterior.omega.check, density.zero.point, correction.factor.posterior, correction.factor.prior, 
     BFr0.description.plus.evidence, median.and.credible.interval,
     mcmc.description.plus.evidence)

## Figure: Posterior Densities of Accuracy Rate for Laypeople per Condition ##
  # dataframe with posterior samples of the groups
  df                   <- data.frame(posterior.omega.description.only, 
                                     posterior.omega.description.plus.evidence)
  df                   <- df %>% gather(condition, omega, 1:2)
  df$condition         <- as.factor(df$condition)
  levels(df$condition) <- c('Description Only',
                            'Description Plus Evidence')
  
  # Save and export image as .eps (8.18 x 4.88 inch)
  if(saveFigures) cairo_ps(file = '../R_Output/Images/AccuracyPerCondition_Revision.eps', onefile = TRUE, fallback_resolution = 600,
           width = 8.18, height = 4.88)
  
  par(cex.main = 2.2, mar = c(5, 5.5, 2, 2) + 0.1, mgp = c(3, 1, 0), 
      cex.lab = 2, font.lab = 1, cex.axis = 1.6, bty = 'n', las = 1)
  plot(x = c(.5, .75), y = c(0, 45), type = 'n', 
       xlab = expression(paste("Accuracy Rate ", omega, " (in %)")), ylab = '', 
       main = '', axes = FALSE)
  title(ylab = 'Density', line = 1)
  axis(2, lwd.tick = 0, labels = FALSE)
  axis(1, at = seq(.5, .75, by = .05), labels = seq(50, 75, by = 5))
  d <- lapply(split(df$omega, df$condition), density)
  Map(function(dens, col) polygon(dens, col = col), 
      dens = d, col = c(figureColors[['col.des.transparent']][1], figureColors[['col.evi.transparent']][1])) 
  
  graphics.off()
  # Image is saved as specified in cairo_ps
  
  par(mfrow=c(1,1)) # reset layout of plots to default
  
  # remove variables
  rm(d, df)
  
#########################################
## 3. Estimate Bayesian Spearman's Rho ##
#########################################
  
  # Research question: Can non-experts capture the evidential value of the study?
  # H0: Rho is 0 
  # Hr: Rho is bigger than 0 
  # If the data suggest strong evidence (BF > 10) in favor for the hypothesis that the Brier scores differ across 
  # conditions we will analyze the correlation between (averaged) confidence ratings and replication effect sizes  
  # for each condition separately; otherwise, we use the full dataset
  #
  # Prior distribution: Rho ~ Uniform(0, 1) for testing; Rho ~ Uniform(-1, 1) for estimation

## 3.1 Condition 'DescriptionOnly' ##

  confidence.description <- aggregate(confidence.rating ~ study, data = dat.description, FUN = mean)$confidence.rating # (range: -1,1)
  
  # generate MCMC chain (assuming wide uniform Beta distribution)
  set.seed(4491)
  output.gibbs.sampler <- spearmanGibbsSampler(confidence.description, replication.effectsizes, progBar = TRUE) 
  posterior.rho.description.only <- output.gibbs.sampler$rhoSamples
  # save output
  # write.table(as.data.frame(posterior.rho.description.only), file = '../R_Output/Samples/posterior_rho_descriptiononly.txt', quote = FALSE, row.names = FALSE)
  
  # hypothesis testing: Savage-Dickey density ratio
  posterior.density           <- logspline(posterior.rho.description.only)
  density.zero.point          <- dlogspline(   0, posterior.density)
  correction.factor.posterior <- 1 -plogspline(0, posterior.density)
  prior.density.zero.point    <- dbeta(0, 1, 1)
  # Bayes factor
  BFr0.rho.description.only   <- prior.density.zero.point/(density.zero.point / correction.factor.posterior)
  BFr0.rho.description.only
  # parameter estimation
  median.and.credible.interval <- quantile(posterior.rho.description.only, c(0.025, 0.5, 0.975))
  median.and.credible.interval
  # plot prior and posterior
  plot(density(posterior.rho.description.only), main = '', xlab = expression(rho), ylab = 'Density', bty = 'n', las = 1)
  # add prior distribution
  curve(dunif(x, -1, 1), -1, 1, lwd = 2, lty = 3, add = TRUE)
  
  rm(median.and.credible.interval, BFr0.rho.description.only, output.gibbs.sampler,
     posterior.density, density.zero.point, correction.factor.posterior, prior.density.zero.point,
     posterior.rho.description.only)

## 3.2 Condition 'DescriptionPlusEvidence'

  confidence.evidence <- aggregate(confidence.rating ~ study, data = dat.evidence, FUN = mean)$confidence.rating # (range: -1,1)
  
  # generate MCMC chain (assuming wide uniform Beta distribution)
  set.seed(4491)
  output.gibbs.sampler <- spearmanGibbsSampler(confidence.evidence, replication.effectsizes, progBar = TRUE)
  posterior.rho.description.evidence <- output.gibbs.sampler$rhoSamples
  # save output
  # write.table(as.data.frame(posterior.rho.description.evidence), file = '../R_Output/Samples/posterior_rho_descriptionEvidence.txt', quote = FALSE, row.names = FALSE)
  
  # hypothesis testing: Savage-Dickey density ratio
  posterior.density           <- logspline(posterior.rho.description.evidence)
  density.zero.point          <- dlogspline(   0, posterior.density)
  correction.factor.posterior <- 1 -plogspline(0, posterior.density)
  prior.density.zero.point    <- dbeta(0, 1, 1)
  # Bayes factor
  BFr0.description.plus.evidence  <- prior.density.zero.point/(density.zero.point / correction.factor.posterior)
  BFr0.description.plus.evidence
  # parameter estimation
  median.and.credible.interval <- quantile(posterior.rho.description.evidence, c(0.025, 0.5, 0.975))
  median.and.credible.interval
  # plot prior and posterior
  plot(density(posterior.rho.description.evidence), main = '', xlab = expression(rho), ylab = 'Density', bty = 'n', las = 1)
  # add prior distribution
  curve(dunif(x, -1, 1), -1, 1, lwd = 2, lty = 3, add = TRUE)
  
  rm(median.and.credible.interval, BFr0.description.plus.evidence, output.gibbs.sampler,
     posterior.density, density.zero.point, correction.factor.posterior, prior.density.zero.point,
     posterior.rho.description.evidence)

  ## Figure 7: Scatterplot of Confidence Ratings and Replication Effect Size ##
  
  ## Figure 7a: description only condition
  if(saveFigures) pdf(file = '../R_Output/Images/Scatterplot_DescriptionOnly_Revision.pdf',
           width = 6.17, height = 4.54)
  
  # plot settings
  par(mar = c(5, 5, .1, 0) + 0.1, mgp = c(3.5, 1, 0))
  
  plot(confidence.description, replication.effectsizes,
       xlim    = c(-1  , 1), 
       ylim    = c(-0.2, 1), 
       axes    = FALSE,
       cex     = 2,
       cex.lab = 1.8,
       pch     = 21, 
       bg      = figureColors[['col.des.transparent']][replication.outcomes],
       xlab    = 'Confidence Rating', 
       ylab    = 'Replication Effect Size (r)')
  # add axes and vertical line
  axis(1, at = seq(-1,1,.5)  , lwd = 2, lwd.ticks = 2, cex.axis = 1.6, las = 1)
  axis(2, at = seq(-.2,1, .2), lwd = 2, lwd.ticks = 2, cex.axis = 1.6, las = 1)
  lines(c(0, 0), c(-.2,1), lty = 2, lwd = 1)
  
  graphics.off()
  
  ## Figure 7b:  description plus evidence
  if(saveFigures) pdf(file = '../R_Output/Images/Scatterplot_DescriptionPlusEvidence_Revision.pdf', 
           width = 6.17, height = 4.54)
  
  # plot settings
  par(mar = c(5, 5, .1, 0) + 0.1, mgp = c(3.5, 1, 0))
  
  plot(confidence.evidence, replication.effectsizes,
       xlim    = c(-1  , 1), 
       ylim    = c(-0.2, 1), 
       axes    = FALSE,
       cex     = 2,
       cex.lab = 1.8,
       pch     = 21, 
       bg      = figureColors[['col.evi.transparent']][replication.outcomes],
       xlab    = 'Confidence Rating', 
       ylab    = 'Replication Effect Size (r)')
  # add axes and vertical line
  axis(1, at = seq(-1,1,.5)  , lwd = 2, lwd.ticks = 2, cex.axis = 1.6, las = 1)
  axis(2, at = seq(-.2,1, .2), lwd = 2, lwd.ticks = 2, cex.axis = 1.6, las = 1)
  lines(c(0, 0), c(-.2,1), lty = 2, lwd = 1)
  
  graphics.off()
  # Image is saved as specified in cairo_ps
  
  par(mfrow=c(1,1)) # reset layout of plots to default
  
  rm(confidence.description, confidence.evidence)
  