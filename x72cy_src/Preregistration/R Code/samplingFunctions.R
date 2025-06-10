#  Estimate Bayesian Spearmans Rho
# 
#  R Code taken from:
#   van Doorn, J., & Ly, A. (2018, December 17). Bayesian Latent-Normal Inference for the Rank Sum Test, 
#   the Signed Rank Test, and Spearman’s rho. Retrieved from osf.io/gny35

spearmanGibbsSampler <- function(xVals, yVals, nSamples=1e4, progBar=FALSE, kappa=1) {
  
  if (progBar) {
    myBar <- txtProgressBar(min = 1, max = nSamples, initial = 1, char = "*", 
                            style = 3, width = 50)
  }
  
  n <- length(xVals)
  xRanks <- rank(xVals)
  yRanks <- rank(yVals)
  
  # Target: posterior samples of rho
  rhoSamples <- numeric(nSamples)
  mySigma <- diag(2)
  
  # Construct container for the auxilairy variables: Data augmentation/Gibbs
  latentX <- matrix(nrow=nSamples, ncol=n)
  latentY <- matrix(nrow=nSamples, ncol=n)
  
  # intitialise latent variables
  # intialise rho that is compatible with xVals, yVals
  #
  currentXVals <- qnorm((xRanks)/(2*n+1))
  currentYVals <- qnorm((yRanks)/(2*n+1))
  currentRho <- cor(currentXVals, currentYVals)
  
  for (j in 1:nSamples) {
    # Metropolis sampler: 
    # currentXVals and currentYVals first get realigned with the underlying current rho
    #
    for (i in sample(1:n)) {
      # Gibbs step go through pairs of z^{x}, z^{y} with current rho fixed
      #   Thus, current is respect to the Gibbs step. Furthermore, 
      #   here align latent variables to the rank
      #
      currentXRank <- xRanks[i]
      currentYRank <- yRanks[i]
      
      regressXOnY <- mean(currentYVals[yRanks==currentYRank])
      regressYOnX <- mean(currentXVals[xRanks==currentXRank])
      
      xBounds <- upperLowerTruncation(ranks=xRanks, values=currentXVals, currentRank=currentXRank, n=n)
      currentXVals[i] <- myTruncNormSim(xBounds[["under"]], xBounds[["upper"]], mu=(currentRho*regressXOnY), sd=sqrt(1-currentRho^2))
      
      yBounds <- upperLowerTruncation(ranks=yRanks, values=currentYVals, currentRank=currentYRank, n=n)
      currentYVals[i] <- myTruncNormSim(yBounds[["under"]], yBounds[["upper"]], mu=(currentRho*regressYOnX), sd=sqrt(1-currentRho^2))
    }    
    
    currentXVals <- (currentXVals-mean(currentXVals))/sd(currentXVals)
    currentYVals <- (currentYVals-mean(currentYVals))/sd(currentYVals)
    
    # Store Gibbs update
    latentX[j, ] <- currentXVals
    latentY[j, ] <- currentYVals
    
    # This is the sufficient statistic to evaluate the likelihood part of the MH
    rObs <- cor(currentXVals, currentYVals)
    
    # Do Metropolis step here
    # vectorise the runif(1) < acceptance
    chanceMechanism <- runif(nSamples)
    
    rhoNew <- metropolisOneStep(rhoCurrent=currentRho, rObs=rObs,n=n, alpha=1/kappa, 
                                chanceMechanism[j])
    
    # Store MH update
    rhoSamples[j] <- rhoNew # add proposal to samples if accepted
    currentRho <- rhoNew # add proposal to samples if accepted
    
    if (progBar){setTxtProgressBar(myBar, j)}
  }
  
  rhoSamples <- pearsonToSpearman(rhoSamples) # Transform Pearson's rho to Spearman's rho
  resultsList <- list(rhoSamples = rhoSamples, latentX = latentX, latentY = latentY)
  return(resultsList)
}

metropolisOneStep <- function (rhoCurrent, rObs, n, alpha=1, chanceMechanism) {
  # chanceMechanism is runif(1) vectorised
  # 
  zCurrent <- atanh(rhoCurrent)
  zCandidate <- rnorm(1, mean=atanh(rhoCurrent), sd=1/sqrt(n-3))
  rhoCandidate <- tanh(zCandidate)
  
  logAcceptance <- (alpha-n/2)*(log(1-rhoCandidate^2)-log(1-rhoCurrent^2))+
    n*((1-rhoCurrent*rObs)/(1-rhoCurrent^2)-(1-rhoCandidate*rObs)/(1-rhoCandidate^2))
  
  if (chanceMechanism <= exp(logAcceptance)) {
    return(rhoCandidate)
  } else {
    return(rhoCurrent)
  }
}

upperLowerTruncation <- function(ranks, values, currentRank, n, ranksAreIndices = FALSE) {
  
  if (currentRank == min(ranks)) {
    under <- -Inf
  } else {
    under <- max(values[ranks < currentRank])
  }
  
  if (currentRank == max(ranks)) {
    upper <- Inf
  } else {
    upper <- min(values[ranks > currentRank])
  }
  
  return(list(under=under, upper=upper))
}

myTruncNormSim <- function(lBound = -Inf, uBound = Inf, mu = 0, sd = 1){
  
  lBoundUni <- pnorm(lBound, mean = mu, sd = sd)
  uBoundUni <- pnorm(uBound, mean = mu, sd = sd)  
  mySample <- qnorm(runif(1, lBoundUni, uBoundUni), mean = mu, sd = sd)
  
  return(mySample)
}

pearsonToSpearman <- function(rho){
  
  mySpear <- (6/pi)*asin(rho/2)
  
  return(mySpear)
}

#  Estimate Hierarchical Theta 
# 
#  R Code taken from:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

# # Check that required packages are installed:
# want = c("parallel","rjags","runjags","compute.es")
# have = want %in% rownames(installed.packages())
# if ( any(!have) ) { install.packages( want[!have] ) }
# 
# # Load rjags. Assumes JAGS is already installed.
# try( library(rjags) )
# # Load runjags. Assumes JAGS is already installed.
# try( library(runjags) )
# try( runjags.options( inits.warning=FALSE , rng.warning=FALSE ) )
# 
# # set default number of chains and parallelness for MCMC:
# library(parallel) # for detectCores().
# nCores = detectCores() 
# if ( !is.finite(nCores) ) { nCores = 1 } 
# if ( nCores > 4 ) { 
#   nChainsDefault = 4  # because JAGS has only 4 rng's.
#   runjagsMethodDefault = "parallel"
# }
# if ( nCores == 4 ) { 
#   nChainsDefault = 3  # save 1 core for other processes.
#   runjagsMethodDefault = "parallel"
# }
# if ( nCores < 4 ) { 
#   nChainsDefault = 3 
#   runjagsMethodDefault = "rjags" # NOT parallel
# }


#------------------------------------------------------------------------------
# Function to generate and summarize MCMC chain

genMCMC = function( data , sName="s" , yName="y" ,  aOmega = 1, bOmega = 1,
                    numSavedSteps=50000 , saveName=NULL , thinSteps=1 ) { 
  require(rjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  # N.B.: This function expects the data to be a data frame, 
  # with one component named y being a vector of integer 0,1 values,
  # and one component named s being a factor of subject identifiers.
  y = data[,yName]
  s = as.numeric(data[,sName]) # ensures consecutive integer levels
  a = aOmega
  b = bOmega
  # Do some checking that data make sense:
  if ( any( y!=0 & y!=1 ) ) { stop("All y values must be 0 or 1.") }
  Ntotal = length(y)
  Nsubj = length(unique(s))
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    s = s ,
    a = a ,
    b = b ,
    Ntotal = Ntotal ,
    Nsubj = Nsubj
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
  for ( i in 1:Ntotal ) {
  y[i] ~ dbern( theta[s[i]] )
  }
  for ( sIdx in 1:Nsubj ) {
  theta[sIdx] ~ dbeta( omega*(kappa-2)+1 , (1-omega)*(kappa-2)+1 ) 
  }
  omega ~ dbeta( a , b ) 
  kappa <- kappaMinusTwo + 2
  kappaMinusTwo ~ dgamma( 0.01 , 0.01 )  # mean=1 , sd=10 (generic vague)
  }
  " # close quote for modelString
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  initsList = function() {
    thetaInit = rep(0,Nsubj)
    for ( sIdx in 1:Nsubj ) { # for each subject
      includeRows = ( s == sIdx ) # identify rows of this subject
      yThisSubj = y[includeRows]  # extract data of this subject
      resampledY = sample( yThisSubj , replace=TRUE ) # resample
      thetaInit[sIdx] = sum(resampledY)/length(resampledY) 
    }
    thetaInit = 0.001+0.998*thetaInit # keep away from 0,1
    meanThetaInit = mean( thetaInit )
    kappaInit = 100 # lazy, start high and let burn-in find better value
    return( list( theta=thetaInit , omega=meanThetaInit , 
                  kappaMinusTwo=kappaInit-2 ) )
  }
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "theta","omega","kappa") # The parameters to be monitored
  adaptSteps = 500             # Number of steps to adapt the samplers
  burnInSteps = 500            # Number of steps to burn-in the chains
  nChains = 2                  # nChains should be 2 or more for diagnostics 
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  # Create, initialize, and adapt the model:
  jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList , 
                          n.chains=nChains , n.adapt=adaptSteps )
  # Burn-in:
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel , n.iter=burnInSteps )
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                              n.iter=nIter , thin=thinSteps )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

# Plotting Functions

plotPosterior <- function(samples, posteriorDensity, xlabName, xLow, xHigh, yHigh){
  par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
      font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
  # set colors
  greycolhist <- rgb(0, 0, 0, alpha = 0.7)
  greycol     <- rgb(red = 190, green = 190, blue = 190, alpha = 280, maxColorValue = 280)  
  # set axis settings
  xlow        <- xLow
  xhigh       <- xHigh
  yhigh       <- yHigh
  nbreaks     <- 80
  # histogram
  y <- hist(samples, nbreaks, plot = F)
  
  plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), col  = greycolhist, type = "S",  lwd  = 2,  lty = 1,  
       xlim = c(xLow, xHigh), ylim = c(0, yHigh), xlab = " ", ylab = " ", main = " ", axes = F) 
  axis(1)
  axis(2)
  mtext(xlabName,  side=1, line = 2.8, cex=2)
  mtext("Density", side=2, line = 2.8, cex=2, las=0)
  # include in log spline density estimation
  par(new=T)
  plot(posteriorDensity, ylim = c(0, yHigh), xlim = c(xLow, xHigh), lty = 1, lwd = 2, axes = F, col = greycol)
}