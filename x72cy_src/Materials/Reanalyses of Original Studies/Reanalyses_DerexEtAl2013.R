# Original Article
# 
# Derex, M., Beugin, M.-P., Godelle, B. & Raymond, M. Experimental evidence for the influence of group size on 
# cultural complexity. Nature 503, 389–391 (2013). DOI: 10.1038/nature12774
# 
# Research question:
#
# Are bigger groups characterized by more cultural diversity?
#
# How was this measured?
# 
# A computer game was set up to be a multiplayer game, consisting of either 2, 4, 8, or 16 players. 
# In this game the players had to collect points by either building arrowheads (a simple task), 
# or building fishing nets (a difficult task).
# 
# What was found?
# 
# The more players were playing the game, the more likely it was that they increased the cultural 
# complexity within the game (i.e., by building both items); chi-squared(16.3), p < .0001.

# Bayesian reanalysis
# 
# This reanalysis is based on the analysis conducted by Maarten Marsman on ordered
# binomial probabilities (2018; https://osf.io/z9mwq/)

rm(list = ls())
set.seed(4491)

#Prior parameters of the Beta prior distributions:
alpha <- 1
beta  <- 1

######################################################################
# Original Data
######################################################################

#Data:
x <- c(3, 4, 10, 11)
n <- c(15, 12, 12, 12)

#Log-Bayes factor H0 vs He

lbf.OE <- 3 * lbeta(alpha, beta)
lbf.OE <- lbf.OE + lbeta(sum(x) + alpha, sum(n) - sum(x) + beta)
lbf.OE <- lbf.OE - sum(lbeta(x + alpha, n - x + beta))

#Log-Bayes factor Hr (order restricted) versus He
I <- 0
for(it in 1:1e6) {
  theta <- rbeta(n = 4, alpha + x, beta + n - x)
  if(!is.unsorted(theta))
    I <- I + 1
}
lbf.RE <- log(I / it * factorial(4))

# Log-Bayes factor H0 versus Hr (order restricted)
lbf.OR <- lbf.OE - lbf.RE
BF.org <- c(lbf.OE, lbf.RE, lbf.OR)
# Bayes factor Hr versus H0 (order restricted)
1/exp(lbf.OR)

# The data suggests extreme evidence for the restricted hypothesis 
# that the probabilities for observing cultural complexity increases with group size;
# BFr0  = 21111. 



