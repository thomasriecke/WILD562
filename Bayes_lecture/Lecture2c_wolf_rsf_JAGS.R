# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script will work use through our first 'path analysis,' a term
# often used to describe structural equation models without latent variables
#
# We'll use resource-selection data (used vs. available) from the Bow Valley
# wolf pack in Banff National Park and the R package piecewiseSEM. 

# Let's run two models... first we'll estimate wolf use as a function of 
# two covariates, deer use and elevation.
#
# Second, we'll estimate deer use as 'mediator' that is a function of
# elevation, and also functions as a covariate for wolf use
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(piecewiseSEM)
library(latex2exp)
library(semEff)
library(reshape2)
library(jagsUI)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data from GitHub
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat <- read.csv("C:/Users/thomas.riecke/Desktop/SEM_workshop/data/wolf_rsf.csv")
x <- url("https://github.com/thomasriecke/SEM_workshop/blob/main/data/wolf_rsf.csv?raw=true")
dat <- read.csv(x)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# subset to Bow Valley pack
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(dat)
dat <- subset(dat, pack == 'Bow valley')
dat <- subset(dat, !is.na(distacc))
table(dat$used)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a few quick maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
# map of elevation
ggplot(dat, aes(x = easting, y = northing, col = elevati)) + 
  geom_point(size = 2) + 
  coord_equal() +
  scale_colour_gradient(low = 'yellow', high = 'red')

# map of deer winter use
ggplot(dat, aes(easting, northing, col = deerwin)) + 
  geom_point(size = 2) + 
  coord_equal() +
  scale_colour_gradient(low = 'yellow', high = 'red')
# what's the relationships between deer winter use and elevation?
# how could we visualize or model that?


# map of used and random points
ggplot(dat, aes(easting, northing)) + 
  geom_point(size = dat$used*2+0.01, shape = dat$used) + 
  coord_equal()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# scale covariates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$z.elev <- as.numeric(scale(dat$elevati))
dat$z.deer <- as.numeric(scale(dat$deerwin))






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# JAGS model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("m_path.jags")
cat("
    model {

    # priors for everything assigned in a loop to minimize unnecessary coding mistakes
    for (j in 1:2){
      alpha[j] ~ dnorm(0, 0.1) # intercept and effect of elevation on deer
      beta[j+1] ~ dnorm(0, 0.1) # effects of elevation and deer on wolf use
    }
    beta[1] ~ dlogis(0,1) # intercept


    # st. dev. and precision for deer model
      sigma ~ dgamma(1,1)
      tau = 1/(sigma * sigma)

    for (i in 1:n){
      # here is the model for deer, lm(deer ~ elevation)
      d[i] ~ dnorm(alpha[1] + alpha[2] * e[i], tau)
      # here the model for wolves, glm(used ~ elevation + deer, family = 'binomial')
      logit(psi[i]) = beta[1] + beta[2] * e[i] + beta[3] * d[i]
      y[i] ~ dbern(psi[i])
    }

    }
    ",fill = TRUE)
sink()

jags.data <- list(n = nrow(dat), d = dat$z.deer, e = dat$z.elev, y = dat$used)
inits <- function(){list()}  
parameters <- c('alpha','beta','sigma')

# number of chains (nc), thinning rate (nt), number of iterations (ni), and number to burn-in
nc <- 4
nt <- 10
ni <- 20000
nb <- 10000



library(jagsUI)
Sys.time()
m <- jags(jags.data, inits, parameters, "m_path.jags", parallel = T, 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use the full posteriors to make plots with uncertainty
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iters <- length(m$sims.list$alpha[,1]) # we saved 4000 iterations
res <- 100

# elevation values for prediction
xE <- seq(min(dat$z.elev),max(dat$z.elev), length.out = res)

# expected deer (Ed), 
# direct effect of elevation (De), 
# total effect of elevation (Te)
Ed <- matrix(NA, iters, res)
De <- matrix(NA, iters, res) 
Te <- matrix(NA, iters, res)

# quantiles
qEd <- matrix(NA, res, 5)
qDe <- matrix(NA, res, 5)
qTe <- matrix(NA, res, 5)

for (j in 1:res){
  # expected value of deer given elevation
  Ed[,j] <- m$sims.list$alpha[,1] + m$sims.list$alpha[,2] * xE[j]
  # expected value of P(use) given direct effect of elevation
  De[,j] <- plogis(m$sims.list$beta[,1] + m$sims.list$beta[,2] * xE[j])  
  # expected value of P(use) given total effect of elevation
  Te[,j] <- plogis(m$sims.list$beta[,1] + m$sims.list$beta[,2] * xE[j] + m$sims.list$beta[,3] * Ed[,j])   
  
  qEd[j, ] <- quantile(Ed[,j], c(0.025,0.05,0.5,0.95,0.975))
  qDe[j, ] <- quantile(De[,j], c(0.025,0.05,0.5,0.95,0.975))
  qTe[j, ] <- quantile(Te[,j], c(0.025,0.05,0.5,0.95,0.975))
}

# melt these massive matrices into long format
mEd <- melt(Ed); names(mEd) <- c('iter','x','d')
mDe <- melt(De); names(mDe) <- c('iter','x','w')
mTe <- melt(Te); names(mTe) <- c('iter','x','w')




# plot of probability of use given total effect of elevation
# smooothscatter provides posterior density + median and 95% CrIs
smoothScatter(mTe$w ~ xE[mTe$x], las = 1, nrpoints = 0,
              ylab = 'Expected density of browse',
              xlab = 'Density of wolves')
lines(qTe[,3] ~ xE, lty = 1, lwd = 3, col = 'white')
lines(qTe[,1] ~ xE, lty = 2, lwd = 3, col = 'white')
lines(qTe[,5] ~ xE, lty = 2, lwd = 3, col = 'white')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot to compare direct vs. total effects
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.x <- seq(1400,3000, length.out = 9)
plot.s <- (plot.x - mean(dat$elevati))/sd(dat$elevati)
plot(qTe[,3] ~ xE, type = 'l', ylim = c(0,1),
     xlab = 'Elevation (m)', xaxt = 'n', lwd = 2,
     ylab = 'P(Used)', las = 1, cex.lab = 2)
lines(qTe[,1] ~ xE, lty = 2)
lines(qTe[,5] ~ xE, lty = 2)

lines(qDe[,3] ~ xE, col = 'grey50', lwd = 2)
lines(qDe[,1] ~ xE, col = 'grey50', lty = 3)
lines(qDe[,5] ~ xE, col = 'grey50', lty = 3)
axis(side = 1, at = plot.s, labels = plot.x)
points(jitter(dat$used, 0.1) ~ dat$z.elev, col = 'grey80')


