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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# glm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
glm1 <- glm(used ~ z.elev + z.deer, data = dat, family = 'binomial')
summary(glm1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot prediction of elevation effect from glm with elevation and deer
# as continuous covariates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res <- 100 # resolution of line
px <- seq(min(dat$z.elev), max(dat$z.elev), length.out = res)

py <- predict(glm1, 
              newdata = data.frame(z.elev = px, z.deer = rep(0, res)), 
              se.fit = T)
ey <- plogis(py$fit)
lci <- with(py, plogis(fit + qnorm(0.025)*se.fit))
uci <- with(py, plogis(fit + qnorm(0.975)*se.fit))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(ey ~ px, type = 'l', xlab = 'Elevation (m)', ylim = c(0,1), 
     ylab = 'P(used)', las = 1, cex.lab = 1.5, cex.axis = 1.25, lwd = 3)
lines(lci ~ px, lty = 2)
lines(uci ~ px, lty = 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# and a path analysis where we do the same thing, and also acknowledge 
# (technically, we assume) that elevation has an effect on deer use
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path1 <- psem(
  # this line of code is exactly the same as our glm above!
  glm(used ~ z.elev + z.deer, data = dat, family = 'binomial'),
  # here we indicate that we think elevation affects deer abundance
  glm(z.deer ~ z.elev, data = dat),
  data = dat
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save the summary as a list
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sump1 <- summary(path1, intercepts = T)
sump1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
piecewiseSEM:::plot.psem(path1, show = 'unstd')




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# using Monte Carlo simulation to make an effects plot 
# from the parameter estimates
# *this is a bit hacky, and I much prefer to work directly with posterior
# distributions, see the semEff package to learn more about how 
# to use bootstrapping to work with psem output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samps <- 3000
summary(px) # range of covariate values
p.wolf <- matrix(NA, samps, length(px))
p.deer <- matrix(NA, samps, length(px))
q.wolf <- matrix(NA, length(px), 3)
q.deer <- matrix(NA, length(px), 3)

# monte carlo sample parameter distributions
b0 <- rnorm(samps, sump1$coefficients$Estimate[1], sump1$coefficients$Std.Error[1])
b1 <- rnorm(samps, sump1$coefficients$Estimate[2], sump1$coefficients$Std.Error[2])
b2 <- rnorm(samps, sump1$coefficients$Estimate[3], sump1$coefficients$Std.Error[3])
a0 <- rnorm(samps, sump1$coefficients$Estimate[4], sump1$coefficients$Std.Error[4])
a1 <- rnorm(samps, sump1$coefficients$Estimate[5], sump1$coefficients$Std.Error[5])

str(path1)

# calculate predicted deer and predicted wolves given uncertainty
for (j in 1:res){
  p.deer[,j] <- a0 + a1 * px[j]
  p.wolf[,j] <- plogis(b0 + b1 * px[j] + b2 * p.deer[,j])
  q.deer[j,] <- quantile(p.deer[,j], c(0.025,0.5,0.975))
  q.wolf[j,] <- quantile(p.wolf[,j], c(0.025,0.5,0.975))
}





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot from piecewise SEM MC simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(q.wolf[,2] ~ px, type = 'l', xlab = 'Elevation (m)', ylim = c(0,1), 
     ylab = 'P(used)', las = 1, cex.lab = 1.5, cex.axis = 1.25, lwd = 3,
     col = 'grey50')
lines(q.wolf[,1] ~ px, lty = 2, col = 'grey50') # lower confidence interval
lines(q.wolf[,3] ~ px, lty = 2, col = 'grey50') # upper confidence interval



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create real covariate values for plotting and compare estimates from
# the two different models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# covariate values (for a subset of elevation)
real.elev <- seq(1400,2200, length.out = 5)
z.pos <- (real.elev - mean(dat$elevati))/sd(dat$elevati)


# plot results from a glm (holding deer constant at mean)
plot(ey ~ px, xaxt = 'n', type = 'l', xlab = 'Elevation (m)', ylim = c(0,1), 
     ylab = 'P(used)', las = 1, cex.lab = 1.5, cex.axis = 1.25, lwd = 3,
     xlim = c(-1,2.1))
lines(lci ~ px, lty = 2)
lines(uci ~ px, lty = 2)
axis(side = 1, at = z.pos, labels = real.elev)

# add results from psem that estimates change in deer as a function of elevation
lines(q.wolf[,2] ~ px, type = 'l', col = 'grey50', lwd = 3)
lines(q.wolf[,1] ~ px, type = 'l', col = 'grey50', lwd = 1, lty = 2)
lines(q.wolf[,3] ~ px, type = 'l', col = 'grey50', lwd = 1, lty = 2)

points(jitter(dat$used,0.1) ~ dat$z.elev, cex = 0.5, col = 'grey70')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# note that this is similar to a simpler model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m.elev <- glm(used ~ z.elev, data = dat, family = 'binomial')
summary(m.elev)
pred.elev <- plogis(predict(m.elev, newdata = data.frame(z.elev = px)))

lines(pred.elev ~ px, col = 'red')


m.elev$aic
glm1$aic

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a bit of foreshadowing, it seems that ungulate use has interspecific correlations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(2,2), mar = c(5.1,5.1,2.21,2.1), family = 'sans')
plot(jitter(elkwin) ~ jitter(deerwin), data = dat, las = 1,
     xlab = 'Deer winter use', ylab = 'Elk winter use')
plot(jitter(moosewin) ~ jitter(deerwin), data = dat, las = 1,
     xlab = 'Deer winter use', ylab = 'Moose winter use')
plot(jitter(moosewin) ~ jitter(elkwin), data = dat, las = 1,
     ylab = 'Moose winter use', xlab = 'Elk winter use')

cor.test(dat$elkwin,dat$deerwin)
cor.test(dat$moosewin,dat$deerwin)
cor.test(dat$moosewin,dat$elkwin)

ungulate <- data.frame(dat$elkwin,dat$deerwin,dat$moosewin)
cov(ungulate)
cor(ungulate)
