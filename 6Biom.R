#############################################
#University of Alaska Fairbanks
#FISH 421: Fisheries Population Dynamics
#Lab 6: Biomass Dynamics
#Instructor/lab author: Milo Adkison
#TA/code author: Phil Ganz
#############################################

# The file you’re going to receive contains catch and 
# effort data for two fish stocks. The first, the Aleutian 
# Island radiofish (Sebastes marconi), was first commercially 
# exploited in 1930. This fish contains an oil invaluable in 
# manufacturing vacuum tubes, so it experienced a rapidly 
# developing fishery. The fishery ceased during World War II, 
# and after the war the development of transistors reduced 
# demand for the oil (vacuum tubes are currently used only 
# in primitive equipment, like the space shuttle), so the 
# fishery resumed at a very low level of effort.

# The second stock is Pacific halibut in Area 2, which is the 
# southern part of their range somewhere. This data will be used 
# for your homework assignment.

## Load data. Make sure that files are in your working directory.
## You can view your current working directory's path with getwd() 
## and set it with setwd().
path <- paste(getwd(),"/BIOMDAT.xls", sep="")

##install.packages("readxl")
library(readxl)

## Create object of radiofish data
rad <- read_excel(path, sheet=1)

## Create object of halibut data
hal <- read_excel(path, sheet=2)

# Step 1. Fit the dynamic version of the Schaefer model to the radiofish data.
# a. Assume that catch data are observed with normal error (not lognormal)
# b. Assume the initial biomass was equal to K
# c. Assume the fishery occurs at the beginning of the year, prior to production of new biomass
# d. Make separate columns for biomass at the beginning of the year, predicted catch, biomass 
#    after the fishery, new production, and the residual of the predicted catch
# e. Estimate r, K, q (and MSY, Emsy) by minimizing SUMSQ(predicted catch – observed)
# f. graph the series of observed and predicted catches versus year

## Create a vector of our parameters. It is important that we have them as a vector
## in addition to individual objects because the optim() function only optimizes over
## one object. First we will optimize over the entire parameter set. Later, we will
## optimize over individual parameters.

r <- 0.4
K <- 500
q <- 0.001

par <- c(r=r, K=K, q=q)

## First, create a function that calculates the various values of interest:
rad.pred <- function(par){
	# d. Make separate columns for biomass at the beginning of the year, predicted catch, biomass 
    #    after the fishery, new production, and the residual of the predicted catch
	rad$Biomass_start[1] <- par["K"] #	b.	Assume the initial biomass was equal to K
	for (i in 1:length(rad$Year)){
		for (j in 2:length(rad$Year)){
		rad$Pred[i] <- par["q"] * rad$Effort[i] * rad$Biomass_start[i]
		rad$Biom_end[i] <- rad$Biomass_start[i] - rad$Pred[i] #	c.	Assume the fishery occurs at the beginning of the year, prior to production of new biomass
		rad$Prod[i] <- par["r"] * rad$Biom_end[i] * (1 - rad$Biom_end[i]/par["K"])
		rad$Biomass_start[j] <- rad$Biom_end[j-1] + rad$Prod[j-1]
	}}
	rad$Resid <- rad$Catch - rad$Pred #	a.	Assume that catch data are observed with normal error (not lognormal)
	return(rad)
}

## Next, create a function that calculates the sum of squared residuals calculated
## by the rad.pred() function:
SSQ <- function(par){
	SSQ <- sum(rad.pred(par)$Resid^2)
	return(SSQ)
}

# e. Estimate r, K, q (and MSY, Emsy) by minimizing the sum of squared catch residuals (sum(rad.pred(par)$Resid^2)).
fit <- optim(par, SSQ)

## Use the estimated parameter values to add predicted values to data 
rad <- rad.pred(fit$par)

# f. graph the series of observed and predicted catches versus year
plot(rad$Year, rad$Catch, type="b", xlab="Year", ylab="Catch")
lines(rad$Year, rad$Pred, type="b", col=2)
legend("topright", legend=c("Observed", "Predicted"), col=c(1,2), lwd=c(1,1))

# Step 2. Use subsets of the radiofish data.
# a. Re-estimate r, K, q, MSY, and Emsy using only the data from 1930-1940
SSQ.early <- function(par){
	SSQ <- sum(subset(rad.pred(par)$Resid, rad.pred(par)$Year<41)^2)
	return(SSQ)
}

fit.early <- optim(par, SSQ.early)

# b. re-estimate using data from 1941-1955.
SSQ.late <- function(par){
	SSQ <- sum(subset(rad.pred(par)$Resid, rad.pred(par)$Year>40)^2)
	return(SSQ)
}

fit.late <- optim(par, SSQ.late)

plot(rad$Year, rad$Catch, type="b", xlab="Year", ylab="Catch", ylim=c(0,170))
lines(rad$Year, rad$Pred, type="b", col=2)
lines(rad.pred(fit.early$par)$Year, rad.pred(fit.early$par)$Pred, type="b", col=3)
lines(rad.pred(fit.late$par)$Year, rad.pred(fit.late$par)$Pred, type="b", col=4)
legend("topright", legend=c("Observed", "Predicted", "1930-1940", "1941-1955"), col=c(1,2,3,4), lwd=c(1,1,1,1))

# c. Explain the differences among your three sets of estimates

# Step 3. Look at estimability of different parameters using data from 1941-1955

## First, let's modify rad.pred() and SSQ.late() so that they are functions of r, K, and q individually.
## This will allow us to fix parameters in optim().
rad.pred <- function(r,K,q){
	rad$Biomass_start[1] <- K 
	for (i in 1:length(rad$Year)){
		for (j in 2:length(rad$Year)){
		rad$Pred[i] <- q * rad$Effort[i] * rad$Biomass_start[i]
		rad$Biom_end[i] <- rad$Biomass_start[i] - rad$Pred[i]
		rad$Prod[i] <- r * rad$Biom_end[i] * (1 - rad$Biom_end[i]/K)
		rad$Biomass_start[j] <- rad$Biom_end[j-1] + rad$Prod[j-1]
	}}
	rad$Resid <- rad$Catch - rad$Pred
	return(rad)
}

SSQ.late <- function(r,K,q){
	SSQ <- sum(subset(rad.pred(r,K,q)$Resid, rad.pred(r,K,q)$Year>40)^2)
	return(SSQ)
}

# a. fix K at 500, q at .001, and estimate r using data from 1941-1955.
# How well do the predicted catches match the observed catches over that period?

## Switch optimization method to "BFGS" because the default "Nelder-Mead"
## method is unreliable when estimating only one parameter.
fit.r <- optim(r, SSQ.late, method="BFGS", K=K, q=q)

# b. fix r at 0.4, q at .001, and estimate K using data from 1941-1955. 
# How well do the predicted catches match the observed catches over that period?
fit.K <- optim(K, SSQ.late, method="BFGS", r=r, q=q)

# c. fix r at 0.4, K at 500, and estimate q using data from 1941-1955. 
# How well do the predicted catches match the observed catches over that period?
fit.q <- optim(q, SSQ.late, method="BFGS", r=r, K=K)

## Plot the predicted catches from the different estimates that used 
plot(rad$Year, rad$Catch, type="b", xlab="Year", ylab="Catch", ylim=c(0,150))
lines(rad.pred(fit$par[1], fit$par[2], fit$par[3])$Year, rad.pred(fit$par[1], fit$par[2], fit$par[3])$Pred, type="b", col=2)
lines(rad.pred(fit.r$par, K, q)$Year, rad.pred(fit.r$par, K, q)$Pred, type="b", col=3)
lines(rad.pred(r, fit.K$par, q)$Year, rad.pred(r, fit.K$par, q)$Pred, type="b", col=4)
lines(rad.pred(r, K, fit.q$par)$Year, rad.pred(r, K, fit.q$par)$Pred, type="b", col=5)
legend("topright", legend=c("Observed", "All params", "r", "K", "q"), col=c(1,2,3,4,5), lwd=c(1,1,1,1,1))

# d. Explain why some of the predicted catches match the observed catches well even 
#    though the parameter estimates differ from the better ones estimated from the full data series.

# e. Which parameters can be reliably estimated from the 1941-1955 time series?

# Biomass Dynamics homework

# 1. Fit the model we used for the radiofish to the Pacific halibut data for the time period 1910-1930.

# 2. Using these parameter estimates, predict the catch over the period 1931-1957. 
#    Explain why the predicted and observed don’t match.

# 3. What features of the 1910-1930 data series lead to poor estimates of parameter values? 
#    Which values are most likely to be poorly estimated?

# 4. Now fit the model to the entire set of halibut data (1910-1957).

# 5. Plot the residuals versus time. What other phenomenon might be causing difficulties in estimating a constant set of parameter values?
