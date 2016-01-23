#############################################
#University of Alaska Fairbanks
#FISH 421: Fisheries Population Dynamics
#Lab 8: Stock and Recruitment and Data Quality
#Instructor/lab author: Milo Adkison
#TA/code author: Phil Ganz
#############################################

# Objectives:
# 1. to fit a non-linear stock recruitment model 
#    to data
# 2. to examine the effects of error and contrast 
#    on estimates of optimal escapement

# Step 1. construct some stock-recruitment data

# The Beverton-Holt stock-recruitment relationship 
# can be written in the following form: 
#		R = a*S/(b+S)
#where:
# R = the number of recruits
# S = the number of spawners
# a = the maximum possible recruitment
# b = the number of spawners that produces 1/2 of 
#     the maximum recruitment

# Enter parameter values. (a=2000, b=1000)
a <- 2000
b <- 1000

# Enter a range of spawner abundance, from 200 to 400 in steps of 25. 
S <- seq(200,400, by=25)

# Use the Beverton-Holt formula above to calculate the number of recruits 
# for each of the 9 spawner abundances.
R <- a*S/(b+S)

# Graph the number of recruits versus the number of spawners. 
# Use a scatterplot with symbols only, no line.
plot(S,R)

# Step 2. fit a nonlinear model to the data

# A stochastic version of the Beverton-Holt model is:

#		R.obs = R*exp(sigma * z)
# where:
# R.obs = the observed recruitment
# R = a*S/(b+S), the deterministic prediction for recruitment 
# exp(sigma * z) = the log-normal random part of recruitment
# z = a normally-distributed random number with mean 0 and variance 1
# sigma = a scaling factor which controls how large random fluctuations are (=CV)	

# With normally-distributed error, you estimate model parameters by minimizing the 
# sum of squares (SSQ), i.e. the squared differences between observed and predicted 
# data. For the log-normal model, the best estimates of  ‘a’ and ‘b’ are those that 
# produce the lowest SSQ for log-transformed data. That is, you minimize:
# SSQ = sum of [ln(R.obs) - ln(R)]2

# SSQ is a function of parameters (a.est and b.est) that we will change until SSQ is 
# optimized (in this case minimized). The optimizing function that we will use in R 
# uses a vector of parameters, so when we write our SSQ function, we will make a.est 
# and b.est elements of a larger parameter vector, par. put in a couple of starting 
# estimates such as 3000 and 500.
par <-c(a.est=3000, b.est=500)

## Right now, our "observed" recruitment is just the deterministic relationship we 
## calculated with the "true" parameter values. This will change later.
R.obs <- R

## Function for calculating predicted recruitment and SSQ
SSQ <- function(par){
R.pred <- par[1]*S/(par[2]+S)
ssq <- sum((log(R.obs)-log(R.pred))^2)
return(ssq)
}

## Try it out!
SSQ(par)

## Now we can optimize this function (our objective function): 
fit <- optim(par,SSQ)

# View parameter estimates
fit$par 

# There are no errors in our observations at this point (R.obs=R 
# because sigma=0), so our parameter estimates are very close to the 
# actual parameter values. It is important to test a model on simulated 
# data for which the true parameter values are known. In this case, our 
# model passes the test. Our next step is to see how our model performs 
# when there is random variation in recruitment.

# Step 3. Examine the effect of random error in recruitment on estimates 
# of ‘a’ and ‘b’.

## So that you get the same random numbers as me:
set.seed(269)

# We stated previosly that z has a standard normal distribution. We need 
# to generate a z value for each observation. In R, the default mean and 
# standard deviation used by the rnorm function are 0 and 1, respectively. 
# We just need to tell it how many random numbers we want:
z <- rnorm(9)

# Now, set sigma(CV) equal to 0.2 and incorporate random error into our new 
# observed recruitment values:
CV <- 0.2
R.obs <- R*exp(CV*z)

# Re-estimate parameter values using optim. What parameter estimates does 
# it find? How wrong are they?
fit <- optim(par,SSQ)
fit$par 

# Generate a new set of random numbers to replace the last ones.
z <- rnorm(9)
R.obs <- R*exp(CV*z)

# Optimize again. 
fit <- optim(par,SSQ)
fit$par 

# Are the estimates the same? Why are the estimates wrong, and why do they 
# differ from the previous ones?

# Increase the amount of randomness in your simulation and see how it affects 
# the estimates. Change the CV from 0.2 to 0.4 and optimize again. 
CV <- 0.4
R.obs <- R*exp(CV*z)

fit <- optim(par,SSQ)
fit$par 

# Are the estimates better? Why or why not?

# Hilborn and Walters say that poor estimates often occur not because the estimation 
# procedure is bad, but rather that the data are uninformative. Change the stock sizes 
# used from 200-400 to 100 to 1700 in steps of 200. Optimize again. 
S <- seq(100,1700, by=200)
R <- a*S/(b+S)
R.obs <- R*exp(CV*z)

fit <- optim(par,SSQ)
fit$par 

# Set CV back to 0.2 and try again
CV <- 0.2
R.obs <- R*exp(CV*z)

fit <- optim(par,SSQ)
fit$par 

# Are the estimates any better? Why or why not?