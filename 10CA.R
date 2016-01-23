#############################################
#University of Alaska Fairbanks
#FISH 421: Fisheries Population Dynamics
#Lab 10: Catch-at-age Analysis
#Instructor/lab author: Milo Adkison
#TA/code author: Phil Ganz
#############################################

# Your data file contains several years of catch data by age for a 
# stock of trout.

## Load data. Make sure that files are in your working directory.
## You can view your current working directory's path with getwd() 
## and set it with setwd().
path <- paste(getwd(),"/CAData.xls", sep="")

##install.packages("readxl")
library(readxl)

## Create object of the data. Skip the first 6 empty rows when importing:
trout <- read_excel(path, sheet=1, skip = 6)

# Part 1a: Fit a statistical catch-at-age model with the following characteristics to the catch data only:

# Use the Baranov catch equation assuming continuous fishing.
# Assume a constant natural mortality rate of 0.4
# Include the initial age structure and each year’s recruitment as parameters
# Assume F_a,t is separable into s_a*f_t, where each f_t is a parameter. Thus:

# N_t+1,a+1 = N_t,a * exp(-f_t*s_a - M)
# C_t,a = N_t,a * f_t*s_a/(f_t*s_a + M) * (1 - exp(-f_t*s_a - M))
 
# Assume vulnerability to the fishery varies by age according to the following equation:  

# s_a = 1/(1 + exp(B * (a50-a)))

# Assume the catch data have a normal error structure.

# Hints:

# Many of the parameter values will be in very different units (thousands versus thousandths). 
# Have optim() adjust a dummy value and set the parameter equal to this dummy value times some constant, 
# so that optim() is adjusting parameter values that are all in the neighborhood of 1.

# Make sure you start with reasonable estimates for the parameter values. 

# You can tell if you’re at a good fit if: 
# (1) The observed total catch and the predicted total catch match pretty closely, and 
# (2) the SSQ(catch) should be close to 230,000.

# Years
t <- c(1968:1979)

# Ages
a  <- c(3:7)

# Parameters
# Initial abundance in year 1968:
N_68   <- c(N_4_68=2, 
	     	N_5_68=2, 
		 	N_6_68=2, 
		 	N_7_68=2)
# Recruitment
R   <- c(R_68=10, 
		 R_69=10, 
		 R_70=10, 
		 R_71=10, 
		 R_72=10, 
		 R_73=10, 
		 R_74=10, 
		 R_75=10, 
		 R_76=10, 
		 R_77=10, 
		 R_78=10, 
		 R_79=10)
# Fishing intensity
f   <- c(f_68=-1,
		 f_69=-1,
		 f_70=-1,
		 f_71=-1,
		 f_72=-1,
		 f_73=-1,
		 f_74=-1,
		 f_75=-1,
		 f_76=-1,
		 f_77=-1,
		 f_78=-1,
		 f_79=-1)
# Selectivity
a50 <- c(a50=5)
B   <- c(B=1)

# Natural mortality
M   <- 0.4

# Compile all parameters to be estimated
par <- c(N_68, R, f, a50, B)

# Now let's write the functions that will ultimately describe the predicted catch:
selectivity <- function(a50,B){
   a50 <- abs(a50)
   B <- abs(B)
   s <- NA
   for (i in 1:length(a)){
   s[i] <- 1/(1 + exp(B * (a50-a[i])))
   }
   return(s)}

fish_mort <- function(s,f){
	f <- abs(f)
	F <- matrix(NA, nrow=length(t), ncol=length(a))
	for (i in 1:length(t)){
		for (j in 1:length(a)){
    F[i,j] <- s[j] * f[i]
		}
	}
	return(F)}

abundance <- function(F,N_68,R){
	N_68 <- abs(N_68) * 1000
	R <- abs(R) *1000
	N <- matrix(NA, nrow=length(t), ncol=length(a))
	for (j in 2:length(a)){
    N[1,j] <- N_68[j-1]}
    for (i in 1:length(t)){
    N[i,1] <- R[i]}
    for (i in 1:(length(t)-1)){
    for (j in 1:(length(a)-1)){	
    N[i+1,j+1] <- N[i,j] * exp(-F[i,j]-M)	
    	}
    }
    return(N)}

catch <- function(par){
	s <- selectivity(par["a50"], par["B"])
	F <- fish_mort(s,par[c("f_68", "f_69", "f_70", "f_71", "f_72", "f_73", "f_74", "f_75", "f_76", "f_77", "f_78", "f_79")])
	N <- abundance(F,par[c("N_4_68", "N_5_68", "N_6_68", "N_7_68")],par[c("R_68", "R_69", "R_70", "R_71", "R_72", "R_73", "R_74", "R_75", "R_76", "R_77", "R_78", "R_79")])
	C <- matrix(NA, nrow=length(t), ncol=length(a))
	for (i in 1:length(t)){
    	for (j in 1:length(a)){
    C[i,j] <- N[i,j] * F[i,j]/(F[i,j]+M) * (1 - exp(-F[i,j]-M))}}
    return(C)
}

# Objective function
SSQ.catch <- function(par){
	C <- catch(par)
	SSQ <- sum((trout[1:12,2:6]-C)^2)
	return(SSQ)
}

# Optimize. We boost the number of iterations to 500 due to the complexity of the 
# function being optimized:
fit <- optim(par, SSQ.catch, method="L-BFGS-B", control=list(maxit=500))

# Check to see if SSQ is close to 230,000:
fit$value

# Graph the total catch by year versus the predicted total by year. How well do they agree?
obs.catch <- trout[1:12,2:6]
pred.catch <- catch(fit$par)

plot(t,rowSums(obs.catch), type="b", col=4, lwd=2, xlab="Year", ylab="Catch")
lines(t,rowSums(pred.catch), type="b", col=2, lwd=2)
legend("topleft", c("Observed","Predicted"), col=c(4,2), lwd=c(2,2))

# Why are the catches of the oldest age in the first year and the youngest age in the last year 
# fit perfectly by your model?

# Change M to 0.1 and then to 0.7. Run optim() with each to minimize the SSQ by adjusting all parameters. 
M <- 0.1
fit.0.1 <- optim(par, SSQ.catch, method="L-BFGS-B", control=list(maxit=500))

M <- 0.7
fit.0.7 <- optim(par, SSQ.catch, method="L-BFGS-B", control=list(maxit=500))

# How much does SSQ change? 
fit$value
fit.0.1$value
fit.0.7$value

# What does this say about how well the catch data determines a particular value of M? 
# How much does changing M affect your estimated abundances and f’s? 
pred.s <- selectivity(fit$par["a50"], fit$par["B"])
s.0.1 <- selectivity(fit.0.1$par["a50"], fit.0.1$par["B"])
s.0.7 <- selectivity(fit.0.7$par["a50"], fit.0.7$par["B"])	

pred.F <- fish_mort(pred.s,fit$par[c("f_68", "f_69", "f_70", "f_71", "f_72", "f_73", "f_74", "f_75", "f_76", "f_77", "f_78", "f_79")])
F.0.1 <- fish_mort(s.0.1,fit.0.1$par[c("f_68", "f_69", "f_70", "f_71", "f_72", "f_73", "f_74", "f_75", "f_76", "f_77", "f_78", "f_79")])
F.0.7 <- fish_mort(s.0.7,fit.0.7$par[c("f_68", "f_69", "f_70", "f_71", "f_72", "f_73", "f_74", "f_75", "f_76", "f_77", "f_78", "f_79")])

pred.N <- abundance(pred.F,fit$par[c("N_4_68", "N_5_68", "N_6_68", "N_7_68")],fit$par[c("R_68", "R_69", "R_70", "R_71", "R_72", "R_73", "R_74", "R_75", "R_76", "R_77", "R_78", "R_79")])
N.0.1 <- abundance(F.0.1,fit.0.1$par[c("N_4_68", "N_5_68", "N_6_68", "N_7_68")],fit.0.1$par[c("R_68", "R_69", "R_70", "R_71", "R_72", "R_73", "R_74", "R_75", "R_76", "R_77", "R_78", "R_79")])
N.0.7 <- abundance(F.0.7,fit.0.7$par[c("N_4_68", "N_5_68", "N_6_68", "N_7_68")],fit.0.7$par[c("R_68", "R_69", "R_70", "R_71", "R_72", "R_73", "R_74", "R_75", "R_76", "R_77", "R_78", "R_79")])

# Plot abundance estimates for each M assumption
plot(t,rowSums(pred.N), type="b", col=4, lwd=2, xlab="Year", ylab="Estimated abundance", ylim=c(0,40000))
lines(t,rowSums(N.0.1), type="b", col=2, lwd=2)
lines(t,rowSums(N.0.7), type="b", col=3, lwd=2)
legend("topright", c("M=0.4","M=0.1","M=0.7"), col=c(4,2,3), lwd=c(2,2,2))

# Plot fishing intensity estimates for each M assumption
plot(t,fit$par[c("f_68", "f_69", "f_70", "f_71", "f_72", "f_73", "f_74", "f_75", "f_76", "f_77", "f_78", "f_79")], type="b", col=4, lwd=2, xlab="Year", ylab="Fishing intensity (f)", ylim=c(-4.5,1))
lines(t,fit.0.1$par[c("f_68", "f_69", "f_70", "f_71", "f_72", "f_73", "f_74", "f_75", "f_76", "f_77", "f_78", "f_79")], type="b", col=2, lwd=2)
lines(t,fit.0.7$par[c("f_68", "f_69", "f_70", "f_71", "f_72", "f_73", "f_74", "f_75", "f_76", "f_77", "f_78", "f_79")], type="b", col=3, lwd=2)
legend("topright", c("M=0.4","M=0.1","M=0.7"), col=c(4,2,3), lwd=c(2,2,2))

# How important is the value of M?

# Part 2. Your spreadsheet also contains the catch rates from a survey that is assumed 
# to be proportional to the abundance of the stock. Use the survey data as auxiliary 
# information to improve your estimates. Assume the survey index has a log-normal error 
# distribution.

# a. Leaving all other parameters values alone, adjust the catchability parameter for your 
# survey to its best estimate. 

# Set a starting value for q
q <- c(q=0.3)

# Write a function for survey SSQ
SSQ.survey <- function(q, pred.N){
	pred.srv <- q * pred.N
	SSQ <- sum((log(trout$survey)-log(pred.srv))^2)
	return(SSQ)
}

# Optimize
fit.q <- optim(q, SSQ.survey, method="Brent", lower=0, upper=1, pred.N=rowSums(pred.N))

# Does your estimate of the abundance of the population agree with the survey index? 
plot(t, rowSums(pred.N), type="b", col=4, lwd=2, xlab="Year", ylab="Abundance", ylim=c(0,16000))
lines(t, trout$survey/fit.q$par, type="b", col=2, lwd=2)
legend("bottomleft", c("Predicted", "Survey/q"), col=c(4,2), lwd=c(2,2))

# Where are the major discrepancies? 

# b. Use optim() to adjust all parameters to fit only the survey data as well as possible, 
# ignoring the catch data. 
par <- c(par, q)

survey <- function(par){
	s <- selectivity(par["a50"], par["B"])
	F <- fish_mort(s,par[c("f_68", "f_69", "f_70", "f_71", "f_72", "f_73", "f_74", "f_75", "f_76", "f_77", "f_78", "f_79")])
	N <- abundance(F,par[c("N_4_68", "N_5_68", "N_6_68", "N_7_68")],par[c("R_68", "R_69", "R_70", "R_71", "R_72", "R_73", "R_74", "R_75", "R_76", "R_77", "R_78", "R_79")])
	I <- NA
	for (i in 1:length(t)){
    I[i] <- abs(par["q"]) * sum(N[i,])}
    return(I)
}

SSQ.survey <- function(par){
	pred.srv <- survey(par)
	SSQ <- sum((log(trout$survey)-log(pred.srv))^2)
	return(SSQ)
}

fit.srv <- optim(par, SSQ.survey, method="L-BFGS-B", control=list(maxit=500))

# Graph the predicted and observed survey index versus year. 
obs.srv  <- trout$survey
pred.srv <- survey(fit.srv$par)

plot(t,obs.srv, type="b", col=4, lwd=2, xlab="Year", ylab="Survey")
lines(t,pred.srv, type="b", col=2, lwd=2, lty=2)
legend("topright", c("Observed","Predicted"), col=c(4,2), lwd=c(2,2))

# How well is it possible to fit the data? 

# How well does the model fit to the survey data only match the catch data? 
pred.catch <- catch(fit.srv$par)

plot(t,rowSums(obs.catch), type="b", col=4, lwd=2, xlab="Year", ylab="Catch", ylim=c(0,10000))
lines(t,rowSums(pred.catch), type="b", col=2, lwd=2)
legend("bottomright", c("Observed","Predicted"), col=c(4,2), lwd=c(2,2))

# Explain your results.

# c. Fit the catch data and the survey data simultaneously, by minimizing the SSQ(catch) + 1,000,000*SSQ(survey). 
SSQ <- function(par){
	SSQ <- SSQ.catch(par) + 1000000 * SSQ.survey(par)
	return(SSQ)
}

# Since our objective fuction is even more complex now that it has two sources of data,
# boost the maximum number of iterations again. Feel free to grab a beer, this will take a while.
fit.both <- optim(par, SSQ, method="L-BFGS-B", control=list(maxit=5000))

# How well does this compromise set of parameter values fit the catch data and the survey data? 
pred.catch <- catch(fit.both$par)

plot(t,rowSums(obs.catch), type="b", col=4, lwd=2, xlab="Year", ylab="Catch")
lines(t,rowSums(pred.catch), type="b", col=2, lwd=2)
legend("bottomright", c("Observed","Predicted"), col=c(4,2), lwd=c(2,2))

pred.srv <- survey(fit.both$par)

plot(t,obs.srv, type="b", col=4, lwd=2, xlab="Year", ylab="Survey")
lines(t,pred.srv, type="b", col=2, lwd=2)
legend("topright", c("Observed","Predicted"), col=c(4,2), lwd=c(2,2))

# Why doesn’t it fit either as well as when you fit that data alone?

# d. Re-fit the data giving a weight of 10,000,000 to the survey SSQ. 
SSQ <- function(par){
	SSQ <- SSQ.catch(par) + 10000000 * SSQ.survey(par)
	return(SSQ)
}

fit.both2 <- optim(par, SSQ, method="L-BFGS-B", control=list(maxit=5000))

# Now what do you see? 
pred.catch <- catch(fit.both2$par)

plot(t,rowSums(obs.catch), type="b", col=4, lwd=2, xlab="Year", ylab="Catch", ylim=c(0,3000))
lines(t,rowSums(pred.catch), type="b", col=2, lwd=2)
legend("bottomright", c("Observed","Predicted"), col=c(4,2), lwd=c(2,2))

pred.srv <- survey(fit.both2$par)

plot(t,obs.srv, type="b", col=4, lwd=2, xlab="Year", ylab="Survey")
lines(t,pred.srv, type="b", col=2, lwd=2)
legend("topright", c("Observed","Predicted"), col=c(4,2), lwd=c(2,2))

# Explain your results.