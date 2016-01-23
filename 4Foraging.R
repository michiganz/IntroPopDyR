#############################################
#University of Alaska Fairbanks
#FISH 421: Fisheries Population Dynamics
#Lab 4: Foraging
#Instructor/lab author: Milo Adkison
#TA/code author: Phil Ganz
#############################################

# Today, we’ll look at when to leave a patch of food 
# to look for another patch, and how predation risk 
# affects foraging strategy. 

# The model consists of two types of food patches. One 
# has a total energy of E1, the other has a smaller 
# energy of E2. In either patch, the amount of energy 
# obtained is an asymptotic function of time spent in 
# the patch, as follows: e(t) = a*t/(1 + a*t/E). After 
# leaving a patch, the forager has a probability P1 of 
# encountering a patch of type 1 in the next time step 
# and a probability of P2 of encountering a patch of 
# type 2.

## Parameter values are:

E <- c(10, 5)

a <- c(2, 2)

p <- c(0.2, 0.2)

# 1. For each patch type, calculate the amount of energy 
# obtained by a forager spending 0,1,2,...,25 minutes in 
# the patch. Graph these relationships.

energy <- function(t,patch){
  i <- patch
  en <- a[i]*t/(1+a[i]*t/E[i])
  return(en)
}

t <- seq(0, 25)

e1 <- c(energy(t, 1))

e2 <- c(energy(t, 2))

plot(t, e1, type="l", col="Red", xlab="Time", ylab="Energy obtained")
lines(t, e2, col="Blue")
legend("topleft", inset = 0.05, c("patch 1", "patch 2"), lty = c(1, 1), col = c("Red", 
    "Blue"))

# 2. Make a guess as to the optimal amount of time a forager 
# should spend in each type of patch when it encounters it.

## We can solve this as an optimization problem. Write a function 
## that calculates the rate of energy intake for a 2-patch system 
## as a function of the time spent in each patch:
energy2 <- function(t){
	p1 <- p[1]/(p[1]+p[2])
	rate <- (p1*energy(t[1],1)+(1-p1)*energy(t[2],2))/(p1*t[1]+(1-p1)*t[2]+(1/(p[1]+p[2])))
	return(rate)
}

## Try it out!
t <- c(8,3)
energy2(t)

## Now optimize. By default, optim() minimizes, so by making 
## fnscale=-1, the function is instead maximized.
fit <- optim(t, energy2, control=list(fnscale=-1))

## Estimated optimum foraging times:
fit$par

# 3. Make a simulation model of your foraging strategy in which for 
# 50 time steps, you simulate: whether a patch ws found (pf), which 
# patch was found (wp), energy intake (en), foraging time (ft), cumulative 
# time (ct), and cumulative energy (ce). The probability of finding 
# a patch of some type is 1 – (1-p1)*(1-p2). 

## Time steps
ts <- 50

## So that you get the same random numbers as me:
set.seed(269)

## Construct the simulation as a function of the quantities of interest (t,E,a,p):
forage.sim <- function(t, E, a, p){
	# Establish vectors we will calculate values for
	pf <- NA
	wp <- NA
	ei <- NA
	ft <- NA
	for (i in 1:ts){ # For each iteration
	if (runif(1)<((1-p[1])*(1-p[2]))){pf[i] <- 0} # A patch is either not found
	else{pf[i] <- 1} # Or found
	if (pf[i]==0){wp[i] <- ""} # If (a patch isn't found), {leave Which_patch blank}
	if (pf[i]==1){ # If (a patch is found)...
		if (runif(1)<(p[1]/(p[1]+p[2]))){wp[i] <- 1} #...patch 1 is found with a probability of p[1]/(p[1]+p[2])
	else{wp[i]<-2}} #...otherwise, {the found patch must be patch 2}
	if (wp[i]==""){ei[i] <- 0} # If (no patch is found), {no energy is taken in}
	if (wp[i]==1){ei[i] <- a[1]*t[1]/(1+a[1]*t[1]/E[1])} # a[i]*t/(1+a[i]*t/E[i]) If (patch 1 is found), {an E[1] amount of energy is taken in}
	if (wp[i]==2){ei[i] <- a[2]*t[2]/(2+a[2]*t[2]/E[2])} # If (patch 2 is found), {an E[2] amount of energy is taken in}
	if (wp[i]==""){ft[i] <- 0} # If (no patch is found), {no time is spent foraging}
	if (wp[i]==1){ft[i] <- t[1]} # If (patch 1 is found), {an t[1] amount of time is spent foraging}
	if (wp[i]==2){ft[i] <- t[2]} # If (patch 2 is found), {an t[2] amount of time is spent foraging}
}
	ct <- ft[1] + 1
	ce <- ei[1]
	for (j in 2:ts){
		ct[j] <- ft[j] + 1 + ct[j-1]
		ce[j] <- ei[j] + ce[j-1]
	}
data.frame(
	Patch_Found=pf, 
	Which_patch=wp, 
	Energy_intake=ei, 
	Forage_time=ft, 
	Cumulative_time=ct, 
	Cumulative_energy=ce)
}

#Calculate the rate of energy intake obtained in your simulation.
## Write a short function that returns (Cumulative_energy/Cumulative_time):
energy.rate <- function(t, I){
	sim <- NA
	for (i in 1:I){
	sim[i] <- forage.sim(t, E, a, p)$Cumulative_energy[50]/forage.sim(t, E, a, p)$Cumulative_time[50]} 
	return(mean(sim))
}

# 4. Adjust the patch foraging times to find the best foraging strategy.

## Run 100 iterations of our 50 time step simulation
I <- 1000

## Test the optimum foraging times we found previously
t <- fit$par
f1 <- energy.rate(t, I)
f1

## Increase time in patch 1, decrease time in patch 2
t <- c(5,1)
f2 <- energy.rate(t, I)
f2

## Increase time in patch 2, decrease time in patch 1
t <- c(3,3)
f3 <- energy.rate(t, I)
f3

## Looks like the foraging times we found through optimization do maximize
## the rate of energy intake. Let's set foraging times to those.
t <- fit$par

# 5. Does the strategy change if the probability of encountering 
# a patch of type 2 increases to 0.5?
p <- c(0.2, 0.5)

f1 <- energy.rate(t,I)
f1

## Increase time in patch 1, decrease time in patch 2
t <- c(5,1)
f2 <- energy.rate(t,I)
f2

## Increase time in patch 2, decrease time in patch 1
t <- c(3,3)
f3 <- energy.rate(t,I)
f3

## It looks like the change in probability doesn't affect our strategy.

# 6. Add a per minute predation risk of 0.1 to patch 1. 
# That is, the probability of surviving t minutes in patch 1 is 0.9^t. 

## Return p to the original values:
p <- c(0.2, 0.2)

## Predation risks for each patch:
pr <- c(0.1,0)

## New function:
forage.sim.pr <- function(t,E,a,p,pr){
	# Establish vectors we will calculate values for
	pf <- NA
	wp <- NA
	ei <- NA
	ft <- NA
	### New quantities ###
	rand <- NA
	pred <- NA
	######################
	for (i in 1:ts){ # For each iteration
	if (runif(1)<((1-p[1])*(1-p[2]))){pf[i] <- 0} # A patch is either not found
	else{pf[i] <- 1} # Or found
	if (pf[i]==0){wp[i] <- ""} # If (a patch isn't found), {leave Which_patch blank}
	if (pf[i]==1){ # If (a patch is found)...
		if (runif(1)<(p[1]/(p[1]+p[2]))){wp[i] <- 1} #...patch 1 is found with a probability of p[1]/(p[1]+p[2])
	else{wp[i]<-2}} #...otherwise, {the found patch must be patch 2}
	if (wp[i]==""){ei[i] <- 0} # If (no patch is found), {no energy is taken in}
	if (wp[i]==1){ei[i] <- a[1]*t[1]/(1+a[1]*t[1]/E[1])} # If (patch 1 is found), {an E[1] amount of energy is taken in}
	if (wp[i]==2){ei[i] <- a[2]*t[2]/(2+a[2]*t[2]/E[2])} # If (patch 2 is found), {an E[2] amount of energy is taken in}
	if (wp[i]==""){ft[i] <- 0} # If (no patch is found), {no time is spent foraging}
	if (wp[i]==1){ft[i] <- t[1]} # If (patch 1 is found), {an t[1] amount of time is spent foraging}
	if (wp[i]==2){ft[i] <- t[2]} # If (patch 2 is found), {an t[2] amount of time is spent foraging}
	### New predation code here ###
	rand[i] <- runif(1)
	if (wp[i]==""){pred[i] <- 0} 
	if (wp[i]==1 & rand[i]<(1-(1-pr[1])^t[1])){pred[i] <- 1} 
	if (wp[i]==2 & rand[i]<(1-(1-pr[2])^t[2])){pred[i] <- 1} 
	if (wp[i]==1 & rand[i]>(1-(1-pr[1])^t[1])){pred[i] <- 0} 
	if (wp[i]==2 & rand[i]>(1-(1-pr[2])^t[2])){pred[i] <- 0} 
	###############################
}
	ct <- ft[1] + 1
	ce <- ei[1]
	for (j in 2:ts){
		ct[j] <- ft[j] + 1 + ct[j-1]
		ce[j] <- ei[j] + ce[j-1]
	}
data.frame(Patch_Found=pf, 
	Which_patch=wp, 
	Energy_intake=ei, 
	Forage_time=ft, 
	Cumulative_time=ct, 
	Cumulative_energy=ce,
	Random=rand,
	Predation=pred)
}

# Run a bunch of simulations with predation. Obviously, the fitness of a 
# forager eaten by a predator is zero. 

## Modify our previous simulation function so that rate of energy intake
## is 0 if predation occurs:
energy.rate.pr <- function(pr, I){
	all.vals <- vector("list") # Useful for creating a vector of data frames
	sim <- NA
	for (i in 1:I){
	all.vals[[i]] <- forage.sim.pr(t,E,a,p,pr) 
	if (sum(all.vals[[i]]$Predation)>0){sim[i] <- 0} # If (predation occurs){rate of energy intake is 0} 
	else{sim[i] <- all.vals[[i]]$Cumulative_energy[50]/all.vals[[i]]$Cumulative_time[50]} # Otherwise, it is calculated as before
	} 
	return(mean(sim))
}

## Run simulations with best foraging times found previously
t <- fit$par
f1 <- energy.rate.pr(pr,I)
f1

## Increase time in patch 1, decrease time in patch 2
t <- c(5,1)
f2 <- energy.rate.pr(pr,I)
f2

## Increase time in patch 2, decrease time in patch 1
t <- c(3,3)
f3 <- energy.rate.pr(pr,I)
f3

# How should you evaluate your foraging strategies now? What’s a good strategy?