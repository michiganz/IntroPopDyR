#############################################
#University of Alaska Fairbanks
#FISH 421: Fisheries Population Dynamics
#Lab 7: Fleet Dynamics
#Instructor/lab author: Milo Adkison
#TA/code author: Phil Ganz
#############################################

# Objectives: 
# 
# 1. To revisit the difference in equilibrium states 
# between an open-access fishery and one operated to 
# maximize economic yield.

# 2.To explore the dynamics of a fish stock and fishing 
# fleet interacting far from equilibrium.

# B_t+1 = B_t + r * B_t * (1-B_t/K) - q * V_t * B_t

# B_t = biomass in year t
# r, K = parameters of Schaefer model
# q = vessel catch efficiency
# V_t = number of fishing vessels in year t

# P_t = p * C_t - c * V_t 
#     = p * q * V_t * B_t - c * V_t 

# P_t = total fisheries profit in year t = revenues – costs 
# p = price per unit of catch
# C_t = catch in year t
# c = annual vessel costs

# 1. Make a model of this fishery. 
# Use the following values: r=0.5, K=1000, q=0.001, 
# p=1.5, c=0.5, B_1=K, and V=100. 

r <- 0.5
K <- 1000
q <- 0.001
p <- 1.5
c <- 0.5
B <- K 
V <- 100

par <- c(r=r, K=K, q=q, p=p, c=c, B=K, V=V)

# Use the following column headings: Year, Vessels, 
# Biomass, Catch, Biomass After Fishing, Production, 
# Profit, and Profit per Vessel. Simulate 100 years 
# of the fishery.

fleet <- data.frame(
	Year=1:100, 
	Vessels=NA, 
	Biomass_start=NA, 
	Catch=NA, 
	Biomass_end=NA, 
	Production=NA,
	Profit=NA,
	Vessel_profit=NA)

## Construct the fishery as a function of our parameters
fleet.sim <- function(par){
	fleet$Biomass_start[1] <- par["K"]
	for (i in fleet$Year){
		for (j in 2:length(fleet$Year)){
		fleet$Vessels[i] <- par["V"]
		fleet$Catch[i] <- par["q"] * par["V"] * fleet$Biomass_start[i]
		fleet$Biomass_end[i] <- fleet$Biomass_start[i] - fleet$Catch[i]
		fleet$Production[i] <- par["r"] * fleet$Biomass_end[i] * (1 - fleet$Biomass_end[i]/par["K"])
		fleet$Biomass_start[j] <- fleet$Biomass_start[j-1] - fleet$Catch[j-1] + fleet$Production[j-1]
		fleet$Profit[i] <- par["p"] * fleet$Catch[i] - par["c"] * fleet$Vessels[i]
		fleet$Vessel_profit[i] <- fleet$Profit[i]/fleet$Vessels[i]
	}}
	return(fleet)
}

## Use this function to fill in values of empty data frame:
fleet <- fleet.sim(par)

# 2. Graph catch and biomass vs. year.
plot(fleet$Year, fleet$Biomass_start, type="b", col=2, xlab="Year", ylab="Biomass")
plot(fleet$Year, fleet$Catch, type="b", col=4, xlab="Year", ylab="Catch")

# 3. Find the stock biomass and the number of vessels that 
# produces the maximum economic yield (MEY). Calculate the 
# number of vessels that produce the MSY. How do these two 
# situations compare?

## Rewrite our fleet.sim() function so that it's a function of V
## and returns profit in the last year so that that value can be 
## maximized
fleet.mey <- function(V){
	fleet$Biomass_start[1] <- K
	for (i in fleet$Year){
		for (j in 2:length(fleet$Year)){
		fleet$Vessels[i] <- V
		fleet$Catch[i] <- q * V * fleet$Biomass_start[i]
		fleet$Biomass_end[i] <- fleet$Biomass_start[i] - fleet$Catch[i]
		fleet$Production[i] <- r * fleet$Biomass_end[i] * (1 - fleet$Biomass_end[i]/K)
		fleet$Biomass_start[j] <- fleet$Biomass_start[j-1] - fleet$Catch[j-1] + fleet$Production[j-1]
		fleet$Profit[i] <- p * fleet$Catch[i] - c * fleet$Vessels[i]
		fleet$Vessel_profit[i] <- fleet$Profit[i]/fleet$Vessels[i]
	}}
	return(fleet$Profit[length(fleet$Year)])
}

## Optimize this new function
fit.mey <- optim(V, method="L-BFGS-B", lower=0, upper=Inf, fleet.mey, control=list(fnscale=-1))

## Set vessel number to the number estimated for MEY and use our previous
## function to simulate the trend in biomass and other values.
par["V"] <- fit.mey$par
fleet.mey <- fleet.sim(par)

## Now use the same approach to estimate the number of vessels that produce
## MSY
fleet.msy <- function(V){
	fleet$Biomass_start[1] <- K
	for (i in fleet$Year){
		for (j in 2:length(fleet$Year)){
		fleet$Vessels[i] <- V
		fleet$Catch[i] <- q * V * fleet$Biomass_start[i]
		fleet$Biomass_end[i] <- fleet$Biomass_start[i] - fleet$Catch[i]
		fleet$Production[i] <- r * fleet$Biomass_end[i] * (1 - fleet$Biomass_end[i]/K)
		fleet$Biomass_start[j] <- fleet$Biomass_start[j-1] - fleet$Catch[j-1] + fleet$Production[j-1]
		fleet$Profit[i] <- p * fleet$Catch[i] - c * fleet$Vessels[i]
		fleet$Vessel_profit[i] <- fleet$Profit[i]/fleet$Vessels[i]
	}}
	return(fleet$Catch[length(fleet$Year)])
}

fit.msy <- optim(V, method="L-BFGS-B", lower=0, upper=Inf, fleet.msy, control=list(fnscale=-1))

par["V"] <- fit.msy$par
fleet.msy <- fleet.sim(par)

## Compare trend in biomass with MEY number of vessels (160) and MSY number of vessels (200)
plot(fleet.msy$Year, fleet.msy$Biomass_start, type="b", col=2, xlab="Year", ylab="Biomass", ylim=c(0,1000))
lines(fleet.mey$Year, fleet.mey$Biomass_start, type="b", col=4)
legend("topright", legend=c("MEY", "MSY"), lty=c(1,1), col=c(4,2))

## 4. Find the stock biomass and the number of vessels that would result 
## from an open access fishery. How does this differ from the MEY case?
fleet.oa <- function(V){
	fleet$Biomass_start[1] <- K
	for (i in fleet$Year){
		for (j in 2:length(fleet$Year)){
		fleet$Vessels[i] <- V
		fleet$Catch[i] <- q * V * fleet$Biomass_start[i]
		fleet$Biomass_end[i] <- fleet$Biomass_start[i] - fleet$Catch[i]
		fleet$Production[i] <- r * fleet$Biomass_end[i] * (1 - fleet$Biomass_end[i]/K)
		fleet$Biomass_start[j] <- fleet$Biomass_start[j-1] - fleet$Catch[j-1] + fleet$Production[j-1]
		fleet$Profit[i] <- p * fleet$Catch[i] - c * fleet$Vessels[i]
		fleet$Vessel_profit[i] <- fleet$Profit[i]/fleet$Vessels[i]
	}}
	return(fleet$Vessel_profit[length(fleet$Year)]) #How does this line compare to that of fleet.mey()? 
}

fit.oa <- optim(V, method="L-BFGS-B", lower=0, upper=Inf, fleet.oa)

par["V"] <- fit.oa$par
fleet.oa <- fleet.sim(par)

## Compare open access biomass trend to that of MEY and MSY graphically
plot(fleet.msy$Year, fleet.msy$Biomass_start, type="b", col=2, xlab="Year", ylab="Biomass", ylim=c(0,1000))
lines(fleet.mey$Year, fleet.mey$Biomass_start, type="b", col=4)
lines(fleet.oa$Year, fleet.oa$Biomass_start, type="b", col=3)
legend("topright", legend=c("MEY", "MSY", "Open access"), lty=c(1,1,1), col=c(4,2,3))

# Now let’s include a dynamic response by fishers to the profitability 
# of the fishery. When excess profits are being made, fishers will enter. 
# When costs exceed revenues, some will leave.

# V_t+1 = V_t + P_t/d

# d = cost of a new vessel

# Note that this equation implies that all excess profits are invested in 
# new fishing boats, and that all losses translate to an equivalent loss 
# of boats. Similar equations can be derived from less restrictive assumptions.

# 5. Change the “Vessels” column to follow this dynamical equation. Assume 
# the number of vessels in year 1 is 1. Set d=2. What happens to catch and 
# biomass over time?

par <- c(r=r, K=K, q=q, p=p, c=c, B=K, V=1, d=2)

fleet.sim.d <- function(par){
	fleet$Biomass_start[1] <- par["K"]
	fleet$Vessels[1] <- par["V"] ### New code ###
	for (i in fleet$Year){
		for (j in 2:length(fleet$Year)){
		fleet$Catch[i] <- par["q"] * fleet$Vessels[i] * fleet$Biomass_start[i]
		fleet$Biomass_end[i] <- fleet$Biomass_start[i] - fleet$Catch[i]
		fleet$Production[i] <- par["r"] * fleet$Biomass_end[i] * (1 - fleet$Biomass_end[i]/par["K"])
		fleet$Biomass_start[j] <- fleet$Biomass_start[j-1] - fleet$Catch[j-1] + fleet$Production[j-1]
		fleet$Profit[i] <- par["p"] * fleet$Catch[i] - par["c"] * fleet$Vessels[i]
		fleet$Vessel_profit[i] <- fleet$Profit[i]/fleet$Vessels[i]
		fleet$Vessels[j] <- fleet$Vessels[j-1] + fleet$Profit[j-1]/par["d"]
	}}
	return(fleet)
}

fleet.d <- fleet.sim.d(par)

# 6. Graph vessels and profit vs. year (use a secondary axis for profits).
par(mar=c(5,4,4,5)+.1)
plot(fleet.d$Year, fleet.d$Profit, type="l", lwd=3, col=2, xaxt="n", yaxt="n", xlab="Year", ylab="Profit")
par(new=TRUE)
plot(fleet.d$Year, fleet.d$Vessels, type="l", lwd=3, col=4, xlab="", ylab="")
axis(4)
mtext("Vessels",side=4,line=3)
legend("topright",col=c(2, 4), lty=1, legend=c("Profit","Vessels"))

# What pattern do you see? Does the system come to the open access equilibrium?

# 7. Graph biomass vs. vessels. What pattern do you see? Is the system coming to 
# equilibrium?
plot(fleet.d$Vessels, fleet.d$Biomass_start, type="l", lwd=3, col=4, xlab="Year", ylab="Biomass")

# 8. Change d from 2 to 5. What happens to the dynamics of the fishery? Why does 
# increasing the cost of boats stabilize the dynamics?
par["d"] <- 5

fleet.d <- fleet.sim.d(par)

par(mar=c(5,4,4,5)+.1)
plot(fleet.d$Year, fleet.d$Profit, type="l", lwd=3, col=2, xaxt="n", yaxt="n", xlab="Year", ylab="Profit")
par(new=TRUE)
plot(fleet.d$Year, fleet.d$Vessels, type="l", lwd=3, col=4, xlab="", ylab="")
axis(4)
mtext("Vessels",side=4,line=3)
legend("topright",col=c(2, 4), lty=1, legend=c("Profit","Vessels"))

plot(fleet.d$Vessels, fleet.d$Biomass_start, type="l", lwd=3, col=4, xlab="Year", ylab="Biomass")

# 9. What would happen if the government wanted to encourage development of this 
# untapped fish resource and subsidized vessel construction so that the boat costs 
# were reduced to 1?
par["d"] <- 1

fleet.d <- fleet.sim.d(par)

par(mar=c(5,4,4,5)+.1)
plot(fleet.d$Year, fleet.d$Profit, type="l", lwd=3, col=2, xaxt="n", yaxt="n", xlab="Year", ylab="Profit", ylim=c(0,1200), xlim=c(0,60))
par(new=TRUE)
plot(fleet.d$Year, fleet.d$Vessels, type="l", lwd=3, col=4, xlab="", ylab="")
axis(4)
mtext("Vessels",side=4,line=3)
legend("topright",col=c(2, 4), lty=1, legend=c("Profit","Vessels"))

plot(fleet.d$Vessels, fleet.d$Biomass_start, type="l", lwd=3, col=4, xlab="Year", ylab="Biomass", ylim=c(0,1200), xlim=c(0,1200))

# 10. Set boat prices back to 2. What would happen if the price of fish dropped from 
# 1.5 to 1? 
par["d"] <- 2
par["p"] <- 1

fleet.d <- fleet.sim.d(par)

par(mar=c(5,4,4,5)+.1)
plot(fleet.d$Year, fleet.d$Profit, type="l", lwd=3, col=2, xaxt="n", yaxt="n", xlab="Year", ylab="Profit")
par(new=TRUE)
plot(fleet.d$Year, fleet.d$Vessels, type="l", lwd=3, col=4, xlab="", ylab="")
axis(4)
mtext("Vessels",side=4,line=3)
legend("topright",col=c(2, 4), lty=1, legend=c("Profit","Vessels"))

plot(fleet.d$Vessels, fleet.d$Biomass_start, type="l", lwd=3, col=4, xlab="Year", ylab="Biomass")

# Why does this happen?