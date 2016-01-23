#############################################
#University of Alaska Fairbanks
#FISH 421: Fisheries Population Dynamics
#Lab 9: Mixed Stocks
#Instructor/lab author: Milo Adkison
#TA/code author: Phil Ganz
#############################################

# Once upon a time in a lake there were three little 
# salmon substocks. The one with the biggest productivity 
# spawned in the outlet. The one with the medium productivity 
# spawned in the inlet. And the one with the littlest 
# productivity spawned on the lakeshore. They all produced 
# offspring according to the following equation:

# R = aS/(1+aS/b)

# 	          a	    b
# Outlet	  12	20,000
# Inlet	      8	    10,000
# Lakeshore	  4	    5,000

# 1. What is the unfished equilibrium spawning abundance for 
#    each stock?

# When a stock is unfished, no recruits are lost to fishing.
# If we define recruits as the number of fish that make it
# back to spawn after experiencing any sources of natural
# mortality, then recruits will be equal to spawners. We can 
# then solve for S:

# S = aS/(1+aS/b)
# 1 = a/(1+aS/b)  
# 1+aS/b = a
# aS/b = a-1
# aS = ab-b
# S = b-b/a

# We can then write a function to solve for the unfished stock
# size based on parameters a and b:

unfished.s <- function(a,b){
	S <- b-b/a
	return(S)
}

# Unfished spawning abundance for the outlet population:
out.f0 <- unfished.s(12,20000)

# Unfished spawning abundance for the inlet population:
in.f0 <- unfished.s(8,10000)

# Unfished spawning abundance for the lakeshore population:
shore.f0 <- unfished.s(4,5000)

# This fishery was managed by the Production Optimization 
# Research Committee under contract from ADF&G. The first PORC 
# biologist found out about three little substocks. He said, 
# “I’m going to harvest all of you to produce MSY from each”. 
# What MSY did he calculate for each, and in total? 

# Write a function for yield
yield <- function(S,a,b){
yield <- a*S/(1+a*S/b)-S
return(yield)
}

# Maximize yield() function to get MSY. Set a starting value for S:
S <- 5000

# MSY for the outlet population:
MSY.out <- optim(S, yield, a=12, b=20000, lower=0, upper=Inf, method="L-BFGS-B", control=list(fnscale=-1))

# MSY for the inlet population:
MSY.in <- optim(S, yield, a=8, b=10000, lower=0, upper=Inf, method="L-BFGS-B", control=list(fnscale=-1))

# MSY for the lakeshore population:
MSY.shore <- optim(S, yield, a=4, b=5000, lower=0, upper=Inf, method="L-BFGS-B", control=list(fnscale=-1))

# Combined MSY
MSY <- sum(MSY.out$value, MSY.in$value, MSY.shore$value)

# What was the combined escapement goal?
SMSY <- sum(MSY.out$par, MSY.in$par, MSY.shore$par)

# 2. Unfortunately, the stocks could only be harvested in a mixed 
# stock fishery. The first PORC biologist, although he diligently 
# monitored total escapement and achieved his combined goal every 
# year, was fired after a big, bad GreenyPeas biologist howled 
# about violating the ESA. 

# Do a 20 year projection to find out why. Start with the unfished 
# equilibrium abundance for each substock. Each year you’ll harvest 
# everything except the combined escapement goal. Escapement to each 
# substock will be proportional to its return relative to the returns 
# of the other substocks (why?).  

# Years of the simulation
years <- c(1:20)

# Write a function to return the number of recruits given S,a, and b:
recruits <- function(S,a,b){
R <- a*S/(1+a*S/b)
return(R)
}

# The parameters we have (1=Outlet, 2=Inlet, 3=Lakeshore):
par <- c(a1=12, b1=20000, 
		 a2=8, b2=10000, 
		 a3=4, b3=5000)

# What we want to simulate:
Outlet.Rec <- NA 	
Inlet.Rec <- NA  	
Lakeshore.Rec <- NA  	
O.Esc.Prop <- NA 	
I.Esc.Prop <- NA 
L.Esc.Prop <- NA 	
Outlet.Spawn <- NA 	
Inlet.Spawn <- NA 	
Lakeshore.Spawn <- NA 
Outlet.Harv <- NA 	
Inlet.Harv <- NA 
Lakeshore.Harv <- NA 	
Total.Harv <- NA 	
Outlet.Harv.Prop <- NA 	
Inlet.Harv.Prop <- NA 	
Lakeshore.Harv.Prop <- NA 

# Write a function to simulate the desired values:
sim.Q2 <- function(par){
# Initial recruits
Outlet.Rec[1] <- out.f0
Inlet.Rec[1] <- in.f0
Lakeshore.Rec[1] <- shore.f0
for (i in years){
for (j in 2:length(years)){
# Escapement proportion
O.Esc.Prop[i] <- Outlet.Rec[i]/sum(Outlet.Rec[i], Inlet.Rec[i], Lakeshore.Rec[i])
I.Esc.Prop[i] <- Inlet.Rec[i]/sum(Outlet.Rec[i], Inlet.Rec[i], Lakeshore.Rec[i])
L.Esc.Prop[i] <- Lakeshore.Rec[i]/sum(Outlet.Rec[i], Inlet.Rec[i], Lakeshore.Rec[i])
# Spawning stock
Outlet.Spawn[i] <- O.Esc.Prop[i]*SMSY
Inlet.Spawn[i] <- I.Esc.Prop[i]*SMSY
Lakeshore.Spawn[i] <- L.Esc.Prop[i]*SMSY
# Harvest
Outlet.Harv[i] <- Outlet.Rec[i]-Outlet.Spawn[i]
Inlet.Harv[i] <- Inlet.Rec[i]-Inlet.Spawn[i]
Lakeshore.Harv[i] <- Lakeshore.Rec[i]-Lakeshore.Spawn[i]
Total.Harv[i] <- sum(Outlet.Harv[i], Inlet.Harv[i], Lakeshore.Harv[i])
# Harvest proportion
Outlet.Harv.Prop <- Outlet.Harv[i]/Outlet.Rec[i]
Inlet.Harv.Prop <- Inlet.Harv[i]/Inlet.Rec[i]
Lakeshore.Harv.Prop <- Lakeshore.Harv[i]/Lakeshore.Rec[i]
# Recruits
Outlet.Rec[j] <- recruits(Outlet.Spawn[j-1], par["a1"], par["b1"])
Inlet.Rec[j] <- recruits(Inlet.Spawn[j-1], par["a2"], par["b2"])
Lakeshore.Rec[j] <- recruits(Lakeshore.Spawn[j-1], par["a3"], par["b3"])
}}
# Return objects as data frame
sim <- data.frame(Outlet.Rec,	
				Inlet.Rec, 	
				Lakeshore.Rec, 	
				O.Esc.Prop,	
				I.Esc.Prop,
				L.Esc.Prop,	
				Outlet.Spawn,	
				Inlet.Spawn,	
				Lakeshore.Spawn,
				Outlet.Harv,	
				Inlet.Harv,
				Lakeshore.Harv,	
				Total.Harv,	
				Outlet.Harv.Prop,	
				Inlet.Harv.Prop,	
				Lakeshore.Harv.Prop
)
return(sim)
}

pop.Q2 <- sim.Q2(par)

# How does the total harvest compare to the sum of the MSYs?
pop.Q2$Total.Harv[length(years)]
MSY

# How does the escapement to each substock compare to the escapement that produces MSY?
pop.Q2$Outlet.Spawn[length(years)]
MSY.out$par

pop.Q2$Inlet.Spawn[length(years)]
MSY.in$par

pop.Q2$Lakeshore.Spawn[length(years)]
MSY.shore$par

# 3. The Commissioner huffed and puffed and then shut down the fishery (annoying the fishermen), 
# and the stocks were left alone for many years. Finally, along came a second PORC biologist. 
# She said, “I’m going to take the mixed stock character of the fishery into account. I will 
# find the combined escapement goal that maximizes the yield from the mixed stock.” This she 
# did, and implemented her policy carefully and successfully. Unfortunately, her management 
# resulted in another successful ESA petition by the son of a wolf named Sierra (Sierra’s cub...), 
# and both she and the Commissioner were fired. 

# Do a 20 year projection to find out why. Start with the unfished equilibrium abundance for each substock. 
# Each year you’ll harvest everything except the combined escapement goal. Use optim() to find the combined 
# escapement goal that maximizes total harvest in the mixed stock fishery at equilibrium (i.e., in the 20th year).

# We are going to optimize to find a new mixed-stock SMSY, so we make two changes to our simulation function 
# (in the first and last lines):
# 1. Make sim.Q3 a function of SMSY instead of par
# 2. Return the last year of total harvest instead of the whole data frame, since that is what we'll optimize for
sim.Q3 <- function(SMSY){ # Changed #
# Initial recruits
Outlet.Rec[1] <- out.f0
Inlet.Rec[1] <- in.f0
Lakeshore.Rec[1] <- shore.f0
for (i in years){
for (j in 2:length(years)){
# Escapement proportion
O.Esc.Prop[i] <- Outlet.Rec[i]/sum(Outlet.Rec[i], Inlet.Rec[i], Lakeshore.Rec[i])
I.Esc.Prop[i] <- Inlet.Rec[i]/sum(Outlet.Rec[i], Inlet.Rec[i], Lakeshore.Rec[i])
L.Esc.Prop[i] <- Lakeshore.Rec[i]/sum(Outlet.Rec[i], Inlet.Rec[i], Lakeshore.Rec[i])
# Spawning stock
Outlet.Spawn[i] <- O.Esc.Prop[i]*SMSY
Inlet.Spawn[i] <- I.Esc.Prop[i]*SMSY
Lakeshore.Spawn[i] <- L.Esc.Prop[i]*SMSY
# Harvest
Outlet.Harv[i] <- Outlet.Rec[i]-Outlet.Spawn[i]
Inlet.Harv[i] <- Inlet.Rec[i]-Inlet.Spawn[i]
Lakeshore.Harv[i] <- Lakeshore.Rec[i]-Lakeshore.Spawn[i]
Total.Harv[i] <- sum(Outlet.Harv[i], Inlet.Harv[i], Lakeshore.Harv[i])
# Harvest proportion
Outlet.Harv.Prop <- Outlet.Harv[i]/Outlet.Rec[i]
Inlet.Harv.Prop <- Inlet.Harv[i]/Inlet.Rec[i]
Lakeshore.Harv.Prop <- Lakeshore.Harv[i]/Lakeshore.Rec[i]
# Recruits
Outlet.Rec[j] <- recruits(Outlet.Spawn[j-1], par["a1"], par["b1"])
Inlet.Rec[j] <- recruits(Inlet.Spawn[j-1], par["a2"], par["b2"])
Lakeshore.Rec[j] <- recruits(Lakeshore.Spawn[j-1], par["a3"], par["b3"])
}}
# Return objects as data frame
sim <- data.frame(Outlet.Rec,	
				Inlet.Rec, 	
				Lakeshore.Rec, 	
				O.Esc.Prop,	
				I.Esc.Prop,
				L.Esc.Prop,	
				Outlet.Spawn,	
				Inlet.Spawn,	
				Lakeshore.Spawn,
				Outlet.Harv,	
				Inlet.Harv,
				Lakeshore.Harv,	
				Total.Harv,	
				Outlet.Harv.Prop,	
				Inlet.Harv.Prop,	
				Lakeshore.Harv.Prop
)
return(sim$Total.Harv[length(years)]) # Changed #
}

# Now optimize to find the new SMSY...
SMSY.Q3 <- optim(SMSY, sim.Q3, method="BFGS", control=list(fnscale=-1))

# ...and use that SMSY in the new simulation
SMSY <- SMSY.Q3$par
pop.Q3 <- sim.Q2(par)

# How does the escapement to each substock compare to the escapement that produces MSY? 
pop.Q3$Outlet.Spawn[length(years)]
MSY.out$par

pop.Q3$Inlet.Spawn[length(years)]
MSY.in$par

pop.Q3$Lakeshore.Spawn[length(years)]
MSY.shore$par

# How does the escapement to each substock compare to the unfished escapement?
pop.Q3$Outlet.Spawn[length(years)]
out.f0

pop.Q3$Inlet.Spawn[length(years)]
in.f0

pop.Q3$Lakeshore.Spawn[length(years)]
shore.f0

# 4. Graph the stock recruitment relationship for all three stocks for escapements 
# in the range of 0-30,000. Add to the graph replacement lines for a 50%, 67%, and 
# 80% harvest rate. How does this help explain why the biologists were fired?
S <- seq(0,30000,by=1000)

outlet <- recruits(S, par["a1"], par["b1"])
inlet <- recruits(S, par["a2"], par["b2"])
lakeshore <- recruits(S, par["a3"], par["b3"])

plot(S, outlet, type="l", col=2, lwd=2, xlab="Stock", ylab="Recruitment")
lines(S, inlet, type="l", col=3, lwd=2)
lines(S, lakeshore, type="l", col=4, lwd=2)
lines(S, 2*S, type="l", lty=3, lwd=2)
lines(S, 3*S, type="l", lty=2, lwd=2)
lines(S, 5*S, type="l", lty=1, lwd=2)
legend(20000, 17000, c("outlet","inlet","lakeshore","50% harvest","67% harvest","80% harvest"), col=c(2,3,4,1,1,1), lwd=c(2,2,2,1,1,1), lty=c(1,1,1,3,2,1))

# 5. Graph R/S vs S for each stock for escapements (S) in the range 0-30,000. 
# Add horizontal lines to these graphs at R/S = 1, 2, 3 and 4. 
par(mfrow=c(1,3))

plot(S, outlet/S, type="l", col=2, lwd=2, main="Outlet", xlab="Stock", ylab="Recruits per spawner", ylim=c(0,6))
abline(h=1, lty=2)
abline(h=2, lty=2)
abline(h=3, lty=2)
abline(h=4, lty=2)

plot(S, inlet/S, type="l", col=3, lwd=2, main="Inlet", xlab="Stock", ylab="", ylim=c(0,6))
abline(h=1, lty=2)
abline(h=2, lty=2)
abline(h=3, lty=2)
abline(h=4, lty=2)

plot(S, lakeshore/S, type="l", col=4, lwd=2, main="Lakeshore", xlab="Stock", ylab="", ylim=c(0,6))
abline(h=1, lty=2)
abline(h=2, lty=2)
abline(h=3, lty=2)
abline(h=4, lty=2)

# What does each of these lines tell you?

# 6. You’ve just been hired as the third PORC biologist. This is absolutely the 
# last straw. The NGO’s are throwing sticks and the fishers are throwing bricks 
# at ADF&G and PORC. It’s up to you – how will you manage this mixed stock fishery?