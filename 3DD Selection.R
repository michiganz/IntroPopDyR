#############################################
#University of Alaska Fairbanks
#FISH 421: Fisheries Population Dynamics
#Lab 3: r vs. K selection
#Instructor/lab author: Milo Adkison
#TA/code author: Phil Ganz
#############################################

# (adapted from J. Roughgarden’s 1979 book)

# Let’s refine the logistic model of population 
# growth to include natural selection and evolution. 

# We’ll start with the population model:
# Nt+1 = Wavg Nt

# Where Wavg is the average fitness of an individual 
# in the population.

# The parameters of the logistic model are affected by 
# the genotype. For an individual of genotype AA, the 
# fitness at time t is:
# WAA = 1+rAA*(1-Nt/KAA)

# The fitness of the other genotypes is similar:
# WAa = 1+rAa*(1-Nt/KAa)
# Waa = 1+raa*(1-Nt/Kaa)

# The average fitness is a function of the genotype 
# frequencies. Assuming H-W equilibrium, the average 
# fitness as a function of pt, the frequency of the 
# A allele, is:
# pt*pt*WAA + 2*pt*(1-pt)*WAa +(1-pt)*(1-pt)*Waa

# and selection will change this frequency as follows:
# pt+1 = [pt*WAA +(1-pt)*(WAa)]*pt/Wavg

# Optimally, populations would evolve r’s and K’s as 
# large as possible. At some point, physiological limits 
# will be reached and there will be a tradeoff where 
# increasing carrying capacity (K) will be at the expense 
# of a slower growth rate (r). In this lab, we will 
# examine the optimal tradeoff between these two 
# parameters, and how environmental disturbance can 
# affect the optimal strategy.

# Step 1. Create the model.

## Parameters:
rAA <- 0.8
rAa <- 0.8
raa <- 0.6
KAA <- 8000
KAa <- 8000
Kaa <- 12000

## Population

## Write a function for projecting the population based on
## p1 and N1 values
project <- function(p1, N1){
	Year <- 1:100
	p <- p1
	N <- N1
	AA <- NA
	Aa <- NA
	aa <- NA
	WAA <- NA
	WAa <- NA
	Waa <- NA
	Wavg <- NA
		for (i in Year){
		AA[i] <- p[i] * p[i]
		Aa[i] <- 2 * p[i] * (1-p[i])
		aa[i] <- (1-p[i]) * (1-p[i])
		WAA[i] <- 1+rAA*(1-N[i]/KAA)
		WAa[i] <- 1+rAa*(1-N[i]/KAa)
		Waa[i] <- 1+raa*(1-N[i]/Kaa)
		Wavg[i] <- p[i]*p[i]*WAA[i] + 2*p[i]*(1-p[i])*WAa[i] +(1-p[i])*(1-p[i])*Waa[i]
		p[i+1] <- (p[i]*WAA[i] +(1-p[i])*(WAa[i]))*p[i]/Wavg[i]
		N[i+1] <- N[i] * Wavg[i]
	}
	round(data.frame(Year, p=p[-101], N=N[-101], AA, Aa, aa, WAA, WAa, Waa, Wavg),2)
}

## Try it out!
pop <- project(0.5, 1000)

plot(pop$p, pop$N, type="b", xlab="p", ylab="N")
points(pop$p[1], pop$N[1], col=2, pch=19) # Color and fill first point

# Step 2. Examine the effect of initial conditions on the outcome.
# Try a variety of different combinations of p1 and N1.
project(0.25, 1000)[100,]
project(0.75, 1000)[100,]
project(0.5, 200)[100,]
project(0.5, 4000)[100,]
project(0.25, 4000)[100,]
project(0.75, 200)[100,]

# Do the initial conditions have any effect of the outcome of selection? 

# Is a higher r or a higher K the favored option?

## Set Ks equal, increase raa
Kaa <- 8000
raa <- 1
project(0.5, 1000)[100,]

## Set rs equal, increase Kaa
Kaa <- 12000
raa <- 0.8
project(0.5, 1000)[100,]

## Return raa to its original value
raa <- 0.6

# Step 3. Change the dominant allele.
# Adjust the values of rAa and KAa to make a the dominant allele.

## Make rAa=raa and KAa=Kaa.
rAa <- 0.6
KAa <- 12000

# Does this change the outcome of selection?
project(0.5, 1000)[100,]

# Does this change the influence of the initial conditions? 
project(0.25, 1000)[100,]
project(0.75, 1000)[100,]
project(0.5, 200)[100,]
project(0.5, 4000)[100,]
project(0.25, 4000)[100,]
project(0.75, 200)[100,]

# Is a higher r or a higher K the favored option?

# Step 4. Heterozygote advantage. 
# Change rAa to 0.7 and KAa to 15,000. 
rAa <- 0.7
KAa <- 15000

# Does this change the outcome of selection? 
project(0.5, 1000)[100,]

# Does this change the influence of the initial conditions? 
project(0.25, 1000)[100,]
project(0.75, 1000)[100,]
project(0.5, 200)[100,]
project(0.5, 4000)[100,]
project(0.25, 4000)[100,]
project(0.75, 200)[100,]

# Is a higher r or a higher K the favored option?

# Step 5. Heterozygote disadvantage. 
# Change rAa to 0.7 and KAa to 5,000. 
KAa <- 5000

# Does this change the outcome of selection? 
project(0.5, 1000)[100,]

# Does this change the influence of the initial conditions? 
project(0.25, 1000)[100,]
project(0.75, 1000)[100,]
project(0.5, 200)[100,]
project(0.5, 4000)[100,]
project(0.25, 4000)[100,]
project(0.75, 200)[100,]

# Is a higher r or a higher K the favored option?

# Step 6. Add sporadic disturbance. 
# We’ll now modify the model to occasionally kill off half of the population.

# Add a new parameter to your called “Disturbance frequency” ("D" in the new function). 
# Set its value to 0.4, which means that 40% of the time we’ll have a disturbance 
# that kills half of the population. (I encourage you to try other values, too).
project.d <- function(p1, N1, D){
	Year <- 1:100
	p <- p1
	N <- N1
	AA <- NA
	Aa <- NA
	aa <- NA
	WAA <- NA
	WAa <- NA
	Waa <- NA
	Wavg <- NA
	rand <- NA
		for (i in Year){
		AA[i] <- p[i] * p[i]
		Aa[i] <- 2 * p[i] * (1-p[i])
		aa[i] <- (1-p[i]) * (1-p[i])
		WAA[i] <- 1+rAA*(1-N[i]/KAA)
		WAa[i] <- 1+rAa*(1-N[i]/KAa)
		Waa[i] <- 1+raa*(1-N[i]/Kaa)
		Wavg[i] <- p[i]*p[i]*WAA[i] + 2*p[i]*(1-p[i])*WAa[i] +(1-p[i])*(1-p[i])*Waa[i]
		p[i+1] <- (p[i]*WAA[i] +(1-p[i])*(WAa[i]))*p[i]/Wavg[i]
		N[i+1] <- N[i] * Wavg[i]
		rand[i] <- runif(1) # Generate 1 random number (between 0 and 1 by default) from the uniform distribution
		if (rand[i]<D){ # If the random number is < our disturbance frequency
			N[i+1] <- N[i+1]/2 # Divide abundance by 2
		}
		else{ # Otherwise
			N[i+1] <- N[i+1] # Leave it as is
		}
	}
	round(data.frame(Year, p=p[-101], N=N[-101], AA, Aa, aa, WAA, WAa, Waa, Wavg, rand),2)
}

# Try it out!
set.seed(269) # Assures that you get the same random numbers as me
pop.d <- project.d(0.5, 1000, 0.4)

plot(pop.d$p, pop.d$N, type="b", xlab="p", ylab="N")
points(pop.d$p[1], pop.d$N[1], col=2, pch=19) # Color and fill first point

# How often is the random numbers added to our code less than the disturbance frequency (D)? 
sum(pop.d$rand<0.4)/length(pop.d$rand)

# What happens when it is? What happens when it’s greater than the disturbance frequency?