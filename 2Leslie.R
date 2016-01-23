#############################################
#University of Alaska Fairbanks
#FISH 421: Fisheries Population Dynamics
#Lab 2: Leslie matrix models
#Instructor/lab author: Milo Adkison
#TA/code author: Phil Ganz
#############################################

# US Forest Service biologists recently discovered 
# a new species of mosquito-eating guppy (Poecilidae 
# chompem) in ponds near Sitka. (This may account for 
# the sunny disposition of the local residents). There 
# is some concern that these small, isolated populations 
# might be vulnerable to disturbance. 

# Field studies indicate that the age of first reproduction 
# in females is age 2, then fish skip a year before 
# reproducing again but then reproduce every year thereafter. 
# Reproduction is stressful and half of all females die after 
# reproducing. Survival and fecundity (divided by 2 – why?) 
# have the following rates:

# Survival at age
s <- c(0.7, 0.8, 0.5, 0.9, 0.5)

# Fecundity at age
m <- c(0, 0, 1, 0, 1)

# Arrange these values into a Leslie matrix
## First fill the diagonal of the matrix with
## survival values.
M <- diag(s[1:4], nrow=4, ncol=5)

## Then add fecundity to the top of the matrix
M <- rbind(m,M)
colnames(M) <- c(0:3, "4+")

## Add age 4+ survival
M[5,"4+"] <- s[5]

# 1. Project the growth of this population over 
# 50 years starting from 100 newborns. 
## Create empty matrix to hold values
N <- matrix(NA, nrow=50, ncol=5)

## Start with 100 newborns
N[1,] <- c(100, 0, 0, 0, 0)

## Project the population 50 years by multipling 
## the leslie matrix (M) by the abundance at age 
## in the previous year.
for (i in 2:NROW(N)){
	N[i,] <- M %*% N[i-1,]
}
N <- round(N)

## Add a column for total abundance
N <- cbind(N, Total=rowSums(N))

## Graph the total abundance through time.
plot(N[,"Total"], type="b", xlab="Year", ylab="Abundance")

## a. Describe the behavior in the first 10 years. 
##    Explain it.

## b. What is the intrinsic growth rate of this 
##    population (lambda)? 
## Approximate lambda by dividing the abundance
## in the last year by the abundace in the second-
## to-last year.
lambda <- N[50,"Total"]/N[49,"Total"]

## c.	What is the stable age distribution?
N[50,1:5]/N[50,"Total"]

## Plot abundance-at-age over time
plot(N[,1], type="b", col=1, xlab="Year", ylab="Abundance")
lines(N[,2], type="b", col=2)
lines(N[,3], type="b", col=3)
lines(N[,4], type="b", col=4)
lines(N[,5], type="b", col=5)
legend("topright", 
	legend=c("Age 0", "Age 1", "Age 2", "Age 3", "Age 4+"), 
	pch=1 , col=c(1:5))

# 2. Reproductive value is defined as the contribution 
# an individual makes to future generations relative to 
# the contribution a newborn makes. How large a population 
# in year 50 did we see when we started with 100 newborns? 
N[50,"Total"]

# How could we calculate the reproductive value of each of 
# the other ages? 

## a. Calculate the reproductive value of individuals of each age.
ev <- eigen(M)
U <- ev$vectors
V <- solve(Conj(U))
v <- abs(Re(V[1,]))
r.val <- v/v[1]

plot(0:4,r.val,pch=16,col="red",type="b", xlab="Age",
ylab="Reproductive Value")

# 3. These fish eat mosquitoes to acquire the protein necessary for 
# reproduction. Sitka’s mosquito control program has the potential to 
# reduce the number of mosquito larvae and thus reduce reproduction. 
# Another threat is stocking ponds with trout for recreational fisheries. 
# This would increase predation and reduce survival. 

## a. How important are each of these threats? Calculate rough values for 
## sensitivity and elasticity to changes in survival versus changes in 
## reproduction.

## Increase survival 10%
s_new <- s + 0.1

## Fill the diagonal of the matrix with the new
## survival values.
M_s <- diag(s_new[1:4], nrow=4, ncol=5)

## Then add fecundity to the top of the matrix
M_s <- rbind(m,M_s)
colnames(M_s) <- c(0:3, "4+")

## Add age 4+ survival
M_s[5,"4+"] <- s_new[5]

## Project abundance over 50 years, starting with 100 newborns.
N_s <- matrix(NA, nrow=50, ncol=5)
N_s[1,] <- c(100, 0, 0, 0, 0)
for (i in 2:NROW(N_s)){
	N_s[i,] <- M_s %*% N_s[i-1,]
}
N_s <- round(N_s)

## Add a column for total abundance
N_s <- cbind(N_s, Total=rowSums(N_s))

## Compute sensitivity and elasticity
Survival <- data.frame(Sensitivity=(N_s[50,5]/N_s[49,5]-lambda)/sum(s_new-s)) 
Survival$Elasticity <- Survival$Sensitivity * mean(s)/lambda

## Increase fecundity 10%
m_new <- m + 0.1

## Put new fecundity into the top of the matrix
M_m <- M
M_m["m",] <- m_new

## Project abundance over 50 years, starting with 100 newborns.
N_m <- matrix(NA, nrow=50, ncol=5)
N_m[1,] <- c(100, 0, 0, 0, 0)
for (i in 2:NROW(N_m)){
	N_m[i,] <- M_m %*% N_m[i-1,]
}
N_m <- round(N_m)

## Add a column for total abundance
N_m <- cbind(N_m, Total=rowSums(N_m))

## Compute sensitivity and elasticity
Fecundity <- data.frame(Sensitivity=(N_m[50,5]/N_m[49,5]-lambda)/sum(m_new-m)) 
Fecundity$Elasticity <- Fecundity$Sensitivity * mean(m)/lambda

## Plot sensitivity and elasticity by coefficient
barplot(c(Survival$Sensitivity, Survival$Elasticity, Fecundity$Sensitivity, Fecundity$Elasticity), xlab="Coefficient", col=c("darkblue","red"), beside=TRUE, space=rep(c(1,0),2))
axis(1, at=1:6, tick=FALSE, labels=c(NA,"Survival",NA,NA,"Fecundity",NA))
legend("topright", legend=c("Sensitivity", "Elasticity"), fill=c("darkblue","red"))

# 4. Modify this population model to incorporate a stochastic 
# density-dependent reproduction function. 

## a. A Beverton-Holt form could be added by dividing reproduction 
## by (1 + 0.0005*N). What effect does this have?

## Create empty matrix to hold values
N_bh <- matrix(NA, nrow=50, ncol=5)

## Start with 100 newborns
N_bh[1,] <- c(100, 0, 0, 0, 0)

## Project the population 50 years, dividing
## fecundity each year by (1 + 0.0005*N).
for (i in 2:NROW(N_bh)){
	M_bh <- M
	M_bh[1,] <- M_bh[1,]/(1 + 0.0005 * sum(N_bh[i-1,]))
	N_bh[i,] <- M_bh %*% N_bh[i-1,]
}
N_bh <- round(N_bh)

## Add a column for total abundance
N_bh <- cbind(N_bh, Total=rowSums(N_bh))

## Graph the total abundance through time.
plot(N[,"Total"], type="b", xlab="Year", ylab="Abundance")
lines(N_bh[,"Total"], type="b", col=2)
legend("topleft", legend=c("Original", "Beverton-Holt"), pch=1, col=c(1,2))

## b. A log-normal environmental effect could be added by multiplying 
## by exp(0.5Z), where Z is a normally-distributed random error. What 
## effect does this have?

## Create empty matrix to hold values
N_en <- matrix(NA, nrow=50, ncol=5)

## Start with 100 newborns
N_en[1,] <- c(100, 0, 0, 0, 0)

## Project the population 50 years, multpilying
## fecundity each year by exp(0.5Z).
set.seed(269)
for (i in 2:NROW(N_en)){
	M_en <- M
	M_en[1,] <- M_en[1,]*exp(0.5*rnorm(1))
	N_en[i,] <- M_en %*% N_en[i-1,]
}
N_en <- round(N_en)

## Add a column for total abundance
N_en <- cbind(N_en, Total=rowSums(N_en))

## Graph the total abundance through time.
plot(N[,"Total"], type="b", xlab="Year", ylab="Abundance", ylim=c(0,150))
lines(N_en[,"Total"], type="b", col=2)
legend("topleft", legend=c("Original", "Environmental flux"), pch=1, col=c(1,2))

## c. What other modifications might you make to this model?