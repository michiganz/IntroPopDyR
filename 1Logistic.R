#############################################
#University of Alaska Fairbanks
#FISH 421: Fisheries Population Dynamics
#Lab 1: Logistic population growth
#Instructor/lab author: Milo Adkison
#TA/code author: Phil Ganz
#############################################

# Step 1. Project the growth of the population over 100 years from 
#         a small initial abundance

## 1b. Enter the two parameters of the model

## Intrinsic population growth rate at low abundance
r<-0.2
## Carrying capacity
K<-1000

## 1c. Calculate the abundance of a population that starts with 1 
## individual over a period of 100 years 
Year<-seq(1:100)

## Create vector to hold abundance (N) values
N<-numeric(100)
## In year 1, the abundance will be 1
N[1]<-1
## Use logistic growth equation to calculate remaining abundances
for (t in 1:(max(Year)-1)){
N[t+1]<-N[t]+r*N[t]*(1-N[t]/K)}

# Step 2. Graph the abundance over time
plot(N~Year, type="l")

## Change the initial abundance (in year 1) to 1500. What does the population do?
N[1]<-1500
for (t in 1:(max(Year)-1)){
N[t+1]<-N[t]+r*N[t]*(1-N[t]/K)}

plot(N~Year, type="l")

# Step 3. Graph N(t+1) for values of N(t) between 0 and 1500.

## 3a. In the vector N(t), place numbers ranging from 0 to 1500 in 
## increments of 15.
N_t<-seq(0,1500,by=15)

## In the N(t+1) vector, use the logistic formula to calculate what 
## abundance the population would have if it started from the value 
## in vector N(t) the year before. 
N_t1<-numeric(length(N_t))
for (t in 1:length(N_t)){
N_t1[t]<-N_t[t]+r*N_t[t]*(1-N_t[t]/K)}

## 3b. Make a scatterplot graph of N(t+1) with N(t) as the X value.
plot(N_t1~N_t, type ="l", col="Blue", ylim=c(0,1200), ylab="N_t+1")

## 3c. Graph the replacement line, the 45 degree line which shows 
## N(t+1) equal to N(t). 
lines(N_t,N_t, type="l", col="Red")

# Step 4. Change the parameters.

## What do you think will happen if you change K? Try it. Try several values.
K<-c(100,500,1000,1500,2000,5000)

## Set up empty matrix to hold values of N for different K's.
N<-matrix(NA, nrow=length(Year), ncol=(length(K)))

## 4a. Change the initial abundance to 1. 
## Brackets = [row,column], so N[1,] indicates that abundance is 1 for the 
## first row and EVERY coulumn, i.e. the first year, every experimental K  
N[1,]<-1

## Use logistic growth equation to calculate remaining abundances for each K
for (k in 1:length(K)){
for (t in 1:(max(Year)-1)){
N[t+1,k]<-N[t,k]+r*N[t,k]*(1-N[t,k]/K[k])}}

## Create multipanel plot for comparison
par(mfrow=c(2,3))
for (k in 1:length(K)){
plot(N[,k]~Year, type="l", ylim=c(0,5000), ylab="N", main=paste("K=",K[k]))}

## After you’ve tried a few, change K back to 1000.
K<-1000

## 4b. What do you think will happen if we double ‘r’ from 0.2 to 0.4? 
## Think about this for a few minutes before you try it.  
## Now try it.
r<-c(0.2,0.4)

## Set up empty matrix to hold values of N for different r's.
N<-matrix(NA, nrow=length(Year), ncol=(length(r)))

## Change the initial abundance to 1.  
N[1,]<-1

## Use logistic growth equation to calculate remaining abundances for each r
for (i in 1:length(r)){
for (t in 1:(max(Year)-1)){
N[t+1,i]<-N[t,i]+r[i]*N[t,i]*(1-N[t,i]/K)}}

## Create multipanel plot for comparison
par(mfrow=c(1,2))
for (i in 1:length(r)){
plot(N[,i]~Year, type="l", ylim=c(0,1000), ylab="N", main=paste("r=",r[i]))}

## What changed? What didn’t change? Why? 

## Look at the graph of N(t+1) vs. N(t). How does it change when ‘r’ is changed 
## from 0.2 to 0.4. What’s the change in the amount of density-dependence?
N_t1<-matrix(NA, nrow=length(N_t), ncol=(length(r)))

## In the r_0.2 and r_0.4 columns, use the logistic formula to calculate what 
## abundance the population would have if it started from the value 
## in column N(t) the year before and had growth rate r.
for (i in 1:length(r)){
for (t in 1:length(N_t)){
N_t1[t,i]<-N_t[t]+r[i]*N_t[t]*(1-N_t[t]/K)}}

## Look at the graph of N(t+1) vs. N(t). How does it 
## change when ‘r’ is changed from 0.2 to 0.4. What’s 
## the change in the amount of density-dependence?
par(mfrow=c(1,2))
for (i in 1:length(r)){
plot(N_t1[,i]~N_t, type ="l", col="Blue", ylim=c(0,1200), ylab="N_t+1", main=paste("r=",r[i]))
lines(N_t,N_t, type="l", col="Red")}

## 4c. Now try ‘r’ = 1.0. What happened? 
r<-c(0.2,0.4,1.0)

## Abundance vs. time
N<-matrix(NA, nrow=length(Year), ncol=(length(r)))
N[1,]<-1
for (i in 1:length(r)){
for (t in 1:(max(Year)-1)){
N[t+1,i]<-N[t,i]+r[i]*N[t,i]*(1-N[t,i]/K)}}

## N_t+1 vs. N_t
N_t1<-matrix(NA, nrow=length(N_t), ncol=(length(r)))
for (i in 1:length(r)){
for (t in 1:length(N_t)){
N_t1[t,i]<-N_t[t]+r[i]*N_t[t]*(1-N_t[t]/K)}}

par(mfrow=c(2,3))
for (i in 1:length(r)){
plot(N[,i]~Year, type="l", ylim=c(0,1000), ylab="N", main=paste("r=",r[i]))}
for (i in 1:length(r)){
plot(N_t1[,i]~N_t, type ="l", col="Blue", ylim=c(0,1200), ylab="N_t+1", main=paste("r=",r[i]))
lines(N_t,N_t, type="l", col="Red")}

## 4d. Try ‘r’ = 2.0. What’s that funny little squiggle there starting about 
## year 15? Look at the abundance for years 8-20. Is the population climbing 
## steadily to its carrying capacity? What is it doing?
r<-c(2.0)

## Abundance vs. time
N<-matrix(NA, nrow=length(Year), ncol=(length(r)))
N[1,]<-1
for (i in 1:length(r)){
for (t in 1:(max(Year)-1)){
N[t+1,i]<-N[t,i]+r[i]*N[t,i]*(1-N[t,i]/K)}}
 
## N_t+1 vs. N_t
N_t1<-matrix(NA, nrow=length(N_t), ncol=(length(r)))
for (i in 1:length(r)){
for (t in 1:length(N_t)){
N_t1[t,i]<-N_t[t]+r[i]*N_t[t]*(1-N_t[t]/K)}}

par(mfrow=c(1,2))
for (i in 1:length(r)){
plot(N[,i]~Year, type="l", ylim=c(0,1200), ylab="N", main=paste("r=",r[i]))}
for (i in 1:length(r)){
plot(N_t1[,i]~N_t, type ="l", col="Blue", ylim=c(0,1200), ylab="N_t+1", main=paste("r=",r[i]))
lines(N_t,N_t, type="l", col="Red")}

## 4e. Try r = 2.1. Does the population reach carrying capacity? Will it ever? 
## How about when r = 2.5? Does this look like any natural population you’ve 
## heard about?
r<-c(2.1,2.5)

## Abundance vs. time
N<-matrix(NA, nrow=length(Year), ncol=(length(r)))
N[1,]<-1
for (i in 1:length(r)){
for (t in 1:(max(Year)-1)){
N[t+1,i]<-N[t,i]+r[i]*N[t,i]*(1-N[t,i]/K)}}

## N_t+1 vs. N_t
N_t1<-matrix(NA, nrow=length(N_t), ncol=(length(r)))
for (i in 1:length(r)){
for (t in 1:length(N_t)){
N_t1[t,i]<-N_t[t]+r[i]*N_t[t]*(1-N_t[t]/K)}}

par(mfrow=c(2,2))
for (i in 1:length(r)){
plot(N[,i]~Year, type="l", ylim=c(0,1200), ylab="N", main=paste("r=",r[i]))}
for (i in 1:length(r)){
plot(N_t1[,i]~N_t, type ="l", col="Blue", ylim=c(0,1200), ylab="N_t+1", main=paste("r=",r[i]))
lines(N_t,N_t, type="l", col="Red")}

## 4f. Try r = 3.0. Describe what you see. Does it look like something you’d 
## expect from a simple equation with 2 parameters? Is it random?
r<-c(3.0)

## Abundance vs. time
N<-matrix(NA, nrow=length(Year), ncol=(length(r)))
N[1,]<-1
for (i in 1:length(r)){
for (t in 1:(max(Year)-1)){
N[t+1,i]<-N[t,i]+r[i]*N[t,i]*(1-N[t,i]/K)}}

## N_t+1 vs. N_t
N_t1<-matrix(NA, nrow=length(N_t), ncol=(length(r)))
for (i in 1:length(r)){
for (t in 1:length(N_t)){
N_t1[t,i]<-N_t[t]+r[i]*N_t[t]*(1-N_t[t]/K)}}
 
par(mfrow=c(1,2))
for (i in 1:length(r)){
plot(N[,i]~Year, type="l", ylim=c(0,1300), ylab="N", main=paste("r=",r[i]))}
for (i in 1:length(r)){
plot(N_t1[,i]~N_t, type ="l", col="Blue", ylim=c(0,1300), ylab="N_t+1", main=paste("r=",r[i]))
lines(N_t,N_t, type="l", col="Red")}

# Step 5. Create a bifurcation graph (this is complicated – follow directions carefully)

## 5a. Construct a table which gives the abundance of the population for each of the last 
## 10 years for a range of different values of ‘r’. 

## Put in values of ‘r’ starting with 1.80 increasing to 3.00 in steps of 0.05 
r <- seq(1.8,3, by=0.05)

## 5b. Graph the final 10 abundances of the population for each value of r. 
N<-matrix(NA, nrow=length(Year), ncol=(length(r)))
N[1,]<-1
for (i in 1:length(r)){
for (t in 1:(max(Year)-1)){
N[t+1,i]<-N[t,i]+r[i]*N[t,i]*(1-N[t,i]/K)}}

## Remove all but the last 10 years (rows) of abundances 
N<-N[-(1:90),]

## Create a scatterplot graph with the first row (the r values) as the X axis and the rest 
## as 10 different Y values for each X.
plot(N[1,]~r, type="p", pch=1, ylab="N", main="Bifurcation graph")
for (t in 2:length(N[,1])){
points(r,N[t,], type="p", pch=t)}

## How does the bifurcation graph you just created relate to the graph of abundance over time? 
## Do you think the bifurcation graph would change if we’d used the abundances from years 191 
## to 200 instead of 91 to 100?

## If we made a bifurcation graph of a larger table, say of all abundances from year 200-300, 
## what would we expect it to look like for values of ‘r’ between 2.7 and 3.0?

# Homework: with these simple models, you can’t have chaos without a delay in the effect of 
# density-dependence. Rewrite the logistic model so that the abundance of the population reduces 
# productivity two years in the future, not one year. At what values of  ‘r’ do you start seeing 
# cyclic and chaotic behavior? 
r<-c(0.2, 0.6, 1, 1.4)

## Create vector to hold abundance (N) values
N<-matrix(NA, nrow=length(Year), ncol=(length(r)))

## In year 1, the abundance will be 1
N[1,]<-1

## Use exponential growth equation to get abundance in second year
N[2,]<-N[1,]+r*N[1,]

## Use logistic growth equation to calculate remaining abundances
for (i in 1:length(r)){
for (t in 1:(max(Year)-2)){
N[t+2,i]<-N[t+1,i]+r[i]*N[t+1,i]*(1-N[t,i]/K)}}

par(mfrow=c(2,2))
for (i in 1:length(r)){
plot(N[,i]~Year, type ="l", ylim=c(0,1600), ylab="N", main=paste("r=",r[i]))}