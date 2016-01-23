#############################################
#University of Alaska Fairbanks
#FISH 421: Fisheries Population Dynamics
#Lab 5: Yield per recruit
#Instructor/lab author: Milo Adkison
#TA/code author: Phil Ganz
#############################################

# Step 1. Create a model of the population and the fishery

## Parameters
F <- 0           ## Fishing mortality
t_c <- 2         ## Age at first capture
M <- 0.2         ## Natural mortality
k <- 0.5         ## Brody growth coefficient
t_0 <- -0.2      ## Hypothetical age at which weight is 0
W_inf <- 10      ## Asymptotic weight
b <- 3           ## Dimension exponent to get from length to volume

## Age
Age <- c(1:9)

## Abundance
N <- numeric(9)
N[1] <- 1000
for (i in 2:length(Age)){
	if ((i-1)>=t_c){
		N[i] <- N[i-1]*exp(-M-F)}
	else {
		N[i] <- N[i-1]*exp(-M)}}

## Weight of an individual
W <- numeric(9)
for (i in 1:length(Age)){
	W[i] <- W_inf*(1-exp(-k*(i-t_0)))^b}

## Total biomass of cohort
B <- N*W

## Catch in numbers (Baronov catch equation)
C <- numeric(9)
for (i in 1:length(Age)){
	if (i>=t_c){
		C[i] <- N[i]*F/(F+M)*(1-exp(-F-M))}
	else {
		C[i] <- 0}}

## Yield in biomass
Y <- C*W

## Plot cohort biomass by age
plot(Age, B, type="l")

# Q: What’s the critical age? What does ‘critical age’ mean?
# Q: Why did you set F=0 to calculate the critical age?
# Q: If your spreadsheet population were oysters on your fish 
# farm and you were paid on a per pound basis, how would these 
# calculations help you?
# Q: What if bigger oysters were worth more per pound than 
# small ones? What would you change in these calculations?

# 2. Examine how the fishery control parameters 
# F (fishing mortality rate) and tc (first age 
# caught) affect the yield per recruit.

## Make a table of Y/R calculations for different 
## values of F and tc.

F <- seq(0,1, by=0.1)
t_c <- seq(1:9)

## Abundance
N <- array(dim=c(length(Age), length(F), length(t_c)))
N[1,,] <- 1000
for (i in 2:length(Age)){                              # Creates a row for each age
	for (j in 1:length(F)){                            # Creates a column for each level of fishing mortality
		for (k in 1:length(t_c)){                      # Creates a separate matrix for each value of T_c
	if ((i-1)>=k){
		N[i,j,k] <- N[i-1,j,k]*exp(-M-F[j])}
	else {
		N[i,j,k] <- N[i-1,j,k]*exp(-M)}}}}

## Total biomass of cohort
B <- N*W                                               

## Catch in numbers (Baronov catch equation)
C <- array(dim=c(length(Age), length(F), length(t_c)))
for (i in 1:length(Age)){
	for (j in 1:length(F)){
		for (k in 1:length(t_c)){
	if (i>=k){
		C[i,j,k] <- N[i,j,k]*F[j]/(F[j]+M)*(1-exp(-F[j]-M))}
	else {
		C[i,j,k] <- 0}}}}

## Yield in biomass
Y <- C*W

## Yield per recruit
YR <- matrix(nrow=length(t_c), ncol=length(F))
for (j in 1:length(F)){
	for (k in 1:length(t_c)){
		YR[k,j] <- sum(Y[,j,k])/1000}}

colnames(YR) <- F                                      # Label columns of YR matrix with F values
rownames(YR) <- t_c                                    # Label rows of YR matrix with t_c values

## Create heatmap of yield per recruit values for 
## different combinations of F and t_c.
## Run install.packages("reshape2") and/or 
## install.packages("ggplot2") if you dont have
## the necessary packages.
library(reshape2)
melted_YR <- melt(YR)                                  # Reshape matrix into data frame (for plotting)

colnames(melted_YR) <- c("t_c","F","YR")               # Name columns 

library(ggplot2)
ggplot(data = melted_YR, aes(factor(F), factor(t_c), fill = YR))+                       
geom_tile(aes(fill = YR)) + 
geom_text(aes(fill = melted_YR$YR, label = round(melted_YR$YR, 2)))+
scale_fill_gradient(low = "white", high = "green")+
ggtitle("Yield per recruit")+ 
xlab("F")+
ylab("t_c")+
ylim(rev(levels(factor(melted_YR$t_c))))                # Reverse the default order of yaxis labels

# Q: What F and tc give the maximum Y/R? Are they on 
# this table?
# Q: If F was 0.8 and you could not control effort, 
# what mesh size limits would you use?
# Q: If you could not effectively prevent fish of 2 years 
# old or older from being caught, how hard should you fish? 
# *** Write this down - this is Fmax ***
# Q: From an enforcement standpoint, is it easier to limit 
# gear or effort?

# 3. Determine how to get the maximum yield from a fishery 
# without reducing spawning biomass per recruit below 35% 
# of its unfished level. 

## Assume fish less than 5 yrs old do not reproduce, 1/2 
## of 6-year olds are mature, and all older fish are mature. 
mat <- c(0,0,0,0,0,0.5,1,1,1)

## Calculate SSB/R for an unfished stock
SSB_0 <- sum(B[,1,1]*mat/N[1,1,1])

## Write a function that calculates the percentage of spawning 
## biomass per recruit (when compared to unfished SSB/R) for 
## any values of F and t_c.
GetSSB <- function(F,t_c){
## Abundance
N <- numeric(9)
N[1] <- 1000
for (i in 2:length(Age)){
	if ((i-1)>=t_c){
		N[i] <- N[i-1]*exp(-M-F)}
	else {
		N[i] <- N[i-1]*exp(-M)}}

## Total biomass of cohort
B <- N*W

## SSB
SSB <- sum(B*mat)

## SSB per recruit
SSB_R <- SSB/N[1]

## Percentage of SSB unfished
SSB_P <- SSB_R/SSB_0
return(SSB_P)}

## Try the function out!
GetSSB(0.4,2)
GetSSB(0.2,4)

## Assume that tc = 2. 
t_c <- 2

## How hard can you fish without driving SSB/R below 35% of 
## unfished level?
## To do this, we create another function of F that is SSB
## minus our target (in this case, 35%). This way, when we minimize
## (set=0) the difference function by adjusting F, our F will be at 
## a level that results in 35% SSB.
d <- 0.35
diff <-function(F){
	GetSSB(F,t_c)-d}

uniroot(diff,interval=c(0,1))

## Create a new table of the yield per recruit values for which SSB 
## remains above 35% of the unfished level.

# Apply GetSSB function to each combination of F and t_c in the data frame
melted_YR$SSB <- apply(melted_YR, 1, function(x){GetSSB(x[2],x[1])})  

# Mark as NA the yield per recruit values for which the associated SSB is < 35%   
is.na(melted_YR$YR) <- melted_YR$SSB<0.35

# Use the exact same plot code as before for plot with YRs omitted when SSB < 35%
ggplot(data = melted_YR, aes(factor(F), factor(t_c), fill = YR))+                       
geom_tile(aes(fill = YR)) + 
geom_text(aes(fill = melted_YR$YR, label = round(melted_YR$YR, 2)))+
scale_fill_gradient(low = "white", high = "green")+
ggtitle("Yield per recruit")+ 
xlab("F")+
ylab("t_c")+
ylim(rev(levels(factor(melted_YR$t_c))))