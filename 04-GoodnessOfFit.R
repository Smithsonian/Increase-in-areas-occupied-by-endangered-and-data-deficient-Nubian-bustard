##############################################################################
# Non-systematic surveys reveal increases in areas occupied by endangered and 
# data-deficient Nubian bustard

# Ramiro D. Crego, Haydée Hernández-Yáñez, Thomas Rabeil, Yves Hingrat, Peter Leimgruber, and Jared A. Stabach

#### Model Checks and Goodness-of-fit
library(MCMCvis)
library(mcmcOutput)
library(dplyr)

# Load Model Data
load("./ModelData/dataCARHydra")
attach(Modeldata)

# Nubian
load("./ModelOutput/mcCARNub")

# Visualy inspect chains
diagPlot(mcparams,c("mu_a0", "sig_a0", "alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "mu_b0", "beta0", "beta1","occtrend"))

# Load nubian mc
load("./ModelOutput/mcCARSummaryNub")
Molsummary



# Goodness-of-fit
## Estimate number of observations in each site based on bernouli trials limited to the number of j replicates per site

# Get sites sampled per season
Y = YNub
repsampled<-list()
for (s in 1:nseasons){
  xA <- Y[,,s]
  xB <- which(apply(xA, 1, function(y) !all(is.na(y))), arr.ind = T) # Extract row name that was surveyed each saason
  repsampled[[s]] <- xB
}

#Estimate number of replicates per site per season
n<-list()
for (s in 1:nseasons){
  xA <- Y[,,s]
  xB <- rowSums(!is.na(xA))
  n[[s]] <- xB
}

# Calculate number of observations per site per season
y<-list()
for (s in 1:nseasons){
  xA <- Y[,,s]
  xB <- rowSums(xA, na.rm=TRUE)
  y[[s]] <- xB
}

## Arrange p
prows <- grepl(colnames(mcparams), pattern = "p\\[")

pss <- mcparams[,prows]

ps <- array(NA, dim = c(nrow(pss), nsites, nreps, nseasons))

seas <-  seq(0,(nsites*nreps*nseasons), (nsites*nreps))
num <- seq(0,(nsites*nreps),nsites)
for(s in 1:nseasons){
  temp <- pss[, (1+seas[s]):(seas[s+1])]
  for(i in 1:nreps){
    ps[,,i,s] <- temp[,(1+num[i]):(num[i+1])]
  }}

psSum <- array(NA, dim = c(nrow(pss), nsites, nseasons))

for(s in 1:nseasons){
  for(i in repsampled[[s]]) {
    for(j in n[[s]][i]){
      psSum[,i,s] <- ifelse(j ==1, ps[,i,1:j,s], rowSums(ps[,i,1:j,s]))
    }}}

# Arrange z
Zrows <- grepl(colnames(mcparams), pattern = "z\\[")
zss <- mcparams[,Zrows]

zs <- array(NA, dim = c(nrow(pss), nsites, nseasons))

numz <- seq(0,(nsites*nseasons), nsites)
for(s in 1:nseasons){
  zs[,,s] <- zss[,(1+numz[s]):(numz[s+1])]
}


# Calculate Freeman-Tukey Bayesian p-value
Tobs2 <- Tsim2 <- array(NA, dim = c(nrow(zs), nseasons))
for(s in 1:nseasons){
  for(iter in 1:nrow(zs)) {
    yTemp <- array(NA, dim = c(nsites, nreps)) # Create empty array to fill with a simulated dataset
    for(i in repsampled[[s]]){
      for(j in 1:n[[s]][i]){
        yTemp[i,j] <- rbinom(1,1,(ps[iter,i,j,s]*zs[iter,i,s])) # Simulate dataset
      }}
    Tobs2[iter,s] <- sum((sqrt(y[[s]][repsampled[[s]]]) - sqrt(psSum[iter,repsampled[[s]],s] * zs[iter,repsampled[[s]],s]))^2) # Freeman-Tukey statistic for observed data
    ySim <- rowSums(yTemp[repsampled[[s]],], na.rm=TRUE) # Calculate total number of observations from simulated data
    Tsim2[iter,s] <- sum((sqrt(ySim) - sqrt(psSum[iter,repsampled[[s]],s] * zs[iter,repsampled[[s]],s]))^2) # Freeman-Tukey statistic for simulated data
  }
}

Tobs22 <- apply(Tobs2, 1, sum)
Tsim22 <- apply(Tsim2, 1, sum)

MASS::eqscplot(Tobs22, Tsim22, xlim=range(Tobs22, Tsim22), ylim=range(Tobs22, Tsim22),
               xlab="Observed data", ylab="Simulated data")
abline(0, 1, lwd=2, col='red')
mean(Tsim22 > Tobs22) # the P value


### Compare naive occupancy
SimNaivOcc <- data.frame()
for(iter in 1:nrow(zs)){
  temp2 <- data.frame()
  for(s in 1:nseasons){
    
    yTemp <- array(NA, dim = c(nsites, nreps))
    for(i in repsampled[[s]]){
      for(j in 1:n[[s]][i]){
        yTemp[i,j] <- rbinom(1,1,(ps[iter,i,j,s]*zs[iter,i,s]))
      }}
    Sum1<-apply(yTemp[repsampled[[s]],], 1, sum,na.rm=TRUE)
    sum2<-replace(Sum1,which(Sum1>0),1)
    sum3<-sum(sum2)
    temp<- data.frame(y = sum3/length(repsampled[[s]]), Season = yearseason[s], Sim = iter)
    temp2 <- rbind(temp2, temp)
  }
  SimNaivOcc <- rbind(SimNaivOcc, temp2)
}

resultssim <- SimNaivOcc %>% group_by(Season, Sim) %>% summarise(mean(y))

finalsim <- resultssim %>% group_by(Season) %>% summarise(mean(`mean(y)`), sd(`mean(y)`)) %>% as.data.frame()
colnames(finalsim) <- c("Season","Mean","SD")

#Nubian Naive occupancy
Nubnaiveocc <- data.frame()
for (t in 1:nseasons){
  x<-YNub[repsampled[[t]],,t]
  (sum1<-apply(x, 1, sum,na.rm=TRUE))
  sum2<-replace(sum1,which(sum1>0),1)
  sum3<-sum(sum2)
  temp<- data.frame(y = sum3/length(repsampled[[t]]), Season = yearseason[t])
  Nubnaiveocc <- rbind(Nubnaiveocc, temp)
}

table <- cbind(Nubnaiveocc,finalsim)
table$dif <- table$y - table$Mean
table
write.csv(table, "SimvsObsNubianData.csv")


#### Check for prior influence on posterior estimates

#simulate data from the prior used in the model
#number of iterations should equal the number of draws times the number of chains 
PR <- rnorm(nrow(mcparams), 0, 10) 

MCMCtrace(mcparams, params = 'mu_a0', priors = PR, pdf = FALSE)
MCMCtrace(mcparams, params = 'alpha0', priors = PR, pdf = FALSE)
MCMCtrace(mcparams, params = 'alpha1', priors = PR, pdf = FALSE)
MCMCtrace(mcparams, params = 'alpha2', priors = PR, pdf = FALSE)
MCMCtrace(mcparams, params = 'alpha3', priors = PR, pdf = FALSE)
MCMCtrace(mcparams, params = 'alpha4', priors = PR, pdf = FALSE)
MCMCtrace(mcparams, params = 'mu_b0', priors = PR, pdf = FALSE)
MCMCtrace(mcparams, params = 'beta0', priors = PR, pdf = FALSE)
MCMCtrace(mcparams, params = 'beta1', priors = PR, pdf = FALSE)
MCMCtrace(mcparams, params = 'occtrend', priors = PR, pdf = FALSE)
MCMCtrace(mcparams, params = 'sig_a0', priors = PR, pdf = FALSE)
MCMCtrace(mcparams, params = 'sig_b0', priors = PR, pdf = FALSE)

# Following Gimenez et al. (2009), overlap should be < 0.35