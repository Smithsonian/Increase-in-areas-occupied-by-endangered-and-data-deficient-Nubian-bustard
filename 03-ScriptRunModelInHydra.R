##############################################################################
# Non-systematic surveys reveal increases in areas occupied by endangered and 
# data-deficient Nubian bustard

# Ramiro D. Crego, Haydée Hernández-Yáñez, Thomas Rabeil, Yves Hingrat, Peter Leimgruber, and Jared A. Stabach


# Script for Hydra (note working directory was specified for the Hydra cluster)
setwd("/scratch/genomics/cregord/bustardproj")
load("./dataCARHydra")
attach(Modeldata)
load("./covariatesCARHydra")
attach(Covariates)

# Load libraries

library(nimble)
library(coda)
library(mcmcOutput)

# Specify model in BUGS language
occmod <- nimbleCode({
  
  # Priors
  #Alpha0
  mu_a0 ~ dnorm(0, 0.01)        #Mean
  sig_a0 ~ dunif(0, 10)      #SD
  tau_a0 <- pow(sig_a0, -2)    #Precision
  
  for (t in 1:nseasons) {
    alpha0[t] ~ dnorm(mu_a0, tau_a0)
  }
  
  #Beta0
  mu_b0 ~ dnorm(0, 0.01)        #Mean
  sig_b0 ~ dunif(0, 10)      #SD
  tau_b0 <- pow(sig_b0, -2)    #Precision
  
  for (t in 1:nseasons) {
    beta0[t] ~ dnorm(mu_b0, tau_b0)
  }
  
  #Alpha1
  alpha1 ~ dnorm(0, 0.01) 
  
  #Alpha2
  alpha2 ~ dnorm(0, 0.01)
  
  #Alpha3
  alpha3 ~ dnorm(0, 0.01)
  
  #Alpha4
  alpha4 ~ dnorm(0, 0.01)
  
  #Beta1
  beta1 ~ dnorm(0, 0.01)
  
  # CAR prior distribution for spatial random effects
  eta[1:nsites] ~ dcar_normal(adj[], weights[], num[], tau, zero_mean=1)
  
  # Parameter expansion for spatial random effects for better mixing
  for (i in 1:nsites){
    S[i] <- Xi * eta[i]
  }
  
  # Hierarchical half-Cauchy prior for random effect standard deviations
  # from Gelman (2006), also see Chelgren et al. (Ecology 2011a)
  Xi ~ dnorm(0, tau.xi)
  tau.xi <- pow(sig.xi, -2)
  sig.xi <- u # Cauchy scale parameter u estimated from data
  u ~ dunif(0, 3) # Prior for the Cauchy scale parameter
  sd.eta <- abs(Xi)/sqrt(tau)
  tau ~ dgamma(0.5, 0.5) # This is a chisquare rv with 1 df
  v.eta <- pow(sd.eta, 2)
  
  # Likelihood
  for(t in 1:nseasons){  
    for (i in 1:nsites) { #start initial loop over the R sites
      # True state model for the partially observed true state
      z[i,t] ~ dbern(psi[i,t])		# True occupancy z at site i
      logit(psi[i,t]) <- lpsi.lim[i,t]
      lpsi.lim[i,t] <- min(1000, max(-1000, lpsi[i,t]))
      lpsi[i,t] <- alpha0[t] + alpha1 * NDVI[i,t] + alpha2 * pow(NDVI[i,t],2) + alpha3 * elev[i] + alpha4 * SR[i] + S[i]
      
      for (j in 1:nreps) { # start a second loop over the T replicates
        # Observation model for the actual observations
        y[i,j,t] ~ dbern(p.eff[i,j,t])	# Detection-nondetection at i and j
        p.eff[i,j,t] <- z[i,t] * p[i,j,t]
        logit(p[i,j,t]) <- lp.lim[i,j,t]
        lp.lim[i,j,t] <- min(1000, max(-1000, lp[i,j,t]))  # 'Stabilize' logit
        lp[i,j,t] <- beta0[t] + beta1 * tranlength[i,j,t]
        
      }#j detection
      
      # Calculate psi residuals
      yres[i,t]<-psi[i,t]-z[i,t]
      
    }# i reps
  }#t seasons
  
})

# Create some constants, data, and initial values to pass to the model builder
Consts <- list(nsites = nsites, nreps = nreps, nseasons=nseasons,  adj = adj, weights = weights, num = num)

Data <- list(y = YNub, NDVI = NDVI, elev = Elev, SR = SR, tranlength = TransLengthSt)

zst <- matrix(1,nrow = nsites, ncol = nseasons)
Inits <- list(
  z = zst, 
  mu_a0 = rnorm(1), 
  mu_b0 = rnorm(1), 
  sig_a0 = runif(1, 0 ,1),
  sig_b0 = runif(1, 0 ,1),
  alpha1 = rnorm(1),
  alpha2 = rnorm(1),
  alpha3 = rnorm(1),
  alpha4 = rnorm(1),
  beta1 = rnorm(1),
  eta = rep(0, nsites)
)

# Parameters to estimate
params <- c("mu_a0", "sig_a0", "alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "mu_b0", "sig_b0", "beta0", "beta1", "yres", "z", "p")


# Use 'foreach' to do a long run in parallel
# ''''''''''''''''''''''''''''''''''''''''''
library(foreach)
library(parallel)
library(doParallel)


ncore <- Sys.getenv("NSLOTS")

cl <- makeCluster(ncore, type = "FORK")
registerDoParallel(cl)

seeds <- 1:ncore

system.time(
  out <- foreach(x = seeds, .packages="nimble",
                  .errorhandling='remove', .inorder=FALSE) %dopar% {
                    set.seed(x)
                    nimbleMCMC(occmod, constants = Consts, data = Data, inits = Inits, monitors=params, nburnin = 100000, niter = 200000, thin = 200, nchains=1, summary = F, samplesAsCodaMCMC=TRUE)
                  } )
stopCluster(cl)

mclist <- coda::mcmc.list(out)
mcparams <- mcmcOutput(mclist)

Molsummary <- summary(mcparams, n.eff=TRUE)

save(mcparams, file = "./mcCARNub")
save(Molsummary, file = "./mcCARSummaryNub")


