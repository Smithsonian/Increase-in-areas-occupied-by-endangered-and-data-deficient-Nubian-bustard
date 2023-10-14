##############################################################################
# Non-systematic surveys reveal increases in areas occupied by endangered and 
# data-deficient Nubian bustard

# Ramiro D. Crego, Haydée Hernández-Yáñez, Thomas Rabeil, Yves Hingrat, Peter Leimgruber, and Jared A. Stabach

## Load libraries
library(bayestestR)
library(dplyr)
library(ggplot2)
library(mcmcOutput)
library(coda)

# Set ggplot theme
theme_set(theme_bw(base_size = 16))

# Load Model Data
load("./ModelData/dataCARHydra")
attach(Modeldata)

load("./ModelData/covariatesCARHydra")

load("./ModelData/datesDry")

Dates <- dates$mean
Years <- format(dates$mean, "%Y")
seasons <- yearseason

# Occupancy Model Nubian

#Nubian Naive occupancy
Nubnaiveocc <- data.frame()
for (t in 1:nseasons){
  x<-YNub[repsampled[[t]],,t]
  sum1<-apply(x, 1, sum,na.rm=TRUE)
  sum2<-replace(sum1,which(sum1>0),1)
  sum3<-sum(sum2)
  temp<- data.frame(y = sum3/length(repsampled[[t]]), Season = yearseason[t])
  Nubnaiveocc <- rbind(Nubnaiveocc, temp)
}


# Load nubian mc
load("./ModelOutput/mcCARNub")

#Params rows 
mualpha0_rowsP <- grepl(colnames(mcparams), pattern = "^mu_a0")
alpha0_rowsP <- grepl(colnames(mcparams), pattern = "^alpha0\\[")
alpha1_rowsP <- grepl(colnames(mcparams), pattern = "^alpha1")
alpha2_rowsP <- grepl(colnames(mcparams), pattern = "^alpha2")
alpha3_rowsP <- grepl(colnames(mcparams), pattern = "^alpha3")
alpha4_rowsP <- grepl(colnames(mcparams), pattern = "^alpha4")
trend_rowsP <- grepl(colnames(mcparams), pattern = "^occtrend")
mubeta0_rowsP <- grepl(colnames(mcparams), pattern = "^mu_b0")
beta0_rowsP <- grepl(colnames(mcparams), pattern = "^beta0\\[")
beta1_rowsP <- grepl(colnames(mcparams), pattern = "^beta1")

SigXi_rowsP <- grepl(colnames(mcparams), pattern = "^sig.xi")
mean(mcparams[,SigXi_rowsP])
SDETA_rowsP <- grepl(colnames(mcparams), pattern = "^sd.eta")
mean(mcparams[,SDETA_rowsP])
VETA_rowsP <- grepl(colnames(mcparams), pattern = "^v.eta")
mean(mcparams[,VETA_rowsP])

# NDVI
mean(mcparams[,alpha1_rowsP]<0)
median(mcparams[,alpha1_rowsP])
ci(mcparams[,alpha1_rowsP], 0.89)

#NDVI2
mean(mcparams[,alpha2_rowsP]<0)
median(mcparams[,alpha2_rowsP])
ci(mcparams[,alpha2_rowsP], 0.89)

#Elev
mean(mcparams[,alpha3_rowsP]>0)
median(mcparams[,alpha3_rowsP])
ci(mcparams[,alpha3_rowsP], 0.89)

#SR
mean(mcparams[,alpha4_rowsP]<0)
median(mcparams[,alpha4_rowsP])
ci(mcparams[,alpha4_rowsP], 0.89)

# Trans length
mean(mcparams[,beta1_rowsP]>0)
median(mcparams[,beta1_rowsP])
ci(mcparams[,beta1_rowsP], 0.89)

# Occ trend
mean(mcparams[,trend_rowsP]>0)
median(mcparams[,trend_rowsP])
ci(mcparams[,trend_rowsP], 0.89)

#########################################
#-Figure caterpiller plot for covariates

params <- cbind(mcparams[,alpha1_rowsP], mcparams[,alpha2_rowsP], mcparams[,alpha3_rowsP], mcparams[,alpha4_rowsP], mcparams[,trend_rowsP], mcparams[,beta1_rowsP])

covnames <- c('NDVI','NDVI^2','Elevation','Surface roughness', "Occupanyc trend",'Transect length' )
paramsdf <- data.frame()
for(i in 1:length(covnames)){
  PSpos <- mean(params[,i] > 0)
  PSneg <- mean(params[,i] < 0)
  temp <- data.frame(Mean = median(params[,i]), ci(params[,i], 0.89),  PS = ifelse(mean(params[,i]) > 0, PSpos, PSneg), Cov = covnames[i])
  paramsdf <- rbind(paramsdf,temp)
}

# visualise beta means & 89% CIs
(FigParams <- ggplot(paramsdf, aes(x=factor(covnames, levels =covnames), y=Mean, colour=PS*100)) +
    geom_hline(yintercept=0,linetype = 'dashed', linewidth=1) + coord_flip() +
    geom_point(size=3) + geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0, size=2) +
    scale_colour_gradient2(name='Posterior\nprobability (%)', low='#FDE725FF', mid = '#22A884FF', high='#440154FF', midpoint = 80) + ylab('Parameter estimate') + xlab(NULL) +
    theme(axis.line=element_line(size=1),
          axis.ticks=element_line(size=1), axis.ticks.length=unit(2, "mm"),
          strip.text=element_text(size=14,face='bold'), strip.background=element_blank(),
          panel.border=element_blank(), panel.grid=element_blank()))

ggsave(file = "./Figures/Fig2.jpg", plot = FigParams)

#########################################################
# Detection probability as a function of transect length 

tl.sim <- seq(0.5, 5, 0.1)
tl.simscal <- (tl.sim - 3.330228)/1.137136

tlPost <- matrix(NA, nrow = nrow(mcparams), ncol = length(tl.simscal))
for (i in 1:nrow(tlPost)){
  tlPost[i,] <- plogis(mcparams[i,mubeta0_rowsP] + mcparams[i,beta1_rowsP] * tl.simscal)
}

tl4plot <- data.frame()
for (i in 1:ncol(tlPost)){
  temp <- data.frame(x = tl.sim[i], y = median(tlPost[,i]), yUpper = ci(tlPost[,i], 0.89)$CI_high, yLower = ci(tlPost[,i], 0.89)$CI_low)
  tl4plot <- rbind(tl4plot, temp)
}

#Plot 

(tl <- ggplot(tl4plot, aes(x = x, y = y)) + geom_ribbon(aes(ymin = yLower, ymax = yUpper), fill = "grey90") + geom_line(size=1, colour = 'black') + labs(y=expression("Detection probability"), x = "Transect length (km)") +
    theme(axis.line=element_line(size=1),
          axis.ticks=element_line(size=1), axis.ticks.length=unit(2, "mm"),
          strip.text=element_text(size=11,face='bold'), strip.background=element_blank(),
          panel.border=element_blank(), panel.grid=element_blank()))

ggsave(file = "./Figures/DetectionProbTLDryCAR.jpg", plot = tl, width = 8, height = 6)

################################
# Detection probability by year

beta0P <- mcparams[,beta0_rowsP]

pPost <- matrix(NA, nrow = nrow(mcparams), ncol = ncol(beta0P))

for(t in 1:nseasons){
  for (i in 1:nrow(pPost)){
    pPost[i,t] <- plogis(beta0P[i,t])
  }}

p4plot <- data.frame()

for (i in 1:nseasons){
  temp <- data.frame(occ = 0, y = median(pPost[,i]), yUpper = ci(pPost[,i], 0.89)$CI_high, yLower = ci(pPost[,i], 0.89)$CI_low, Season = Years[i])
  p4plot <- rbind(p4plot, temp)
}

(FigureDetProb <- ggplot(p4plot, aes(x=Season, y=y)) + 
    labs(y="Detection probability (89% CI)", x = "") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1),
          strip.background=element_blank(), 
          panel.grid=element_blank()) + geom_errorbar(aes(ymin=yLower, ymax=yUpper), width=0, size=1, color = 'darkgrey') + geom_point(size=2) +
    theme(axis.line=element_line(size=1),
          axis.ticks=element_line(size=1), axis.ticks.length=unit(2, "mm"),
          strip.text=element_text(size=14,face='bold'), strip.background=element_blank(),
          panel.border=element_blank(), panel.grid=element_blank())) 


ggsave(file = "./Figures/Fig1.jpg", plot = FigureDetProb, width = 10, height = 6)


################################
# Occupancy probability by year

#Plot 
z_rowsP <- grepl(colnames(mcparams), pattern = "^z\\[")
z <- mcparams[,z_rowsP]

out <- matrix(NA, nrow = nrow(mcparams), ncol = nseasons)

for(i in 1:nrow(mcparams)){
  temp <- data.frame(x = z[i,], season = rep(1:nseasons, each = nsites))
  for(t in 1:nseasons){
    temp1 <- temp %>% filter(season == t)
    out[i,t] <- sum(temp1[,1])/nsites
  }
}

occ4plot <- data.frame()
for (i in 1:nseasons){
  temp <- data.frame(y = median(out[,i]), yUpper = ci(out[,i], 0.89)$CI_high, yLower = ci(out[,i], 0.89)$CI_low, Season = Years[i])
  occ4plot <- rbind(occ4plot, temp)
}

occ4plot$Years<- as.integer(Years)
Nubnaiveocc$Years <- as.integer(Years)

# Plot Occ figure

(Figure3 <- ggplot(occ4plot, aes(x=Years, y=y)) + 
    labs(y="Proportion of areas occupied (89% CI)", x = "") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1),
          strip.background=element_blank(), 
          panel.grid=element_blank()) + geom_errorbar(aes(ymin=yLower, ymax=yUpper), width=0, size=1, color = 'darkgrey') + geom_point(size=2) + geom_point(data = Nubnaiveocc, aes(x= Years, y=y), shape = 17, size=2) +
    theme(axis.line=element_line(size=1),
          axis.ticks=element_line(size=1), axis.ticks.length=unit(2, "mm"),
          strip.text=element_text(size=11,face='bold'), strip.background=element_blank(), panel.border=element_blank(), panel.grid=element_blank()) + scale_x_discrete(limits = as.integer(Years)))

ggsave(file = "./Figures/Fig3.jpg", plot = Figure3, width = 10, height = 6)
