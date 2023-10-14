##############################################################################
# Non-systematic surveys reveal increases in areas occupied by endangered and 
# data-deficient Nubian bustard

# Ramiro D. Crego, Haydée Hernández-Yáñez, Thomas Rabeil, Yves Hingrat, Peter Leimgruber, and Jared A. Stabach

## Load libraries
library(dplyr)
library(plyr)
library(ggplot2)
library(sf)

## Load data
load("./Data/DataNubian")
load("./Data/GridIDs")
Grid <- st_read("./Data/Grid.shp")

### Arrange dataset with observations at each secondary visit per season
yearseason <- sort(unique(Track_cellNub$YearSeason))
ncells <- nrow(Grid)

YNub <- array(NA, dim = c(ncells,3,length(yearseason)))
for(i in 1:length(yearseason)){
  temp1 <- filter(Track_cellNub, YearSeason == yearseason[i])
  temp2 <- temp1 %>% group_by(CellID) %>% dplyr::reframe(RepID = seq_along(CellID), Count = Count)
  W <- data.frame(rep(1:ncells, rep(10, ncells)), rep(1:10, ncells))
  colnames(W) <- c("CellID", "RepID")
  J <- full_join(W, temp2, by = c("CellID", "RepID"))
  X <- split(J$Count, f = J$CellID)
  X <- do.call(rbind, X)
  YNub[,,i] <- X[,1:3] #Truncate to 3 secondary sampling
} #y

YNub #This data shows the secondary replicates per season_year to loop over in the model

# Transfrom into presence/absence
YNub[YNub >= 1] <- 1

##########################
# Create IDs and indices


# Site IDs
site <- seq(1,ncells,1)

# Number of sites
nsites <- as.integer(length(site))

# Replicate ID
rep <- seq(1:3)

# Number of reps
nreps <- as.integer(length(rep))

# Season ID
season <- seq(1:length(yearseason))

# Number of seasons
nseasons <- as.integer(length(season))

# Reps sampled per season
repsampled<-list()
for (s in 1:nseasons){
  xA <- YNub[,,s]
  xB <- which(apply(xA, 1, function(y) !all(is.na(y))), arr.ind = T) # Extract row name that was surveyed each saason
  repsampled[[s]] <- xB
}

nXseason <- as.data.frame(sapply(repsampled, length))
colnames(nXseason) <- "Nr cells sampled"
rownames(nXseason) <- yearseason
nXseason$PropSampled <- nXseason$`Nr cells sampled`/ncells

#write.csv(nXseason, "NrCellsSampledPerSeason.csv")


####################################
# Covariates detection probability 

# Transect length
load("./Data/Track_cellTra")

TransLength <- array(NA, dim = c(ncells,3,length(yearseason)))
for(i in 1:length(yearseason)){
  temp1 <- filter(Track_cellTra, YearSeason == yearseason[i])
  temp2 <- temp1 %>% group_by(CellID) %>% dplyr::summarise(RepID = seq_along(CellID), Length = Length)
  W <- data.frame(rep(1:ncells, rep(10, ncells)), rep(1:10, ncells))
  colnames(W) <- c("CellID", "RepID")
  J <- full_join(W, temp2, by = c("CellID", "RepID"))
  X <- split(J$Length, f = J$CellID)
  X <- do.call(rbind, X)
  TransLength[,,i] <- X[,1:3] #Truncate to 3 secondary sampling
} #y

# Complete length of missing replicates with the mean length of transects
for(s in 1:length(yearseason)){
  for(i in 1:ncells){
    TransLength[i,which(is.na(TransLength[i,,s])==TRUE),s] = mean(TransLength[i,,s],na.rm = T)
  }}

MeanTrL <- mean(TransLength, na.rm = T)
# 3.075959
SDTrL <- sd(TransLength, na.rm = T)
# 1.327599

TransLength[is.na(TransLength)] <- MeanTrL

#Standarize transect length
TransLengthSt <- TransLength
for(s in 1:length(yearseason)){
  for(i in 1:ncells){
    for(j in 1:3){
      TransLengthSt[i,j,s] = (TransLength[i,j,s] - MeanTrL)/SDTrL
    }}}

#Data for autocorrelation structure
library(spdep)
centroid <- st_centroid(Grid)

neigh <- dnearneigh(centroid, d1 = 0, d2 = 7000)
bugsnb <- nb2WB(neigh)
bugsnb$num

# Combine data
Modeldata <- list(YNub, nsites, nreps, nseasons, repsampled, yearseason, TransLengthSt, MeanTrL, SDTrL, bugsnb$adj, bugsnb$weights, bugsnb$num)
heads <- c("YNub", "nsites", "nreps", "nseasons", "repsampled", "yearseason", "TransLengthSt", "MeanTrL", "SDTrL", "adj", "weights", "num")

Modeldata <- setNames(Modeldata, nm = heads)
#Save arrays
save(Modeldata, file = "./ModelData/dataCARHydra")



#####################################
## Calculate survey dates per season

TransDates <- data.frame()
for(i in 1:length(yearseason)){
  temp1 <- filter(Track_cellTra, YearSeason == yearseason[i])
  temp2 <- temp1 %>% group_by(CellID) %>% dplyr::summarise(RepID = seq_along(CellID), Date = Date)
  W <- data.frame(rep(1:ncells, rep(10, ncells)), rep(1:10, ncells))
  colnames(W) <- c("CellID", "RepID")
  J <- full_join(W, temp2, by = c("CellID", "RepID"))
  X <- J %>% filter(RepID <= 3)
  I <- data.frame(min = min(X$Date, na.rm = T), max = max(X$Date, na.rm = T), mean = mean(X$Date, na.rm = T))
  TransDates <- rbind(TransDates,I) #Truncate to 3 secondary sampling
} #y


dates <- TransDates
dates$Diff <- dates$max - dates$min
save(dates, file = "./ModelData/datesDry")

#Estimate number of replicates per site per season
n<-list()
for (s in 1:nseasons){
  xA <- YNub[,,s]
  xB <- rowSums(!is.na(xA))
  n[[s]] <- count(xB)
}

replic <- matrix(NA, nrow = 10, ncol = 1)
for (s in 1:nseasons){
  temp <- n[[s]]
  replic[s,]    <- temp[3,2] + temp[4,2]
}
replic <- as.data.frame(replic)
colnames(replic) <- "Nr cell with secondary visits"

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

Summary <- cbind(nXseason,replic,dates,Nubnaiveocc[,1])
#write.csv(Summary, "./Summary.csv")