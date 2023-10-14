##############################################################################
# Non-systematic surveys reveal increases in areas occupied by endangered and 
# data-deficient Nubian bustard

# Ramiro D. Crego, Haydée Hernández-Yáñez, Thomas Rabeil, Yves Hingrat, Peter Leimgruber, and Jared A. Stabach

# Load libraries
library(tidyverse)
library(sf)
library(tmap)
library(stars)

# Load grid
Grid <- st_read("./Data/Grid.shp")

#Reproject to WGS84 lat long
GridWGS84 <- st_transform(Grid, crs = "+proj=longlat +datum=WGS84")
GridWGS84$const <- rep(1:2, each=4000)[1:nrow(GridWGS84)] 
#Extract time for each survey and attach to transects
load("./ModelData/datesDry")
dates

## Extract covariates 
library(rgee)
ee_Initialize()
ee_check()


# Extract Elevation
elev = ee$Image('USGS/SRTMGL1_003')

meanElev <- data.frame()
for(i in 1:2){
  Gridsub <- filter(GridWGS84, const == i)
  temp = ee_extract(x=elev, y=Gridsub, fun=ee$Reducer$mean(), sf=F, scale=30)
  meanElev <- rbind(meanElev, temp)
}#i
hist(meanElev$elevation)

(meanE <- mean(meanElev$elevation))
#Dry 421.9051   #Full grid 415.2194
(sdE <- sd(meanElev$elevation))
#Dry 35.86525   #Full grid 41.96866
meanElevSt <- (meanElev$elevation - meanE)/sdE

save(meanElev, file = "./ModelData/ElevationCAR")
save(meanElevSt, file = "./ModelData/ElevationStCAR")

# Extract Surface Roughness

sr = elev$reduceNeighborhood(reducer = ee$Reducer$stdDev(), kernel = ee$Kernel$circle(radius = 10, units = "pixels"), skipMasked = T)$rename("SR")

meanSR <- data.frame()
for(i in 1:2){
  Gridsub <- filter(GridWGS84, const == i)
  temp = ee_extract(x=sr, y=Gridsub, fun=ee$Reducer$mean(), sf=F, scale=30)
  meanSR <- rbind(meanSR, temp)
}#i

(meansr <- mean(meanSR$SR))
#Dry 6.323302   #Full grid 5.055114
(sdsr <- sd(meanSR$SR))
#Dry 2.441323 #Full grid 2.125241
meanSRSt <- (meanSR$SR - meansr)/sdsr

save(meanSR, file = "./ModelData/SRCAR")
save(meanSRSt, file = "./ModelData/SRStCAR")


# Extract NDVI for each stratum in Google Earth Engine
dates <- as.data.frame(dates)

out <- data.frame()
for (t in 1:nrow(dates)){
  tmin <- factor(dates[t,1] - 8)
  tmax <- factor(dates[t,2] + 8)
  date1 = ee$String(tmin)
  date2 = ee$String(tmax)
  
  # NDVI
  modisMaxNDVI = ee$ImageCollection('MODIS/006/MOD13Q1')$filterDate(date1,date2)$select('NDVI')$max()$multiply(0.0001)
  for(i in 1:2){
    Gridsub <- filter(GridWGS84, const == i)
    modisNDVI = ee_extract(x=modisMaxNDVI, y=Gridsub, fun=ee$Reducer$mean(), sf=F, scale=250)
    modisNDVI$iter <- t
    out <- rbind(out, modisNDVI)
  }#i
} #t


# Arrange data 

MaxNDVI <- matrix(NA, nrow = nrow(GridWGS84), ncol = nrow(dates))
for (t in 1:nrow(dates)){
  temp <- out %>% filter(iter == t)
  MaxNDVI[,t] <-temp[,3]
}

(meanndvi <- mean(MaxNDVI))
#Dry 0.1027951  #Full grid 0.1115856
(sdndvi <- sd(MaxNDVI))
#Dry 0.0298537  #Full grid 0.02248806

MaxNDVIST <- matrix(NA, nrow = nrow(GridWGS84), ncol = nrow(dates))
for (t in 1:nrow(dates)){
  MaxNDVIST[,t] <- (MaxNDVI[,t] - meanndvi)/sdndvi
}
save(MaxNDVI, file = "./ModelData/NDVICAR")
save(MaxNDVIST, file = "./ModelData/NDVIStCAR")

##Check correlation
for(s in 1:10){
  print(cor(data.frame(meanElevSt, meanSRSt, MaxNDVIST[,s])))
}


# Map of all covariates
library(sf)
library(tmap)

load("./ModelData/NDVICAR")
load("./ModelData/ElevationCAR")
load("./ModelData/SRCAR")
Grid <- st_read("./Data/FullGrid.shp")

GridCov <- cbind(Grid,data.frame(Elevation = meanElev$elevation),data.frame(SR = meanSR$SR),MaxNDVI)

M1 <- tm_shape(GridCov) + tm_fill(col = "Elevation", palette = "inferno", legend.show = T, title = "", title.size = 0.6) + tm_layout("Elevation (m)", title.size = 1.2, title.position = c('center','top'), title.color = "white", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

M2 <- tm_shape(GridCov) + tm_fill(col = "SR", palette = "inferno", legend.show = T, title = "", title.size = 0.6) + tm_layout("Surface roughness", title.size = 1.2, title.position = c('center','top'), title.color = "white", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

M3 <- tm_shape(GridCov) + tm_fill(col = "X1", palette = "Greens", legend.show = T, title = "", title.size = 0.6) + tm_layout("NDVI 2008", title.size = 1.2, title.position = c('center','top'), title.color = "black", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

M4 <- tm_shape(GridCov) + tm_fill(col = "X2", palette = "Greens", legend.show = T, title = "", title.size = 0.6) + tm_layout("NDVI 2009", title.size = 1.2, title.position = c('center','top'), title.color = "black", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

M5 <- tm_shape(GridCov) + tm_fill(col = "X3", palette = "Greens", legend.show = T, title = "", title.size = 0.6) + tm_layout("NDVI 2010", title.size = 1.2, title.position = c('center','top'), title.color = "black", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

M6 <- tm_shape(GridCov) + tm_fill(col = "X4", palette = "Greens", legend.show = T, title = "", title.size = 0.6) + tm_layout("NDVI 2011", title.size = 1.2, title.position = c('center','top'), title.color = "black", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

M7 <- tm_shape(GridCov) + tm_fill(col = "X5", palette = "Greens", legend.show = T, title = "", title.size = 0.6) + tm_layout("NDVI 2012", title.size = 1.2, title.position = c('center','top'), title.color = "black", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

M8 <- tm_shape(GridCov) + tm_fill(col = "X6", palette = "Greens", legend.show = T, title = "", title.size = 0.6) + tm_layout("NDVI 2013", title.size = 1.2, title.position = c('center','top'), title.color = "black", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

M9 <- tm_shape(GridCov) + tm_fill(col = "X7", palette = "Greens", legend.show = T, title = "", title.size = 0.6) + tm_layout("NDVI 2014", title.size = 1.2, title.position = c('center','top'), title.color = "black", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

M10 <- tm_shape(GridCov) + tm_fill(col = "X8", palette = "Greens", legend.show = T, title = "", title.size = 0.6) + tm_layout("NDVI 2015", title.size = 1.2, title.position = c('center','top'), title.color = "black", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

M11 <- tm_shape(GridCov) + tm_fill(col = "X9", palette = "Greens", legend.show = T, title = "", title.size = 0.6) + tm_layout("NDVI 2016", title.size = 1.2, title.position = c('center','top'), title.color = "black", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

M12 <- tm_shape(GridCov) + tm_fill(col = "X10", palette = "Greens", legend.show = T, title = "", title.size = 0.6) + tm_layout("NDVI 2017", title.size = 1.2, title.position = c('center','top'), title.color = "black", legend.title.size=1, legend.text.size = 0.4, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

mapCov <- tmap_arrange(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, M11, M12, ncol = 4)

tmap_save(mapCov, "./Figures/CovariablesMap.jpg", width=8, height=6, units = "in")