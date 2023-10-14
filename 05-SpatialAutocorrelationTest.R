##############################################################################
# Non-systematic surveys reveal increases in areas occupied by endangered and 
# data-deficient Nubian bustard

# Ramiro D. Crego, Haydée Hernández-Yáñez, Thomas Rabeil, Yves Hingrat, Peter Leimgruber, and Jared A. Stabach

# Spatial autocorrelation test

# Load data
library(sf)
library(dplyr)
library(spdep)

load("./ModelData/dataCARHydra")
attach(Modeldata)

Grid <- st_read("./Data/Grid.shp")
centroid <- st_centroid(Grid)
neigh <- dnearneigh(centroid, d1 = 0, d2 = 7000)
lw <- nb2listw(neigh,style="B")

load("./ModelOutput/mcCARNub")
yres <- grepl(colnames(mcparams), pattern = "yres\\[")
yres_Post <- mcparams[,yres]
yresMean <-data.frame(yres = apply(yres_Post, MARGIN = 2, mean)) 

yresMean$Site <- rep(1:nsites, length(yresMean))
yresMean$Season <- rep(1:nseasons, each=nsites)

library(pgirmess)
list <- list()
for(s in 1:10){
  temp1 <- filter(yresMean, Season == 1)
  temp1 <- as.vector(temp1$yres)
  corrpgir <- correlog(st_coordinates(centroid)[repsampled[[s]],], temp1[repsampled[[s]]], method = "Moran")
  plot(corrpgir)
  list[[s]] <- corrpgir
}
list
