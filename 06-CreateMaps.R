##############################################################################
# Non-systematic surveys reveal increases in areas occupied by endangered and 
# data-deficient Nubian bustard

# Ramiro D. Crego, Haydée Hernández-Yáñez, Thomas Rabeil, Yves Hingrat, Peter Leimgruber, and Jared A. Stabach

# Load libraries
library(dplyr)
library(sf)
library(tmap)
library(ceramic)

# Occupancy map (mean and sd)

# Load nubian mc
load("./ModelData/datesDry")
seasons <- yearseason

load("./ModelOutput/mcCARNub")
z_rowsP <- grepl(colnames(mcparams), pattern = "^z\\[")
z <- mcparams[,z_rowsP]
meanz <- apply(z, 2, mean)
sdz <- apply(z, 2, sd)
meanz <- data.frame(occ = meanz, season = rep(1:nseasons, each = nsites))
sdz <- data.frame(occ = sdz, season = rep(1:nseasons, each = nsites))

out <- matrix(NA, nrow = nsites, ncol = nseasons)
for(t in 1:nseasons){ 
  temp1 <- meanz %>% filter(season == t)
  out[,t] <-temp1[,1]
}
out <- as.data.frame(out)
colnames(out) <- seasons

outsd <- matrix(NA, nrow = nsites, ncol = nseasons)
for(t in 1:nseasons){ 
  temp1 <- sdz %>% filter(season == t)
  outsd[,t] <-temp1[,1]
}
outsd <- as.data.frame(outsd)
colnames(outsd) <- seasons

Grid <- st_read("./Data/Grid.shp")
GridZ <- cbind(Grid, out)
GridSD <- cbind(Grid, outsd)
BustData <- st_read("./Data/BustData.shp")
Nubian <- BustData %>% filter(Species == "Nubian")
Tracks <- st_read("./Data/Tracks.shp")

#Display Tracks and Species data
m1 <- tm_shape(GridZ) + tm_fill(col = "X2008_Dry", palette = "viridis", legend.show = F) + tm_layout("2008", title.size = 1.2, title.position = c('center','top'), title.color = "white") + tm_shape(filter(Nubian, Year == 2008)) + tm_dots(col = "red", size = 0.1)
m2 <- tm_shape(GridZ) + tm_fill(col = "X2009_Dry", palette = "viridis", legend.show = F) + tm_layout("2009", title.size = 1.2, title.position = c('center','top'), title.color = "white")  + tm_shape(filter(Nubian, Year == 2009)) + tm_dots(col = "red", size = 0.1)
m3 <- tm_shape(GridZ) + tm_fill(col = "X2010_Dry", palette = "viridis", legend.show = F) + tm_layout("2010", title.size = 1.2, title.position = c('center','top'), title.color = "white")   + tm_shape(filter(Nubian, Year == 2010)) + tm_dots(col = "red", size = 0.1)
m4 <- tm_shape(GridZ) + tm_fill(col = "X2011_Dry", palette = "viridis", legend.show = F) + tm_layout("2011", title.size = 1.2, title.position = c('center','top'), title.color = "white")  + tm_shape(filter(Nubian, Year == 2011)) + tm_dots(col = "red", size = 0.1)
m5 <- tm_shape(GridZ) + tm_fill(col = "X2012_Dry", palette = "viridis", legend.show = F) + tm_layout("2012", title.size = 1.2, title.position = c('center','top'), title.color = "white")  + tm_shape(filter(Nubian, Year == 2012)) + tm_dots(col = "red", size = 0.1)
m6 <- tm_shape(GridZ) + tm_fill(col = "X2013_Dry", palette = "viridis", legend.show = F) + tm_layout("2013", title.size = 1.2, title.position = c('center','top'), title.color = "white")  + tm_shape(filter(Nubian, Year == 2013)) + tm_dots(col = "red", size = 0.1)
m7 <- tm_shape(GridZ) + tm_fill(col = "X2014_Dry", palette = "viridis", legend.show = F)  + tm_layout("2014", title.size = 1.2, title.position = c('center','top'), title.color = "white")  + tm_shape(filter(Nubian, Year == 2014)) + tm_dots(col = "red", size = 0.1)
m8 <- tm_shape(GridZ) + tm_fill(col = "X2015_Dry", palette = "viridis", legend.show = F) + tm_layout("2015", title.size = 1.2, title.position = c('center','top'), title.color = "white")  + tm_shape(filter(Nubian, Year == 2015)) + tm_dots(col = "red", size = 0.1)
m9 <- tm_shape(GridZ) + tm_fill(col = "X2016_Dry", palette = "viridis", legend.show = F) + tm_layout("2016", title.size = 1.2, title.position = c('center','top'), title.color = "white")  + tm_shape(filter(Nubian, Year == 2016)) + tm_dots(col = "red", size = 0.1)
m10 <- tm_shape(GridZ) + tm_fill(col = "X2017_Dry", palette = "viridis", legend.show = F) + tm_layout("2017", title.size = 1.2, title.position = c('center','top'), title.color = "white")   + tm_shape(filter(Nubian, Year == 2017)) + tm_dots(col = "red", size = 0.1)
m11 <- tm_shape(GridZ) + tm_fill(col = "X2017_Dry", palette = "viridis", title = "Occupancy probability", legend.reverse = T) + tm_shape(Nubian) + tm_dots(col = "red", size = 0.1) + tm_add_legend("symbol", shape = 16, size = 0.5, col = "red", labels = "Occurrence") + tm_layout(legend.only = T)


data("World")
bbox <- st_bbox(st_buffer(Grid, 550000))

Sys.setenv(MAPBOX_API_KEY="pk.eyJ1IjoicmFtaXJvY3JlZ28iLCJhIjoiY2tjbHg2MDd2MWxiOTJxbzNwdnZ4cW8yciJ9.uCHnVQRN8GtQAkoRQbUwKw")
basemap <- cc_location(loc = c(12.78384, 14.69689), buffer = 1500000, max_tiles = 4, type = "mapbox.satellite")
studyarea <- st_as_sf(st_as_sfc(st_bbox(Grid)))
studyarea$Name <- c("Study area")

m12 <- tm_shape(basemap, bbox = bbox) + tm_rgb() + tm_shape(World) + tm_borders() + tm_text("name", size = 0.7, ymod = 1) +
  tm_shape(studyarea) + tm_fill(col = 'grey', alpha = 0.5, title = "Occupancy probability") + tm_shape(Tracks) + tm_lines() + tm_scale_bar() + tm_add_legend("symbol", shape = 15, size = 0.5, col = "grey", labels = "Study area") + tm_add_legend("line", size = 0.5, col = "black", labels = "Survey tracks") + tm_layout(legend.position = c("right","top"), legend.bg.color = "white", legend.bg.alpha = 0.9, legend.text.size = 0.4)

map <- tmap_arrange(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, ncol = 3)

tmap_save(map, "./Figures/Fig4.jpg", width=5, height=7.51, units = "in")

#Display Tracks and Species data
msd1 <- tm_shape(GridSD) + tm_fill(col = "X2008_Dry", palette = "viridis", legend.show = F) + tm_layout("2008", title.size = 1.2, title.position = c('center','top'), title.color = "white")+ tm_shape(filter(Nubian, Year == 2008)) + tm_dots(col = "red", size = 0.1)
msd2 <- tm_shape(GridSD) + tm_fill(col = "X2009_Dry", palette = "viridis", legend.show = F) + tm_layout("2009", title.size = 1.2, title.position = c('center','top'), title.color = "white") + tm_shape(filter(Nubian, Year == 2009)) + tm_dots(col = "red", size = 0.1)
msd3 <- tm_shape(GridSD) + tm_fill(col = "X2010_Dry", palette = "viridis", legend.show = F) + tm_layout("2010", title.size = 1.2, title.position = c('center','top'), title.color = "white") + tm_shape(filter(Nubian, Year == 2010)) + tm_dots(col = "red", size = 0.1)
msd4 <- tm_shape(GridSD) + tm_fill(col = "X2011_Dry", palette = "viridis", legend.show = F) + tm_layout("2011", title.size = 1.2, title.position = c('center','top'), title.color = "white")  + tm_shape(filter(Nubian, Year == 2011)) + tm_dots(col = "red", size = 0.1)
msd5 <- tm_shape(GridSD) + tm_fill(col = "X2012_Dry", palette = "viridis", legend.show = F) + tm_layout("2012", title.size = 1.2, title.position = c('center','top'), title.color = "white")  + tm_shape(filter(Nubian, Year == 2012)) + tm_dots(col = "red", size = 0.1)
msd6 <- tm_shape(GridSD) + tm_fill(col = "X2013_Dry", palette = "viridis", legend.show = F) + tm_layout("2013", title.size = 1.2, title.position = c('center','top'), title.color = "white")  + tm_shape(filter(Nubian, Year == 2013)) + tm_dots(col = "red", size = 0.1)
msd7 <- tm_shape(GridSD) + tm_fill(col = "X2014_Dry", palette = "viridis", legend.show = F)  + tm_layout("2014", title.size = 1.2, title.position = c('center','top'), title.color = "white")  + tm_shape(filter(Nubian, Year == 2014)) + tm_dots(col = "red", size = 0.1)
msd8 <- tm_shape(GridSD) + tm_fill(col = "X2015_Dry", palette = "viridis", legend.show = F) + tm_layout("2015", title.size = 1.2, title.position = c('center','top'), title.color = "white") + tm_shape(filter(Nubian, Year == 2015)) + tm_dots(col = "red", size = 0.1)
msd9 <- tm_shape(GridSD) + tm_fill(col = "X2016_Dry", palette = "viridis", legend.show = F) + tm_layout("2016", title.size = 1.2, title.position = c('center','top'), title.color = "white") + tm_shape(filter(Nubian, Year == 2016)) + tm_dots(col = "red", size = 0.1)
msd10 <- tm_shape(GridSD) + tm_fill(col = "X2017_Dry", palette = "viridis", legend.show = F) + tm_layout("2017", title.size = 1.2, title.position = c('center','top'), title.color = "white") + tm_shape(filter(Nubian, Year == 2017)) + tm_dots(col = "red", size = 0.1)
msd11 <- tm_shape(GridSD) + tm_fill(col = "X2017_Dry", palette = "viridis", title = "Standard Deviation", legend.reverse = T) + tm_shape(Nubian) + tm_dots(col = "red", size = 0.1) + tm_add_legend("symbol", shape = 16, size = 0.5, col = "red", labels = "Occurrence") + tm_layout(legend.only = T)

mapSD <- tmap_arrange(msd1, msd2, msd3, msd4, msd5, msd6, msd7, msd8, msd9, msd10, msd11, m12, ncol = 3)

tmap_save(mapSD, "./Figures/NubianOccupancySD.jpg", width=5, height=7.51, units = "in")

# Plot spatial effect
# Arrange S
srows <- grepl(colnames(mcparams), pattern = "S")
Sss <- mcparams[,srows]

Smean <- apply(Sss, 2, mean)
hist(Smean)
mean(Smean)


GridS <- cbind(Grid, Smean)
tm_shape(GridS) + tm_fill(col = "Smedian", palette = "viridis", legend.show = T)

mapS <- tm_shape(GridS) + tm_fill(col = "Smean", palette = "viridis", legend.show = T, title = "Spatial effect", title.size = 0.6) + tm_layout( legend.title.size=1, legend.text.size = 0.7, legend.position = c("right","bottom"), legend.bg.color = "white", outer.bg.color = "white", legend.bg.alpha = 1)

tmap_save(mapS, "./Figures/SpatialEffectMap.jpg", width=5, height=5, units = "in")
