library(geodata)
library(MODIS)
library(rgdal)
library(sp)
library(sf)
library(gdalUtils)
library(tidyverse)


rast = raster("PopDensity/pop_density.tif")
blank_rast <- raster("blankRaster.grd")

rast = projectRaster(rast, blank_rast, method="bilinear", na.rm=TRUE)

##############################################################################
# Train Test Sample Extract
##############################################################################
filt = readRDS("pm_data.rds")
locations = filt %>% distinct(loc, .keep_all=TRUE) %>% dplyr::select(loc, lat, lon)
loc = SpatialPoints(coords = locations[,3:2], proj4string = CRS("EPSG:4326"))
loc2 = spTransform(loc, CRS("EPSG:2154"))


values = raster::extract(rast, loc2, method = "bilinear")
locations$pop_density = values

locations = locations[,-c(2:3)]
filt = left_join(filt, locations, by=c("loc" = "loc"))

saveRDS(filt, "pm_data.rds")




###########################################################################
# Prediction Grid
###########################################################################

grid = readRDS("prediction_grid.rds")

locations = grid %>% distinct(lonL, latL, .keep_all=TRUE) %>% dplyr::select(lonL, latL)
loc = SpatialPoints(coords = locations, proj4string = CRS("EPSG:2154"))

rast = raster("PopDensity/pop_density.tif")

values = raster::extract(rast, loc, method = "bilinear")
locations$pop_density = values
grid = left_join(grid, locations, by=c("lonL" = "lonL", "latL" = "latL"))
