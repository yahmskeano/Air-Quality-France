library(geodata)
library(MODIS)
library(rgdal)
library(sp)
library(sf)
library(gdalUtils)
library(tidyverse)
library(foreach)
library(doParallel)

france = st_read("Data/fr.geojson")
france <- subset(france, !nom %in% c("Guadeloupe","Martinique", "Guyane", "La Réunion","Mayotte" ))
france = st_transform(france, projection(rast))


files = list.files("Landcover")
names = c("water", "trees", "grass", "veg", "crops", "shrub", "built", "bare", "snow", "clouds")

################################################################################
## Aggregate tiles from 10mx10m to 5kmx5km
################################################################################
# 
# for(i in 2:6) {
#   file = paste0("Landcover/",files[i])
#   rast = raster(file)
#   rast = crop(rast, france)
#   print("cropping done")
#   
#   UseCores <- 10
#   cl       <- makeCluster(UseCores)
#   registerDoParallel(cl)
#   
#   foreach (j = 1:10) %dopar% {
#     library(raster)
#     rast_temp = mask(rast, rast == j , maskvalue = 0)
#     rast_temp <- aggregate(x = rast_temp, fact = 500, fun = sum, na.rm = T)
#     rast_temp = rast_temp / (250000 * j)
#     writeRaster(rast_temp, paste0("Landcover/Reduced/Raster_",i,"_",names[j]))
#   }
#   print(Sys.time())
#   print(i)
#   stopCluster(cl)
#   gc()
# }


################################################################################
## Merge tiles and map to Metro France blocks
################################################################################

blank_rast <- raster("blankRaster.grd")


for (i in 1:10){
  rast_list = list()
  for (j in 1:6){
    rast = raster(paste0("Landcover/Reduced/Raster_", j, "_", names[i],".gri"))
    rast = projectRaster(rast, crs = projection(blank_rast), method="ngb", na.rm=FALSE)
    rast = resample(rast, blank_rast, method = "ngb", na.rm=FALSE)
    rast_list[[j]] = rast
  }
  full_rast = mosaic(rast_list[[1]],
                     rast_list[[2]],
                     rast_list[[3]],
                     rast_list[[4]],
                     rast_list[[5]],
                     rast_list[[6]],
                     fun = max,
                     na.rm = TRUE)
  full_rast = projectRaster(full_rast, crs=CRS("EPSG:2154"), method="ngb", na.rm=T)
  writeRaster(full_rast, filename = paste0("Landcover/", names[i],".gri"), overwrite=T)
}




###############################################
## Train Test sample 
###############################################


filt = readRDS("pm_data.rds")
locations = filt %>% distinct(loc, .keep_all=TRUE) %>% dplyr::select(loc, lat, lon)
loc = SpatialPoints(coords = locations[,3:2], proj4string = CRS("EPSG:4326"))
loc2 = spTransform(loc, CRS("EPSG:2154"))


for (i in 1:10){
  rast = raster(paste0("Landcover/", names[i],".gri"))
  
  values = raster::extract(rast, loc2, method = "bilinear")
  locations$value = values
  colnames(locations)[(3+i)] = paste0(names[i],"_prop")
}
locations$veg_prop = locations$veg_prop + locations$shrub_prop + locations$grass_prop + locations$bare_prop
locations = locations[,-c(6,9,11,12,13)]
filt = left_join(filt, locations, by=c("lon" = "lon", "lat" = "lat"))

saveRDS(filt, "pm_data.rds")


grid = readRDS("prediction_grid.rds")

locations = grid %>% distinct(lonL, latL, .keep_all=TRUE) %>% dplyr::select(lonL, latL)
loc = SpatialPoints(coords = locations, proj4string = CRS("EPSG:2154"))


for (i in 1:10){
  rast = raster(paste0("Landcover/", names[i],".gri"))
  
  values = raster::extract(rast, loc, method = "bilinear")
  locations$value = values
  colnames(locations)[(2+i)] = paste0(names[i],"_prop")
}
grid = left_join(grid, locations, by=c("lonL" = "lonL", "latL" = "latL"))
