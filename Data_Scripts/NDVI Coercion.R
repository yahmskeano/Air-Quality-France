library(ncdf4)
library(ncdf4.helpers)
library(PCICt)
library(tidyverse)
library(raster)
library(dtplyr)
library(sf)


blank_rast <- raster("blankRaster.grd")
france = st_read("Data/fr.geojson")
france <- subset(france, !nom %in% c("Guadeloupe","Martinique", "Guyane", "La Réunion","Mayotte" ))


filename = list.files("NDVI")

nc_data <-nc_open(paste0("NDVI/",filename))

lons = ncvar_get(nc_data, varid="lon")
lats = ncvar_get(nc_data, varid="lat")
times = as.Date(ncvar_get(nc_data, varid="time"), origin = as.Date("2000-01-01"))
values = ncvar_get(nc_data, varid="_1_km_16_days_NDVI")

df <- expand.grid(lons, lats, times) %>% 
  rename(lon = Var1, lat = Var2, date = Var3) %>%
  mutate(value = as.vector(values))

nc_close(nc_data)

#########################################
# Loop for creating final prediction grid
#########################################
# r_stack <- stack()
# for (j in 1:length(times)){
#   rast_date = times[j]
#   rast <- rasterFromXYZ(df %>% filter(date == rast_date) %>% dplyr::select(-date), crs="+proj=longlat +datum=WGS84 +no_defs ")
#   rast <- aggregate(rast, fact=5, fun=mean, na.rm=T)
#   rast <- projectRaster(rast, crs=CRS("EPSG:2154"))
#   rast <- raster::resample(rast, blank_rast, method="bilinear")
#   rast <- mask(rast, france)
#   r_stack = stack(r_stack, rast)
# }
# 
# names(r_stack) = times
# 
# writeRaster(r_stack, "NDVI/ndvi_stack.grd", overwrite=TRUE)

#################################################
# Loop for creating train/test set
#################################################
# 
# filt = readRDS("pm_data.rds")
# locations = filt %>% distinct(loc, .keep_all=TRUE) %>% dplyr::select(loc, lat, lon)
# loc = SpatialPoints(coords = locations[,3:2], proj4string = CRS("EPSG:4326"))
# loc2 = spTransform(loc, CRS("EPSG:2154"))
# 
# final <- data.frame()
# for (j in 1:length(times)){
#   rast_date = times[j]
#   if (j < length(times)){
#     dates = seq(times[j], times[(j+1)]-1, 1)
#   }
#   else{
#     dates = seq(times[j], as.Date("2020-12-31"), 1)
#   }
#   n_dates = length(dates)
#   rast <- rasterFromXYZ(df %>% filter(date == rast_date) %>% dplyr::select(-date), crs="+proj=longlat +datum=WGS84 +no_defs ")
#   rast <- projectRaster(rast, crs = CRS("EPSG:2154"))
#   ndvi = raster::extract(rast, loc2, method = "bilinear")
#   temp_df = data.frame(ndvi = ndvi,
#                        lon = rep(coordinates(loc)[,1], n_dates),
#                        lat = rep(coordinates(loc)[,2], n_dates),
#                        date = sort(rep(dates, length(loc))))
#   final = rbind(final, temp_df)
# }
# 
# filt = left_join(filt, final, by=c("lat" = "lat", "lon" = "lon", "date"="date"))
# saveRDS(filt,"pm_data.rds")



##############################################################################
# Prediction grid building
##############################################################################

# grid = readRDS("prediction_grid.rds")
# ncells=  nrow(distinct(grid, lonL,latL))
# 
# final <- data.frame()
# for (j in 1:length(times)){
#   date = times[j]
#   if (j < length(times)){
#     dates = seq(times[j], times[(j+1)]-1, 1)
#   }
#   else{
#     dates = seq(times[j], as.Date("2020-12-31"), 1)
#   }
#   temp_df = grid %>% filter(date == times[j])
#   temp_df = do.call("rbind", replicate(length(dates), temp_df, simplify = FALSE))
#   temp_df$dates = sort(rep(dates, ncells))
#   final = rbind(final, temp_df)
# }
# 
# saveRDS(grid, file="prediction_grid.rds")
