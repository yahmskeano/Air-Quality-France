library(geodata)
library(MODIS)
library(rgdal)
library(sp)
library(gdalUtils)
library(stringr)

france = st_read("Data/fr.geojson")
france <- subset(france, !nom %in% c("Guadeloupe","Martinique", "Guyane", "La Réunion","Mayotte" ))


#####################################
## Create AOD Rasters
#####################################
# (1) Store each hdf as raster
# (2) Average all raster bands for daily mean
# (3) Merge rasters by day 
# (4) Crop rasters to just France
# (5) Save final averaged raster
# 
# for (julian in 1:366){
#   print(Sys.time())
#   date = as.Date(julian-1, origin=as.Date("2020-01-01"))
#   files = list.files("AOD", pattern= paste0(".A2020", str_pad(julian, 3, pad="0"), "."))
# 
#   means = list()
# 
#   for (i in 1:length(files)){
#     sds = get_subdatasets(paste0("AOD/", files[i]))
#     gdal_translate(sds[1], dst_dataset = "AOD/aod_55.tif")
# 
#     rast <- stack("AOD/aod_55.tif")
#     mn = raster::calc(rast, fun = mean, na.rm=TRUE)
# 
#     means[[i]] = mn
#   }
# 
#   merged = means[[1]]
# 
#   for(i in 2:length(files)){
#     merged = mosaic(merged, means[[i]], fun=mean)
#   }
#   merged = projectRaster(merged, crs = CRS(SRS_string = "EPSG:4326"), method="ngb")
#   merged = crop(merged, france)
#   writeRaster(merged, filename=paste0("AOD/MergedAOD/47", date,".grd"))
# }


#####################################
## Create NA Raster
#####################################
# rast = raster("AOD/MergedAOD/472020-01-01.grd")
# 
# comp_rate = vector()
# 
# na_rast = !is.na(rast)
# comp_rate = c(comp_rate, sum(!is.na(getValues(rast)))/prod(dim(rast))) 
# 
# dates = seq(as.Date("2020-01-01"), as.Date("2020-12-31"), 1)
# 
# for (i  in 2:366){
#   
#   date = dates[i]
#   rast = raster(paste0("AOD/MergedAOD/47", date, ".grd"))
#   na_rast = na_rast + (!is.na(rast))
#   comp_rate = c(comp_rate, sum(!is.na(getValues(rast)))/prod(dim(rast))) 
# }


###############################################
## Upsample to 5km grid cells
###############################################

# blank_rast <- raster("blankRaster.grd")

# r_stack = stack()
# 
# for (julian in 1:366){
#   rast_date = as.Date(julian-1, origin=as.Date("2020-01-01"))
#   rast <- raster(paste0("AOD/MergedAOD/47",rast_date,".grd"))
#   rast <- raster::aggregate(rast, fact=5, fun=mean, na.rm=TRUE)
#   r_stack = stack(r_stack, rast)
# }
# 
# writeRaster(r_stack, paste0("AOD/47Stack.grd"))




###############################################
## Train Test sample 
###############################################


filt = readRDS("pm_data.rds")
locations = filt %>% distinct(loc, .keep_all=TRUE) %>% dplyr::select(loc, lat, lon)
loc = SpatialPoints(coords = locations[,3:2], proj4string = CRS("EPSG:4326"))
loc2 = spTransform(loc, CRS("EPSG:2154"))

final <- data.frame()
for (julian in 1:366){
  rast_date = as.Date(julian-1, origin=as.Date("2020-01-01"))
  rast <- raster(paste0("AOD/MergedAOD/",rast_date,".grd"))
  rast <- projectRaster(rast, crs = CRS("EPSG:2154"))
  aod47 = raster::extract(rast, loc2, method = "bilinear")
  temp_df = data.frame(aod47 = aod47,
                       lon = coordinates(loc)[,1],
                       lat = coordinates(loc)[,2],
                       date = rep(rast_date, length(loc)))
  final = rbind(final, temp_df)
}

filt_upd = left_join(filt, final, by=c("lat" = "lat", "lon" = "lon", "date"="date"))


#
# 
# 
# r_stack <- stack("AOD/47Stack.grd")
# r_stack <- projectRaster(r_stack, crs = CRS("EPSG:2154"))
# aod47 = raster::extract(r_stack, loc2, method = "bilinear")
# five_km = data.frame(aod47_5km = as.vector(aod47),
#                      lon = rep(coordinates(loc)[,1],366),
#                      lat = rep(coordinates(loc)[,2], 366),
#                      date = sort(rep(seq(as.Date("2020-01-01"), as.Date("2020-12-31"),1), length(loc))))
# 
# filt_upd = left_join(filt_upd, five_km, by=c("lat" = "lat", "lon" = "lon", "date"="date"))
# saveRDS(filt_upd,"pm_data.rds")
# 



###########################################################################
# Prediction Grid
###########################################################################

grid = readRDS("prediction_grid.rds")

locations = grid %>% distinct(lonL, latL, .keep_all=TRUE) %>% dplyr::select(lonL, latL)
loc = SpatialPoints(coords = locations, proj4string = CRS("EPSG:2154"))

r_stack <- stack("AOD/47Stack.grd")
r_stack <- projectRaster(r_stack, crs = CRS("EPSG:2154"))
r_stack = raster::resample(r_stack, blank_rast, method = "bilinear", na.rm=T)
temp_df = raster::extract(r_stack, loc, method="bilinear", na.rm=T)
temp_df = as.data.frame(temp_df)

colnames(temp_df) = c(seq(as.Date("2020-01-01"), as.Date("2020-12-31"), by=1))
temp_df$lonL = locations$lonL
temp_df$latL = locations$latL
temp_df = temp_df %>% pivot_longer(1:366)
colnames(temp_df)[3:4] = c("date", "aod47")
temp_df$date = as.Date(temp_df$date)

grid= left_join(grid, temp_df, by=c("lonL" = "lonL", "latL" = "latL", "date" = "date"))

saveRDS(grid,"prediction_grid.rds")



r_stack <- stack("AOD/CAMS_aod469.gri")
r_stack <- projectRaster(r_stack, crs = CRS("EPSG:2154"))

temp_df = raster::extract(r_stack, loc, method="bilinear", na.rm=T)
temp_df = as.data.frame(temp_df)

night_df = temp_df[,367:732]
temp_df = temp_df[,1:366]
colnames(temp_df) = c(seq(as.Date("2020-01-01"), as.Date("2020-12-31"), by=1))
colnames(night_df) = c(seq(as.Date("2020-01-01"), as.Date("2020-12-31"), by=1))
temp_df$lonL = locations$lonL
temp_df$latL = locations$latL
night_df$lonL = locations$lonL
night_df$latL = locations$latL
temp_df = temp_df  %>% pivot_longer(1:366)
night_df = night_df  %>% pivot_longer(1:366)
colnames(temp_df)[3:4] = c("date", "aod47_CAMS")
colnames(night_df)[3:4] = c("date", "aod47_CAMS_n")
temp_df$date = as.Date(temp_df$date)
night_df$date = as.Date(night_df$date)

grid= left_join(grid, temp_df, by=c("lonL" = "lonL", "latL" = "latL", "date" = "date"))
grid= left_join(grid, night_df, by=c("lonL" = "lonL", "latL" = "latL", "date" = "date"))

grid$aod47_imp[is.na(grid$aod47_imp)] = rowMeans(grid[is.na(grid$aod47_imp),7:8])
saveRDS(grid,"prediction_grid.rds")
