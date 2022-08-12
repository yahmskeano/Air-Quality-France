library(ncdf4)
library(ncdf4.helpers)
library(PCICt)
library(stringr)
library(raster)


vars = c("aod469", "aod550", "aod670", "aod865", "aod1240")

filepath = paste0("AOD/cams.nc")
nc_data <-nc_open(filepath)
blank_rast <- raster("blankRaster.grd")

lons = ncvar_get(nc_data, varid="longitude")

lats = ncvar_get(nc_data, varid="latitude")

for (i in 1:5){
  cur_var = vars[i]
  time_series <- nc.get.time.series(nc_data, v = cur_var,
                                    time.dim.name = "time")
  dates <- as.Date(substring(as.character(time_series),1,10))
  
  values <- ncvar_get(nc_data, cur_var)
  
  df <- expand.grid(lons, lats, dates) %>% 
    rename(lon = Var1, lat = Var2, date = Var3) %>%
    mutate(value = as.vector(values))

  df$value[df$value < -99] = NA
  df$time = rep(c("midday", "midnight"), nrow(df)/2)
  r_stack = stack()
  dates <- unique(dates)
  for (j in 1:length(dates)){
    rast_date = dates[j]
    rast <- rasterFromXYZ(df %>% filter(date == rast_date, time == "midday") %>% select(-date, -time), digits=4)
    rast <- raster::resample(rast, blank_rast, method="bilinear")
    r_stack = stack(r_stack, rast)
  }
  for (j in 1:length(dates)){
    rast_date = dates[j]
    rast <- rasterFromXYZ(df %>% filter(date == rast_date, time == "midnight") %>% select(-date, -time), digits=4)
    rast <- raster::resample(rast, blank_rast, method="bilinear")
    r_stack = stack(r_stack, rast)
  }
  names(r_stack) = c(dates, dates)
  writeRaster(r_stack, paste0("AOD/CAMS_", cur_var,".stk"), overwrite=TRUE)
  print(paste0(cur_var, " done \n"))
  
}

nc_close(nc_data)

#########################################
## Test/Train
########################################

filt = readRDS("pm_data.rds")
locations = filt %>% distinct(loc, .keep_all=TRUE) %>% dplyr::select(loc, lat, lon)
loc = SpatialPoints(coords = locations[,3:2], proj4string = CRS("EPSG:4326"))
loc2 = spTransform(loc, CRS("EPSG:2154"))



for (i in 1:5){
  cur_var = vars[i]
  

  r_stack <- stack(paste0("AOD/CAMS_", cur_var,".grd"))
  r_stack <- projectRaster(r_stack, crs = CRS("EPSG:2154"))
  aod = raster::extract(r_stack, loc2, method = "bilinear")
  five_km = data.frame(aod = as.vector(aod[,1:366]),
                       aod2 = as.vector(aod[,367:732]),
                       lon = rep(coordinates(loc)[,1],366),
                       lat = rep(coordinates(loc)[,2], 366),
                       date = sort(rep(seq(as.Date("2020-01-01"), as.Date("2020-12-31"),1), length(loc))))
  
  filt = readRDS("pm_data.rds")
  colnames(five_km)[1:2] = c(paste0("CAMS_", cur_var, "_day"), paste0("CAMS_", cur_var, "_night"))
  filt_upd = left_join(filt, five_km, by=c("lat" = "lat", "lon" = "lon", "date"="date"))
  saveRDS(filt_upd,"pm_data.rds")

}