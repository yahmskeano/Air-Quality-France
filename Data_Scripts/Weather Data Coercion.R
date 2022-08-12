library(ncdf4)
library(ncdf4.helpers)
library(PCICt)
library(sf)


vars = c("rr", "hu", "tg", "pp", "fg")
names = c("precipitation", "humidity", "temp", "pressure", "wind")

blank_rast <- raster("blankRaster.grd")

france = st_read("Data/fr.geojson")
france <- subset(france, !nom %in% c("Guadeloupe","Martinique", "Guyane", "La Réunion","Mayotte" ))

for (i in 1:5){
  cur_var = vars[i]
  filepath = paste0("Weather/", cur_var, "_ens_mean_0.1deg_reg_v24.0e.nc")
  nc_data <-nc_open(filepath)
  
  
  lons = ncvar_get(nc_data, varid="longitude")
  
  lats = ncvar_get(nc_data, varid="latitude")
  
  times = ncvar_get(nc_data, varid="time")
  
  dates = as.Date(times, origin = as.Date("1950-01-01"))
  
  
  lon_index <- which(lons >= extent(france)[1] & lons <= extent(france)[2])
  
  lat_index <- which(lats >= extent(france)[3] & lats <= extent(france)[4])
  
  time_index <- which(dates >= "2020-01-01" & dates <= "2020-12-31")
  
  time_series <- nc.get.time.series(nc_data, v = cur_var,
                                    time.dim.name = "time")
  
  values <- nc.get.var.subset.by.axes(nc_data, cur_var,
                                      axis.indices = list(X = lon_index,
                                                          Y = lat_index,
                                                          T = time_index))
  
  df <- expand.grid(lons[lon_index], lats[lat_index], dates[time_index]) %>% 
    rename(lon = Var1, lat = Var2, date = Var3) %>%
    mutate(value = as.vector(values))

  df$value[df$value < -99] = NA
  df$date = as.Date(df$date)
  r_stack = stack()
  
  for (j in 1:length(dates[time_index])){
    rast_date = dates[time_index][j]
    rast <- rasterFromXYZ(df %>% filter(date == rast_date) %>% select(-date))
    rast <- raster::resample(rast, blank_rast, method="bilinear")
    rast <- mask(rast, france)
    r_stack = stack(r_stack, rast)
  }
  names(r_stack) = dates[time_index]
  writeRaster(r_stack, paste0("Weather/", names[i],".stk"))
  print(paste0(names[i], " done \n"))
}



filt = readRDS("pm_data.rds")
locations = filt %>% distinct(loc, .keep_all=TRUE) %>% dplyr::select(loc, lat, lon)
loc = SpatialPoints(coords = locations[,3:2], proj4string = CRS("EPSG:4326"))
loc2 = spTransform(loc, CRS("EPSG:2154"))



for (i in 1:5){
  cur_var = names[i]
  
  
  r_stack <- stack(paste0("Weather/", cur_var,".grd"))
  r_stack <- projectRaster(r_stack, crs = CRS("EPSG:2154"))
  weather = raster::extract(r_stack, loc2, method = "bilinear")
  temp = data.frame(value = as.vector(weather),
                       lon = rep(coordinates(loc)[,1],366),
                       lat = rep(coordinates(loc)[,2], 366),
                       date = sort(rep(seq(as.Date("2020-01-01"), as.Date("2020-12-31"),1), length(loc))))
  
  filt = readRDS("pm_data.rds")
  colnames(temp)[1] = cur_var
  filt_upd = left_join(filt, temp, by=c("lat" = "lat", "lon" = "lon", "date"="date"))
  saveRDS(filt_upd,"pm_data.rds")
  
}



grid = readRDS("prediction_grid.rds")
locations = grid %>% distinct(lonL, latL, .keep_all=TRUE) %>% dplyr::select(lonL, latL)
loc = SpatialPoints(coords = locations, proj4string = CRS("EPSG:2154"))

for (i in 2:5){
  cur_var = names[i]
  
  
  r_stack <- stack(paste0("Weather/", cur_var,".grd"))
  r_stack <- projectRaster(r_stack, crs = CRS("EPSG:2154"))
  temp_df = raster::extract(r_stack, loc, method="bilinear", na.rm=T)
  temp_df = as.data.frame(temp_df)
  
  colnames(temp_df) = c(seq(as.Date("2020-01-01"), as.Date("2020-12-31"), by=1))
  temp_df$lonL = locations$lonL
  temp_df$latL = locations$latL
  temp_df = temp_df %>% pivot_longer(1:366)
  colnames(temp_df)[3:4] = c("date", cur_var)
  temp_df$date = as.Date(temp_df$date)
  
  grid= left_join(grid, temp_df, by=c("lonL" = "lonL", "latL" = "latL", "date" = "date"))
  
}