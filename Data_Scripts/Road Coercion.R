library(raster)
library(tidyverse)
library(dtplyr)
library(sf)
library(rgdal)
library(rgeos)
library(terra)
library(maptools)
blank_rast <- raster("blankRaster.grd")
france <- st_read("Data/fr.geojson")
france = st_transform(france, CRS("EPSG:2154"))
blank_rast <- projectRaster(blank_rast, crs=CRS("EPSG:2154"))


roads <- st_read("Data/roads.shp")
rds <- st_transform(roads, crs = CRS("+init=epsg:2154"))

#rds = rds %>% filter(type %in% c("primary", "trunk", "motorway", "secondary", "tertiary"))
rds = rds %>% filter(type == "residential")
rs <- raster(extent(rds), crs=projection(rds))
res(rs) = res(blank_rast)
rs[] <- 1:ncell(rs)
#rs = mask(rs, france)

rsp <- rasterToPolygons(rs)
rsp <- st_as_sf(rsp)
rsp = st_transform(rsp, st_crs(rds))
ints = st_intersection(rds, rsp)
ints$len = st_length(ints)
join = st_join(rsp, ints)
out = group_by(join, layer.y) %>%
  summarize(length = sum(len))

rast = rasterize(out, blank_rast, "length")
# writeRaster(rast,"Roads/hway_density.grd")
writeRaster(rast,"Roads/res_density.grd")






filt = readRDS("pm_data.rds")
locations = filt %>% distinct(loc, .keep_all=TRUE) %>% dplyr::select(loc, lat, lon)
loc = SpatialPoints(coords = locations[,3:2], proj4string = CRS("EPSG:4326"))
loc2 = spTransform(loc, CRS("EPSG:2154"))

names = c("res", "hway")

replaceNA <- function(x, na.rm, ...){ 
  if(is.na(x[1]))
    return(0)
  else
    return(x)
} 

for (i in 1:2){
  cur_var = names[i]
  
  
  r_stack <- stack(paste0("Roads/", cur_var,"_density.grd"))
  r_stack <- projectRaster(r_stack, crs = CRS("EPSG:2154"))
  r_stack <- calc(r_stack, fun = replaceNA)
  
  
  roads = raster::extract(r_stack, loc2, method = "bilinear")
  temp = data.frame(value = as.vector(roads),
                    lon = rep(coordinates(loc)[,1],366),
                    lat = rep(coordinates(loc)[,2], 366),
                    date = sort(rep(seq(as.Date("2020-01-01"), as.Date("2020-12-31"),1), length(loc))))
  
  filt = readRDS("pm_data.rds")
  colnames(temp)[1] = paste0(cur_var, "_density")
  filt_upd = left_join(filt, temp, by=c("lat" = "lat", "lon" = "lon", "date"="date"))
  saveRDS(filt_upd,"pm_data.rds")
  
}





rds <- st_transform(roads, crs = CRS("+init=epsg:2154"))
points = SpatialPoints(coords = locations[,c("lon","lat")], CRS("+init=epsg:4326"))
points = as(points, "sf")
points = st_transform(points, CRS("+init=epsg:2154"))


get_nearest <- function(points, roads, type_list){
  n = nrow(points)
  shortest.dists <- numeric(n)
  n_roads <- numeric(n)
  temp = roads %>% filter(type %in% type_list)
  
  retvals = list()
  
  for (i in seq_len(n)) {
    
    distances <- st_distance(points[i,], temp)
    shortest.dists[i] = min(distances)
    n_roads[i] = sum(as.numeric(distances) < 1000)
  }
  retvals[[1]] = shortest.dists
  retvals[[2]] = n_roads
  
  return(retvals)
}

residentials = get_nearest(points, rds, "residential")
print("res done")
hway1 = get_nearest(points, rds, c("primary", "trunk", "motorway"))
print("hwy done")
hway2 = get_nearest(points, rds, "secondary")
print("hwy done")
hway3 = get_nearest(points, rds, "tertiary")
print("hwy done")

locations$nearest_res = residentials[[1]]
locations$num_res_close =  residentials[[2]]
locations$nearest_maj = hway1[[1]]
locations$num_maj = hway1[[2]]
locations$nearest_sec = hway2[[1]]
locations$num_sec = hway2[[2]]

locations$nearest_ter = hway3[[1]]
locations$num_ter = hway3[[2]]

filt = left_join(filt, locations, by=c("loc" = "loc"))

saveRDS(filt, "pm_data.rds")
