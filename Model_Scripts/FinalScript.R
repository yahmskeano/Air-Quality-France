library(ape)
library(readr)
library(dplyr)
library("sp")
library("spacetime")
library("animation")
library("ggplot2")
library("maps")
library("STRbook")
library("grid")
library("gridExtra")
library(lubridate)
library(skimr)
library(zoo)
library(tidyverse)
library(ggthemes)
library(sf)
library(rgdal)
library(RcmdrMisc)
library("gstat")
library(TSA)
library(forecast)
library(caret)
library(parallel)
library(viridis)
library(grid)
source("GstatFunctions.R")
source("Prediction_Functions.R")

plot_preds <- function(preds, data){
  truth = data$value
  rmse = round(sqrt(get_bt_mse(preds@data$bt_preds, data)),3)
  nnas = !is.na(truth)
  truth = truth[nnas]
  preds = preds$bt_preds[nnas]
  df = data.frame("preds" = preds, "actual" = truth)
  
  title = paste0("Cor: ", round(cor(preds, truth),3), ", RMSE: ", rmse)
  g1 <- ggplot(df, aes(x=actual, y=preds)) + geom_point(alpha=.1, color="blue") + 
    theme_tufte() + labs(x = "Observed", y = "Predicted") +
    geom_abline(slope=1, intercept = 0, color="grey") +theme_hc() + 
    ylim(0, 70)
  
  g1
  
  # g2 <- ggplot(df) + geom_histogram(aes(x= preds - actual), bins=75) + theme_tufte() + labs(x = "Residuals", y="") 
  # 
  # grid.arrange(g1, g2, nrow=1, top = textGrob(title ,gp=gpar(fontsize=20,font=3)), widths=2:1)
}

get_Data <- function(data, locations, type){
  crs_wgs84 <- CRS(SRS_string = "+init=EPSG:2154")
  wkt_wgs84 <- wkt(crs_wgs84)
  
  temp_part <- seq(as.Date(min(data$date)), by = "days", length.out = length(unique(data$date)))
  
  spat_part <- SpatialPoints(coords=  locations[,c("lon","lat")],proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  spat_part = spTransform(spat_part, CRS("+init=epsg:2154"))
  
  data[,c("lonL", "latL")] = data[,c("lonL", "latL")]/10000 
  long_data <- data %>% arrange(date,loc)
  
  
  
  STData <- STFDF(sp = spat_part,
                  time = temp_part,
                  data = long_data)
  
  proj4string(STData) <-  crs_wgs84
  
  return(STData)
  
}


get_parts <- function(start, end, n_days, n_cores){
  dates = seq(start, end, n_days/n_cores)
  parts = paste0(dates, "::", c(dates[-1]-1, end))
  return(parts)
}




get_preds <- function(trn, val, covs, formula, nmax){
  errors = list()
  mses = vector()
  mses_trans = vector()
  maes = vector()
  sds = vector()
  sdsds = vector()
  trn = as(trn, "STIDF") 
  trn = subset(trn, !is.na(trn$value))
  
  
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  parts <- get_parts(min(trn@data$date), max(trn@data$date), max(trn@time), no_cores)
  
  
  
  clusterExport(cl = cl, varlist = c("formula", "trn", "val", "covs", "parts", "nmax"), envir = .GlobalEnv)
  clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
  
  parallelX <- parLapply(cl = cl, X = 1:no_cores, fun = function(x) krigeST(formula , # latitude trend
                                                                            data = trn, # data set w/o 14 July
                                                                            newdata = val[, parts[[x]]], # prediction grid
                                                                            modelList = covs, # semivariogram
                                                                            nmax = nmax,
                                                                            stAni = covs$stAni,
                                                                            nbuffer=2,
                                                                            computeVar = TRUE,
                                                                            progress = TRUE))
  stopCluster(cl)
  
  preds = rbind(parallelX[[1]]@data, parallelX[[2]]@data)
  p_dates = c(index(parallelX[[1]]@time),index(parallelX[[2]]@time))
  for (i in 3:no_cores){
    preds = rbind(preds, parallelX[[i]]@data)
    p_dates = c(p_dates, index(parallelX[[i]]@time))
  }
  
  ST_pred = STFDF(sp = val@sp,
                  time = p_dates,
                  data = preds)
    
  mses = get_mse(ST_pred$var1.pred, val@data)
  ST_pred@data$bt_preds = exp(ST_pred@data$var1.pred + ST_pred@data$var1.var/2)
  ST_pred@data$bt_var = (exp(ST_pred@data$var1.var) - 1) * exp(2*ST_pred@data$var1.pred + ST_pred@data$var1.var/2)
    
  errors[[1]] = mses
  errors[[2]] = ST_pred
  print(paste0("Predictions done ", Sys.time()))
  
  return(errors)
}




filt = readRDS("pm_data.rds")
locations = filt %>% distinct(loc, .keep_all=TRUE) %>% dplyr::select(loc, lat, lon)
filt$trans_value = log(filt$value)

#filt = filt %>% filter(date < "2020-04-01")

filt$month = as.factor(month(filt$date))
filt$weekday = weekdays(filt$date)
filt = filt[,-c(8:19)]

filt$elevation = log(filt$elevation)
filt[,15:16] = sqrt(filt[,15:16])
filt[,19:22] = sqrt(filt[,19:22])

## Train-Test Split
set.seed(111) #formerly 111
# closest = c(8, 24, 29, 32, 38, 43, 49, 54, 63, 79, 84, 83, 91, 98)
# options = 1:nrow(locations)
# options = options[-closest]
# test_samp = c(sample(options, size=12, replace=FALSE), closest)
test_samp = sample(1:nrow(locations), size=26, replace=FALSE)
Test = filt %>% filter(loc %in% locations$loc[test_samp])



train_loc = locations[-test_samp,]

Train = filt %>% filter(loc %in% train_loc$loc)

Train$resids = fitted(fit) - log(Train$value + 1)

STTrain = get_Data(Train, train_loc,"")
STTest = get_Data(Test, locations[test_samp,],"")

vv_nn <-  variogram(resids ~ 1,
                    data = STTrain, # July data # consider pts < 1000 km apart
                    cutoff = 4.2e5,
                    tlags = 0:7,
                    cores = 8,
                    na.omit=T) # 0 days to 6

covs = get_variograms(vv_nn)
get_cov_mse(covs)

plot(vv_nn, covs, all=T, wireframe=T)

trn = STTrain
trn@data = data.frame(trn@data)
val = STTest
val@data = data.frame(val@data)
nmax = 1000
covs = covs[[4]][[4]]
errors = get_preds(trn, val, covs, formula, nmax)
errors = rk_preds(trn, val, covs, formula, nmax, model)

plot_preds(errors[[5]], val)


de_pred=  readRDS("prediction_grid.rds")
preds = predict(glasso_fit, de_pred, se.fit=TRUE)

de_pred$preds = preds$fit
de_pred$lm_var = preds$se.fit^2


temp_pred_grid <- c(as.Date("2020-01-25"), as.Date("2020-02-29"), as.Date("2020-05-14"), as.Date("2020-07-14"), as.Date("2020-09-22"), as.Date("2020-11-27"))


de_pred_dat = de_pred %>% filter(date %in% temp_pred_grid)

de_pred_dat = de_pred_dat %>% arrange(date, lonL, latL)


sp = de_pred_dat %>% distinct(lonL, latL)
spat_pred_grid <- SpatialPoints(coords=  sp[,c("lonL","latL")],proj4string = CRS("EPSG:2154"))



crs_wgs84 <- CRS(SRS_string = "EPSG:2154")
wkt_wgs84 <- wkt(crs_wgs84)

DE_pred <- STFDF(sp = spat_pred_grid, # spatial part
                 time = temp_pred_grid,
                 data = de_pred_dat) # temporal part
proj4string(DE_pred) <-  crs_wgs84


get_pred_grid <- function(trn, val, covs, formula, nmax){
  errors = list()
  mses = vector()
  mses_trans = vector()
  maes = vector()
  sds = vector()
  sdsds = vector()
  trn = as(trn, "STIDF") 
  trn = subset(trn, !is.na(trn$value))
  
  parts = c(seq(1,23667, 1578), 23667)
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  
  
  
  clusterExport(cl = cl, varlist = c("formula", "trn", "val", "covs", "parts", "nmax"), envir = .GlobalEnv)
  clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
  
  parallelX <- parLapply(cl = cl, X = 1:no_cores, fun = function(x) krigeST(resids~1, # latitude trend
                                                                            data = trn, # data set w/o 14 July
                                                                            newdata = val[parts[x]:parts[x+1], ], # prediction grid
                                                                            modelList = covs, # semivariogram
                                                                            nmax = nmax,
                                                                            stAni = covs$stAni,
                                                                            nbuffer=2,
                                                                            computeVar = TRUE,
                                                                            progress = TRUE))
  stopCluster(cl)
  
  preds = rbind(parallelX[[1]]@data, parallelX[[2]]@data)
  p_dates = index(parallelX[[1]]@time)
  for (i in 3:no_cores){
    preds = rbind(preds, parallelX[[i]]@data)
  }
  preds = arrange(preds, date, lonL, latL)
  preds = distinct(preds)
  
  sp = preds %>% distinct(lonL, latL)
  spat_pred_grid <- SpatialPoints(coords=  sp[,c("lonL","latL")],proj4string = CRS("EPSG:2154"))
  ST_pred = STFDF(sp = spat_pred_grid,
                  time = p_dates,
                  data = preds)
  
  proj4string(ST_pred) <-  crs_wgs84

  ST_pred@data$rk_preds = ST_pred@data$preds - ST_pred@data$var1.pred
  ST_pred@data$var = ST_pred@data$var1.var + ST_pred@data$lm_var
  
  ST_pred@data$bt_preds = exp(ST_pred@data$rk_preds + ST_pred@data$var/2)
  ST_pred@data$bt_var = (exp(ST_pred@data$var) - 1) * exp(2*ST_pred@data$rk_preds + ST_pred@data$var/2)

  print(paste0("Predictions done ", Sys.time()))
  
  return(ST_pred)
}


trn = STTrain
trn@data = data.frame(trn@data)
val = DE_pred
val@data = data.frame(val@data)
grid_ok = get_pred_grid_ok(trn, val,covs, formula, nmax)
grid_rk = get_pred_grid(trn, val,covs, formula, nmax)

# idw formula = log(value) ~  elevation  + NDVI + idw_wind + idw_humidity + idw_precip + idw_temp + C1.366+S1.366+S1.7+C1.7 + nearest_maj + nearest_res + num_maj)

grid$trees_prop[(grid$elevation < 6) & (grid$trees_prop < 0.02) & (grid$lonL < 600000) & (grid$lonL > 500000) &(grid$latL <6600000)] = NA

grid$crops_prop[(grid$elevation < 6) & (grid$crops_prop < 0.02) & (grid$lonL < 600000) & (grid$lonL > 500000) &(grid$latL <6600000)] = NA

grid$built_prop[(grid$elevation < 6) & (grid$built_prop < 0.01) & (grid$lonL < 600000)& (grid$lonL > 500000) & (grid$latL <6600000)] = NA
