library(gstat)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(raster)

#####################################################
# Fitted Covariances
#####################################################

covs_rk = readRDS("thesis_output/covs_rk.rds")
vv = readRDS("thesis_output/vv_rk.rds")
plot(vv, list(covs_rk[[1]][[1]],covs_rk[[2]][[1]], covs_rk[[3]][[1]], covs_rk[[4]][[2]], covs_rk[[5]][[4]]), all=T, wireframe=T, xlab ="distance(m)", ylab="time (days)")

covs_ok = readRDS("thesis_output/covs_ok.rds")
vv_ok = readRDS("thesis_output/vv_ok.rds")
plot(vv_ok, list(covs_ok[[1]][[1]],covs_ok[[2]][[1]], covs_ok[[3]][[2]], covs_ok[[4]][[4]], covs_ok[[5]][[6]]), all=T, wireframe=T, xlab ="distance(m)", ylab="time (days)")



tseq <- seq(0,7,1)
grid <- expand.grid(sort(unique(vv$spacelag)),tseq)
dg <- data.frame(spacelag=grid[,1],timelag=grid[,2])
table_ok <- variogramSurface(covs[[4]][[2]],dg)
table_rk <- variogramSurface(covs_rk[[5]][[6]],dg)


df = table_ok
df$gamma_rk = table_rk$gamma
df$vv_ok = c(0,vv$gamma)
df$vv_rk = c(0,vv_rk$gamma)

df %>% filter(timelag==0) %>% ggplot(aes(x=spacelag)) + 
  geom_point(aes(y=vv_ok, color="Data")) + geom_line(aes(y=gamma+.01, color="Data")) +
  geom_point(aes(y=vv_rk, color="Residuals")) + geom_line(aes(y=gamma_rk, color="Residuals")) +
  scale_color_manual(name="Type", values=c("Data" = "green", "Residuals" = "red"))


filt = filt %>%
  group_by(loc) %>%
  mutate(value = case_when(is.na(value) ~ mean(value, na.rm=TRUE),
                           TRUE ~ as.numeric(value) 
  )
  )

filt = filt %>%
  group_by(loc) %>%
  mutate(trans_value = case_when(is.na(trans_value) ~ mean(trans_value, na.rm=TRUE),
                           TRUE ~ as.numeric(trans_value) 
  )
  )

mean_set = filt %>% filter(value <= 200, trans_value > -2.3) %>% group_by(date) %>% summarise(mn_val = mean(value, na.rm=T), mn_trans = mean(trans_value, na.rm=T))

g1 <- ggplot() + geom_line(data= filt %>%filter(value <= 125), aes(x=date, y=value, group=loc),alpha=0.1) + 
  theme_hc() +  labs(x="", y="PM2.5")
g2 <- ggplot() + geom_line(data= filt %>%filter(value <= 125, trans_value >-2.3), aes(x=date, y=log(value+1), group=loc),alpha=0.1) + 
  theme_hc() +  labs(x="Date", y="logPM2.5")
grid.arrange(g1, g2, ncol=1)

#####################################################
# NDVI
#####################################################

ndvi = readRDS("NDVI/reduced_ndvi_df.rds")
ggplot(ndvi) + 
  geom_line(aes(x=date, y=median), color="forestgreen") +
  geom_line(aes(x=date, y=q3), color="lightgreen") + 
  geom_line(aes(x=date, y=q1), color="lightgreen") +
  theme_hc() + labs(x="Date", y="NDVI") +
  ylim(0,1)


#####################################################
# Local kriging cross-validation
#####################################################

mses_ok = readRDS("thesis_output/mses_ok_better.rds")
mses_rk = readRDS("thesis_output/mses_RK.rds")
nmaxes = c(10, 50, 100, 250, 500, 1000)
kriging = c(rep("OK", 10), rep("RK", 10))

cv_df = data.frame(rbind(mses_ok, mses_rk))
colnames(cv_df) = nmaxes
cv_df$kriging = kriging
cv_df = cv_df %>% pivot_longer(1:6, names_to = "NbrSize", values_to="MSE")
cv_df %>% group_by(kriging, NbrSize) %>% summarize(mn = mean(MSE), se = sd(MSE)) %>%
  ggplot(aes(x=NbrSize, y=mn)) + geom_point() + geom_errorbar(aes(ymin=mn-se, ymax=mn+se)) + facet_grid(~kriging)


get_ci <- function(preds){
  upper = preds@data$bt_preds + 2* sqrt(preds@data$var_bt)
  lower = preds@data$bt_preds - 2* sqrt(preds@data$var_bt)
  
  is_in_pi = (preds@data$value >= lower) & (preds@data$value <= upper)
  
  return(is_in_pi)
  
}


rk_fr06 = preds@data[preds@data$loc == "FR06003",]
ok_fr06 = preds_ok@data[preds_ok@data$loc == "FR06003",]
ols_fr06 = preds_ok@data[preds_ols@data$loc == "FR06003",]
rk_fr06 %>% ggplot() +geom_line(aes(x=date, y=bt_preds))
rk_fr06 %>% ggplot() +geom_line(aes(x=date, y=bt_preds, color = "RK")) + 
  geom_line(aes(x=date, y=value, color="True")) + 
  geom_line(data =ok_fr06, aes(x=date, y=bt_preds, color="OK")) +
  scale_color_manual(name="Type", values=c("True" = "black", "RK" = "red", "OK" = "blue"))



##############################################################################
# Anisotropy
##############################################################################

gammas <- matrix(data=NA, nrow = 366, ncol=96)
dates = seq(as.Date("2020-01-01"), as.Date("2020-12-15"), 1)
np = vector(length=96)
for (i in 1:366){
  temp = Train %>% filter(date == dates[i])
  temp = temp[!is.na(temp$value),]
  if (nrow(temp) == 0){ break }
  temp = SpatialPointsDataFrame(coords = temp[,2:3], data=temp[,-c(2:3)], proj4string = CRS("EPSG:4326"))
  temp = spTransform(temp, CRS("EPSG:2154"))
  e <- variogram(trans_value~1, data=temp, alpha=c(0,45,90,135, 180, 225, 270, 315), cutoff = 4.2e5, width = 4.2e5/12)
  if (length(e$gamma) == 96){
    gammas[i, ] = e$gamma * e$np
    np = np + e$np
    }

}

 e$gamma = colSums(gammas, na.rm=T) / np
 e$np = np
# saveRDS(e, file="thesis_output/data_anis.rds")

 
g1 <- e %>% ggplot(aes(x=dist, y=gamma, color=direction)) +geom_point() + geom_line() +theme_hc() +labs(x="Distance(m)", y="Semivariance", color="Direction")
g2 <- resids_e %>% ggplot(aes(x=dist, y=gamma, color=direction)) +geom_point() + geom_line() +theme_hc() +labs(x="Distance(m)", y="Semivariance", color="Direction")

preds_rk@data$var_bt = (exp(preds_rk@data$var) - 1) * exp(2*preds_rk@data$rk_pred+ preds_rk@data$var)


filt$is_train = filt$loc %in% train_loc$loc
blank_rast = raster("blankRaster.grd")
blank_rast = projectRaster(blank_rast, crs=sp::CRS(sf::st_crs(2154)[[2]]))

is_train = filt$loc %in% train_loc$loc
grouped = filt %>% group_by(loc, lonL, latL) %>% summarise(is_train = mean(is_train))
df = data.frame(rasterToPoints(blank_rast))
df %>% ggplot() + geom_tile(aes(x=x, y=y, color=value), fill="transparent") + theme_void() + 
  geom_jitter(data=grouped %>% filter(is_train == TRUE),aes(x=lonL*10000, y=latL*10000), color="forestgreen") + 
  geom_jitter(data=grouped %>% filter(is_train == FALSE),aes(x=lonL*10000, y=latL*10000), color="red") + 
  theme(legend.position = "none") + coord_sf(crs=st_crs(2154))



g1 = grid %>% filter(date == "2020-04-16") %>% ggplot(aes(x=lonL, y=latL, fill=aod47)) + geom_tile() + theme_void() + labs(fill="AOD47")

g2 = grid %>% filter(date == "2020-04-16") %>% ggplot(aes(x=lonL, y=latL, fill=aod47_imp)) + geom_tile() + theme_void() + labs(fill="AOD47")
grid.arrange(g1,g2)



stplot(grid_rk[,,"rk_preds"], main="Residual Kriging Predictions")
stplot(grid_rk[,,"se"], main="Residual Kriging Prediction Std Errors")


full_metrics = readRDS("thesis_outputs/full_cv_metrics_rk.RDS")
ngbrs = rep(c(10,50,100,250,500,1000), 5)
rmspe = c(colMeans(sqrt(full_metrics[[1]][[1]])), 
          colMeans(sqrt(full_metrics[[2]][[1]])),
          colMeans(sqrt(full_metrics[[3]][[1]])),
          colMeans(sqrt(full_metrics[[4]][[1]])),
          colMeans(sqrt(full_metrics[[5]][[1]])))
model = c(rep("Metric", 6),
          rep("ProdSum", 6),
          rep("Separable", 6),
          rep("SumMetric", 6),
          rep("sSumMetric", 6))
df = data.frame(cbind(ngbrs, rmspe, model))
df$ngbrs = as.numeric(df$ngbrs)
df$rmspe = as.numeric(df$rmspe)
g1 = df %>% ggplot(aes(x=ngbrs, y=rmspe, color=model)) + geom_line() + labs(x="Number of Neighbors", y="CV RMSPE", color="Model", title="Residual Local Kriging Cross-Validation") + theme_tufte()


met = readRDS("thesis_output/cv_ok_metric.rds")
ps = readRDS("thesis_output/cv_ok_ps.rds")
sep = readRDS("thesis_output/cv_ok_sep.rds")
smm = readRDS("thesis_output/cv_ok_smm.rds")
ssm = readRDS("thesis_output/cv_ok_ssm.rds")

df$rmspe = as.numeric(c(sqrt(colMeans(met)),
             sqrt(colMeans(ps)),
             sqrt(colMeans(sep)),
             sqrt(colMeans(smm)),
             sqrt(colMeans(ssm))))
g2 = df %>% ggplot(aes(x=ngbrs, y=rmspe, color=model)) + geom_line() + labs(x="Number of Neighbors", y="CV RMSPE", color="Model", title="Residual Local Kriging Cross-Validation") + theme_tufte()
