library(randomForest)
library(doParallel)
library(missForest)
library(tidyverse)
library(caret)
library(parallel)
filt = readRDS("pm_data.rds")
locations = filt %>% distinct(loc, .keep_all=TRUE) %>% dplyr::select(loc, lat, lon)

## Train-Test Split
set.seed(123) #formerly 111

test_samp = sample(1:nrow(locations), size=26, replace=FALSE)
Test = filt %>% filter(loc %in% locations$loc[test_samp])

train_loc = locations[-test_samp,]


filt$is_na = is.na(filt$aod47)

control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
tunegrid <- expand.grid(mtry=c(3,5,7,10))
modellist <- list()
dataset <-filt[!filt$is_na,c("lonL",
                                      "latL",
                                      "t",
                                      "aod47", 
                                      "CAMS_aod469_day",
                                      "CAMS_aod469_night",
                                      "CAMS_aod550_day",
                                      "CAMS_aod550_night",
                                      "CAMS_aod670_day",
                                      "CAMS_aod670_night",
                                      "CAMS_aod865_day",
                                      "CAMS_aod865_night",
                                      "CAMS_aod1240_day",
                                      "CAMS_aod1240_night")]
for (ntree in c(10, 50, 100)) {
  set.seed(111)
  fit <- train(aod47~., data=dataset, method="rf", tuneGrid=tunegrid, trControl=control, ntree=ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
  print(paste(ntree, " done"))
}
# compare results
results <- resamples(modellist)
summary(results)
dotplot(results)



saveRDS(results, "results.rds")


aod.rf <- randomForest(trans_value~ elevation + ndvi + aod47_imp+ humidity  +pressure + wind + temp + elevation + humidity + ndvi + trees_prop + crops_prop + precipitation + nearest_res + nearest_maj + hway_density + nearest_sec + water_prop + built_prop + t + pop_density + nearest_ter + res_density + month + weekday,
                       data=Train[!is.na(Train$value),],
                   ntree = 500,
                   niter = 3,
                   do.Trace = TRUE,
                   )




library("VSURF")
set.seed(414)

pm_surf <- VSURF(trans_value~ elevation + ndvi + aod47_imp+ humidity  +pressure + wind + temp + elevation + humidity + ndvi + trees_prop + crops_prop + precipitation + nearest_res + nearest_maj + hway_density + nearest_sec + water_prop + built_prop + t + pop_density + nearest_ter + res_density,
                 parallel=TRUE, verbose=TRUE, data=Train[!is.na(Train$value),], ntree=300, ncores = (detectCores()-1))
print(Sys.time())

saveRDS(pm_surf, "pm_surf")
