
rk_preds <- function(trn, val, covs, formula, nmax, model){
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
  
  stAni = covs$stAni
  if(is.null(stAni)){
    stAni = 245225
  }
  
  clusterExport(cl = cl, varlist = c("stAni", "trn", "val", "covs", "parts", "nmax"), envir = .GlobalEnv)
  clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
  
  parallelX <- parLapply(cl = cl, X = 1:no_cores, fun = function(x) krigeST(resids~1 , # latitude trend
                                                                            data = trn, # data set w/o 14 July
                                                                            newdata = val[, parts[[x]]], # prediction grid
                                                                            modelList = covs, # semivariogram
                                                                            nmax = nmax,
                                                                            stAni = stAni,
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
  
  
  ST_pred@data$rk_preds = ST_pred@data$preds + ST_pred@data$var1.pred
  ST_pred@data$var = ST_pred@data$var1.var + lm_var
  
  ST_pred@data$bt_preds = exp(ST_pred@data$rk_preds + ST_pred@data$var/2) - 1
  ST_pred@data$bt_var = (exp(ST_pred@data$var) - 1) * exp(2*ST_pred@data$rk_preds + ST_pred@data$var/2)
  
  errors[[1]] = mses
  errors[[2]] = maes
  errors[[3]] = sds
  errors[[4]] = mses_trans
  errors[[5]] = ST_pred
  print(paste0("Predictions done ", Sys.time()))
  
  return(errors)
}



get_pred_grid_ok <- function(trn, val, covs, formula, nmax){
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
  
  parallelX <- parLapply(cl = cl, X = 1:no_cores, fun = function(x) krigeST(trans_value~1, # latitude trend
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
  ST_pred@data$bt_preds = exp(ST_pred@data$var1.pred + ST_pred@data$var1.var/2)
  ST_pred@data$bt_var = (exp(ST_pred@data$var1.var) - 1) * exp(2*ST_pred@data$var1.pred + ST_pred@data$var1.var/2)
  
  print(paste0("Predictions done ", Sys.time()))
  
  return(ST_pred)
}


all_variograms = function(vv){
  ret = list()
  sph_m = vgm(psill = 0.05, "Sph", range = 3e5, nugget = 0)
  exp_m = vgm(psill = 0.05, "Exp", range = 3e5, nugget = 0)
  
  sph_t = vgm(psill = 0.05, "Sph", range = 5, nugget = 0)
  exp_t = vgm(psill = 0.05, "Exp", range = 5, nugget = 0)
  
  sph_met = vgm(psill = 0.35, "Sph", range = 1e5, nugget = 0)
  exp_met = vgm(psill = 0.35, "Exp", range = 1e5, nugget = 0)
  
  sph_joint = vgm(psill = 0.2, "Sph", range = 4e5, nugget = 0)
  exp_joint = vgm(psill = 0.2, "Exp", range = 4e5, nugget = 0)
  
  
  metrics = list(get_metric(vv, exp_joint), get_metric(vv, sph_joint))
  print("metrics done")
  ps = list(get_ps(vv, exp_m, exp_t), get_ps(vv, exp_m, sph_t), get_ps(vv, sph_m, exp_t), get_ps(vv, sph_m, sph_t))
  print("ps done")
  sep = list(get_sep(vv, exp_m, exp_t), get_sep(vv, exp_m, sph_t), get_sep(vv, sph_m, exp_t), get_sep(vv, sph_m, sph_t))
  print("sep done")
  smm = list(get_sMM(vv, exp_m, exp_t, exp_joint), get_sMM(vv, exp_m, sph_t,exp_joint), get_sMM(vv, sph_m, exp_t,exp_joint), get_sMM(vv, sph_m, sph_t,exp_joint),
             get_sMM(vv, exp_m, exp_t, sph_joint), get_sMM(vv, exp_m, sph_t,sph_joint), get_sMM(vv, sph_m, exp_t,sph_joint), get_sMM(vv, sph_m, sph_t,sph_joint))
  print("smm done")
  ssm = list(get_sSM(vv, exp_m, exp_t, exp_joint), get_sSM(vv, exp_m, sph_t,exp_joint), get_sSM(vv, sph_m, exp_t,exp_joint), get_sSM(vv, sph_m, sph_t,exp_joint),
             get_sSM(vv, exp_m, exp_t, sph_joint), get_sSM(vv, exp_m, sph_t,sph_joint), get_sSM(vv, sph_m, exp_t,sph_joint), get_sSM(vv, sph_m, sph_t,sph_joint))
  print("ssm done")
  
  return(list(metrics, ps, sep, smm, ssm))
}


get_metric <- function(vv, joint){
  

  jE <- vgmST(stModel = "metric",
              joint = joint,
              stAni = 100000)
  fit_a <- fit.StVariogram(vv, jE,lower=c(1e-10,1e-10,1e-10,1e-10),
                           control = list(parscale = c(0.1, 500, 0.01, 1000)),
                           upper = c(sill = 0.5, range = 8e5, nugget = 0.1, anis = 500000),
                           method = "L-BFGS-B",
                           fit.method = 2)

  
  return(fit_a)
}

get_ps <- function(vv, space, time){


  prodSum <- vgmST("productSum", space = space,
                   time = time,
                   k=2)
  pS = fit.StVariogram(vv, prodSum, lower=c(sill.s = 0.01,range.s = 5000, nugget.s = 0,
                                            sill.t = 0.01, range.t = 0.5, nugget.t = 0, k.k = 0.01),
                       upper=c(sill.s = 0.5,range.s = 8e5, nugget.s = 0.25, 
                               sill.t = 0.5, range.t = 15, nugget.t = 0.25, k.k = 30),
                       method = "L-BFGS-B",
                       fit.method = 2)

  return(pS)
}


get_sep <- function(vv, space, time){
  
  separable <- vgmST(stModel = "separable",
                     space = space,
                     time = time,
                     sill = 0.35)
  separable <- fit.StVariogram(vv, separable, lower=c(range.s = 5000, nugget.s = 0, range.t = 1,
                                                      nugget.t = 0, sill = 0.25),
                               upper=c(range.s = 8e5, nugget.s = 0.2, range.t = 15,
                                       nugget.t = 0.2, sill = 0.5),
                               control = list(parscale = c(500, 0.1, 1 , 0.1, 10)),
                               method = "L-BFGS-B",
                               fit.method = 2, stAni= 1.5)
  
  return(separable)
}




get_sMM <- function(vv, space, time, joint){
  
  sMM <- vgmST("sumMetric",
               space = space,
               time = time,
               joint = joint,
               stAni = 100000)
  
  sMM <- fit.StVariogram(vv, sMM, fit.method = 2,
                         method = "L-BFGS-B",
                         control = list(parscale = c(.1, 1000, .01, .1, 0.5, .01,
                                                     .1, 1000, .1, 10000),
                                        maxit=10000),
                         lower = c(sill.s = 0, range.s = 10, nugget.s = 0,
                                   sill.t = 0, range.t = 0.1, nugget.t = 0,
                                   sill.st = 0, range.st = 10, nugget.st = 0,
                                   anis = 40),
                         upper = c(sill.s = 0.5, range.s = 8e5, nugget.s = 0.25,
                                   sill.t = 0.5, range.t = 15, nugget.t = 0.25,
                                   sill.st = 0.5, range.st = 10e5, nugget.st = 0.25,
                                   anis = 800000))
  print(Sys.time())
  
  return(sMM)
  
}

get_sSM <- function(vv, space, time, joint){
  
  sSM <- vgmST("simpleSumMetric",
               space = space,
               time = time,
               joint = joint,
               nugget = 0,
               stAni = 100000)
  
  sSM <- fit.StVariogram(vv, sSM, fit.method = 2,
                         method = "L-BFGS-B",
                         control = list(parscale = c(.1, 1000, .1, 0.5,
                                                     .1, 1000, .01, 1000),
                                        maxit=10000),
                         lower = c(sill.s = 0, range.s = 10,
                                   sill.t = 0, range.t = 0.1,
                                   sill.st = 0, range.st = 10, nugget = 0,
                                   anis = 40),
                         upper = c(sill.s = 0.5, range.s = 8e5,
                                   sill.t = 0.5, range.t = 10,
                                   sill.st = 0.5, range.st = 10e5, nugget = 0.25,
                                   anis = 500000))
  
  return(sSM)
  
}






library(caret)
library(mlbench)
library(Hmisc)
library(doMc)

ctrl <- rfeControl(functions = lmFuncs,
                   rerank=TRUE,
                   method = "repeatedcv",
                   repeats = 3,
                   verbose = TRUE)

rfe(x=Train[,-c(1:5,29)], y=Train$trans_value, subsets <- c(1:5, 10, 15, 20, 25), na.action=na.omit)


library(glinternet)


trn = Train[!is.na(Train$trans_value), -c(1:5,12:13)]
y = trn$trans_value
trn = trn %>% select(-trans_value)
trn = trn %>% mutate_if(is.numeric, scale)
i_num <- sapply(trn, is.numeric)
numLevels <- trn %>% sapply(nlevels)
numLevels[23] = 7
numLevels[numLevels==0] <- 1
trn[, !i_num] <- apply(trn[, !i_num], 2, function(col) as.integer(as.factor(col)) - 1)


set.seed(1001)
cv_fit <- glinternet.cv(trn, y, numLevels, verbose=TRUE)
plot(cv_fit)

num_covs = vector()
mains = vector()
ints = vector()
for (i in 1:50){
  coefs <- coef(cv_fit$glinternetFit)[[i]]
  
  m_effects = length(coefs$mainEffects[[1]]) + length(coefs$mainEffects[[2]])
  int_effects = length(coefs$interactions[[1]]) + length(coefs$interactions[[2]])  +length(coefs$interactions[[3]]) 
  
  mains[i] = m_effects
  ints[i] = int_effects
  num_covs[i] = m_effects + int_effects
}

cv_fit <- readRDS("thesis_output/olsInteractionCV.rds")

idx_num <- (1:length(i_num))[i_num]
idx_cat <- (1:length(i_num))[!i_num]
names(numLevels)[idx_cat[coefs$mainEffects$cat]]

basic_fit <- lm(trans_value~elevation + ndvi + temp+pressure+wind+aod47_imp+precipitation+built_prop+pop_density+trees_prop+ ndvi:trees_prop + temp:pressure + temp:wind + temp:precipitation, data=Train)
