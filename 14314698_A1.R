################ GEOG71922 Assessment 1
################ Meles meles species distribution modelling

################ libraries

library(terra)
library(sf)
library(dplyr)
library(ggplot2)

################ set working directory

setwd("C:/Users/lenovo/Downloads/GEOG71922_Assessment_1")

################ read in data

meles = read.csv("Melesmeles.csv")
scot = st_read("scotSamp.shp")
LCM = rast("LCMUK.tif")
demScot = rast("demScotland.tif")

################ inspect data

head(meles)
names(meles)
dim(meles)

print(scot)
LCM
demScot

################ check verification field

table(meles$Identification.verification.status, useNA = "ifany")

################ remove unconfirmed records

meles = meles[meles$Identification.verification.status != "Unconfirmed", ]

################ remove records with coordinate uncertainty > 1000 m

summary(meles$Coordinate.uncertainty_m)

meles = meles[meles$Coordinate.uncertainty_m <= 1000 | 
                is.na(meles$Coordinate.uncertainty_m), ]

################ remove records with missing coordinates

meles = meles[!is.na(meles$Longitude), ]
meles = meles[!is.na(meles$Latitude), ]

################ check cleaned data

dim(meles)
head(meles)

################ convert to sf points

meles_sf = st_as_sf(meles,
                    coords = c("Longitude", "Latitude"),
                    crs = 4326)

################ project to British National Grid

meles_sf = st_transform(meles_sf, 27700)

################ clip to study area

meles_sf = meles_sf[scot, ]

################ inspect result

plot(st_geometry(scot))
plot(st_geometry(meles_sf), add = TRUE, col = "red", pch = 16, cex = 0.3)

print(meles_sf)

################ create broadleaf raster

broadleaf = LCM == 1
plot(broadleaf)

################ create presence variable

meles_sf$Pres = 1

################ create pseudo absence points

set.seed(1)

absence = st_sample(scot,
                    size = nrow(meles_sf),
                    type = "random")

absence = st_sf(Pres = 0,
                geometry = absence)

################ combine presence and absence

pres_sf = meles_sf[, "Pres"]

all_points = rbind(pres_sf, absence)

################ check result

plot(st_geometry(scot))
plot(st_geometry(all_points),
     add = TRUE,
     col = c("blue","red")[all_points$Pres + 1],
     pch = 16)

################ buffer sizes

radii = seq(100,2000,100)

radii

################ characteristic scale analysis

best_scale = 700
best_scale

################ extract covariates using best scale

buf_best = st_buffer(all_points, dist = best_scale)

broad_best = terra::extract(broadleaf, vect(buf_best),
                            fun = mean,
                            na.rm = TRUE)

all_points$broad = broad_best[,2]

################ extract covariates using best scale

buf_best = st_buffer(all_points, dist = best_scale)

broad_best = terra::extract(broadleaf, vect(buf_best),
                            fun = mean,
                            na.rm = TRUE)

all_points$broad = broad_best[,2]

################ extract elevation

elev = terra::extract(demScot, vect(all_points))

all_points$elev = elev[,2]

################ create modelling dataframe

coords = st_coordinates(all_points)

all.cov = st_drop_geometry(all_points)

all.cov$x = coords[,1]
all.cov$y = coords[,2]

all.cov = na.omit(all.cov)

################ check data

dim(all.cov)
head(all.cov)

################ GLM model

glm1 = glm(Pres ~ broad + elev,
           data = all.cov,
           family = "binomial")

summary(glm1)

################ predicted probability

all.cov$glm_pred = predict(glm1,
                           type = "response")

################ broadleaf response curve

new.broad = data.frame(
  broad = seq(min(all.cov$broad),
              max(all.cov$broad),
              length.out = 100),
  elev = mean(all.cov$elev)
)

new.broad$pred = predict(glm1,
                         newdata = new.broad,
                         type = "response")

plot(new.broad$broad,
     new.broad$pred,
     type = "l",
     lwd = 2,
     xlab = "Broadleaf proportion",
     ylab = "Predicted probability")

new.elev = data.frame(
  broad = mean(all.cov$broad),
  elev = seq(min(all.cov$elev),
             max(all.cov$elev),
             length.out = 100)
)

new.elev$pred = predict(glm1,
                        newdata = new.elev,
                        type = "response")

plot(new.elev$elev,
     new.elev$pred,
     type = "l",
     lwd = 2,
     xlab = "Elevation",
     ylab = "Predicted probability")

library(randomForest)

rf1 = randomForest(
  factor(Pres) ~ broad + elev,
  data = all.cov,
  ntree = 500
)

rf1

all.cov$rf_pred = predict(
  rf1,
  type = "prob"
)[,2]

varImpPlot(rf1)

install.packages("pROC")
library(pROC)

################ ROC and AUC

roc_glm = roc(all.cov$Pres, all.cov$glm_pred)
roc_rf = roc(all.cov$Pres, all.cov$rf_pred)

auc(roc_glm)
auc(roc_rf)

plot(roc_glm, col = "black", lwd = 2)
plot(roc_rf, add = TRUE, col = "red", lwd = 2)

legend("bottomright",
       legend = c("GLM", "RF"),
       col = c("black", "red"),
       lwd = 2)

broad_raster = focal(
  LCM == 1,
  w = matrix(1, 7, 7),
  fun = mean,
  na.rm = TRUE
)

elev_raster = demScot

nlyr(broad_raster)
nlyr(demScot)

broad1 = LCM[[1]] == 1

broad_raster = focal(
  broad1,
  w = matrix(1, 7, 7),
  fun = mean,
  na.rm = TRUE
)

broad_raster = terra::resample(broad_raster, demScot)

predictors = c(broad_raster, demScot)
names(predictors) = c("broad", "elev")

glm_map = terra::predict(
  predictors,
  glm1,
  type = "response"
)

plot(glm_map)

rf_map = terra::predict(
  predictors,
  rf1,
  type = "prob",
  index = 2
)
plot(rf_map)

################ spatial cross-validation
install.packages("cowplot")
library(mlr)
library(cowplot)
library(ranger)

################ prepare task for mlr

task = all.cov[, c("broad", "elev", "Pres", "x", "y")]
task$Pres = as.factor(task$Pres)

task = makeClassifTask(
  data = task[, c("broad", "elev", "Pres")],
  target = "Pres",
  positive = "1",
  coordinates = task[, c("x", "y")]
)

################ learners

lrnBinomial = makeLearner(
  "classif.binomial",
  predict.type = "prob",
  fix.factors.prediction = TRUE
)

lrnRF = makeLearner(
  "classif.ranger",
  predict.type = "prob",
  fix.factors.prediction = TRUE
)

################ resampling schemes

perf_levelCV = makeResampleDesc(
  method = "RepCV",
  predict = "test",
  folds = 5,
  reps = 5
)

perf_level_spCV = makeResampleDesc(
  method = "SpRepCV",
  folds = 5,
  reps = 5
)

################ GLM conventional CV

cvBinomial = resample(
  learner = lrnBinomial,
  task = task,
  resampling = perf_levelCV,
  measures = mlr::auc,
  show.info = FALSE
)

print(cvBinomial)

################ GLM spatial CV

sp_cvBinomial = resample(
  learner = lrnBinomial,
  task = task,
  resampling = perf_level_spCV,
  measures = mlr::auc,
  show.info = FALSE
)

print(sp_cvBinomial)

################ RF conventional CV

cvRF = resample(
  learner = lrnRF,
  task = task,
  resampling = perf_levelCV,
  measures = mlr::auc,
  show.info = FALSE
)

print(cvRF)

################ RF spatial CV

sp_cvRF = resample(
  learner = lrnRF,
  task = task,
  resampling = perf_level_spCV,
  measures = mlr::auc,
  show.info = FALSE
)

print(sp_cvRF)

auc_values <- c(0.783, 0.767, 0.799, 0.736)

models <- c(
  "GLM random",
  "GLM spatial",
  "RF random",
  "RF spatial"
)

barplot(
  auc_values,
  names.arg = models,
  col = c("steelblue","steelblue","orange","orange"),
  ylab = "AUC",
  main = "Model performance comparison",
  ylim = c(0,1)
)



plotsSP_GLM = createSpatialResamplingPlots(
  task,
  resample = sp_cvBinomial,
  crs = terra::crs(demScot),
  datum = terra::crs(demScot),
  color.test = "red",
  point.size = 1
)

png("GLM_spatial_CV_folds.png", width = 1800, height = 2400)

cowplot::plot_grid(
  plotlist = plotsSP_GLM[["Plots"]],
  ncol = 1,
  labels = plotsSP_GLM[["Labels"]]
)

dev.off()


plotsSP_RF = createSpatialResamplingPlots(
  task,
  resample = sp_cvRF,
  crs = terra::crs(demScot),
  datum = terra::crs(demScot),
  color.test = "red",
  point.size = 1
)

png("RF_spatial_CV_folds.png", width = 1800, height = 2400)

cowplot::plot_grid(
  plotlist = plotsSP_RF[["Plots"]],
  ncol = 1,
  labels = plotsSP_RF[["Labels"]]
)

dev.off()
