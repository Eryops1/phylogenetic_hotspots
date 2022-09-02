# old school spatial models 


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
library(raster)
library(sf)
library(ggplot2)
theme_set(theme_bw()+theme(text=element_text(size=9), panel.grid=element_blank()))
library(cowplot)
#library(rgdal)
library(mgcv)
library(caret)
library(gbm)
library(parallel)
library(spatialRF)
library(scico)
if(!dir.exists("figures"))dir.create("figures")

# DATA --------------------------------------------------------------------

shp <- readRDS("fin_shp.rds")

# remove not needed data
#shp <- cbind(PD_obs=shp$PD_obs_ts, shp)
shp <- shp[,-grep("obs_p|obs_rank|reps|_cw$|LEVEL2|LEVEL1|LEVEL_3_CO|LEVEL_NAME|ID|\\.3|_rw|CONTI|REGION", names(shp))]
names(shp)<- gsub("\\.1", "_mean", names(shp))
names(shp)<- gsub("\\.2", "_sd", names(shp))
# remove BOU that has no data
shp <- shp[!shp$LEVEL3_COD=="BOU",]

# tmp <- shp[,-grep("rand_|hfp|deforest|bio|change", names(shp))]
# # tmp$centroids <- st_centroid(tmp) %>% 
# #   st_coordinates()
# # tmp$y <- tmp$centroids[,1] # x=lng
# # tmp$x <- tmp$centroids[,2]
# #tmp <- st_drop_geometry(tmp)
# tmp <- tmp[,-grep("centroids|PD_obs|PE_sd|PE|WE|mdr", names(tmp))]

names(shp)
dat <- shp[,c(1,6,13,16:44)]
names(dat)
dat <- na.omit(dat)
dim(dat)


# SpatialRF ---------------------------------------------------------------
## following https://blasbenito.github.io/spatialRF/
library(spatialRF)

#coordinates of the cases
dat$centroids <- st_centroid(dat) %>% 
  st_coordinates()
dat$y <- dat$centroids[,1] # x=lng
dat$x <- dat$centroids[,2]
xy <- st_drop_geometry(dat[, c("x", "y")])

#distance matrix
dm <- dist(xy, method = "euclidean", diag = TRUE, upper = TRUE)
dm <- as.matrix(dm)

#distance thresholds (same units as distance_matrix)
distance.thresholds <- 100000*c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)

## SES.PD ----------------------------
#random seed for reproducibility
random.seed <- 1

dat.ns <- st_drop_geometry(dat)
predictor.variable.names <- names(dat.ns)
predictor.variable.names <- predictor.variable.names[-c(1,2,33,34,35)]
dependent.variable.name <- "SES.PD"

spatialRF::plot_training_df(
  data = dat.ns,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  ncol = 5,
  point.color = viridis::viridis(100, option = "F"),
  line.color = "gray30"
)

## assess autocorrelation
spatialRF::plot_training_df_moran(
  data = dat.ns,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = dm,
  distance.thresholds = distance.thresholds,
  fill.color = viridis::viridis(
    100,
    option = "F",
    direction = -1
  ),
  point.color = "gray40"
)
ggsave("figures/morans_I_in_variables.png", width=7, height=5, dpi=300)


# Finding promising variable interactions
interactions <- spatialRF::the_feature_engineer(
  data = dat.ns,
  dependent.variable.name = "SES.PD",
  predictor.variable.names = predictor.variable.names,
  xy = xy,
  importance.threshold = 0.50, #uses 50% best predictors
  cor.threshold = 0.60, #max corr between interactions and predictors
  seed = random.seed,
  repetitions = 100,
  verbose = TRUE
)
# 3 interactions with area + biome (medit fws, mat_lgm_ano, desert_x_shrub)
## not really intuitive why they should be interacting, leave them for now


## non spatial RF
model.non.spatial <- spatialRF::rf(
  data = dat.ns,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = dm,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  verbose = FALSE
)

## residuals non spatial:
spatialRF::plot_residuals_diagnostics(
  model.non.spatial,
  verbose = FALSE
)

## Variable importance
spatialRF::plot_importance(
  model.non.spatial,
  verbose = FALSE
)
importance.df <- randomForestExplainer::measure_importance(
  model.non.spatial,
  measures = c("mean_min_depth", "no_of_nodes", "times_a_root", "p_value")
)
model.non.spatial <- spatialRF::rf_importance(
  model = model.non.spatial
)

model.non.spatial$importance$per.variable %>%
  ggplot2::ggplot() +
  ggplot2::aes(
    x = importance.oob,
    y = importance.cv
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::theme_bw() +
  ggplot2::xlab("Importance (out-of-bag)") +
  ggplot2::ylab("Contribution to transferability") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "red4")

## Local var importance
local.importance <- spatialRF::get_importance_local(model.non.spatial)

local.importance <- cbind(
  xy,
  local.importance
)
#colors
color.low <- viridis::viridis(
  3,
  option = "F"
)[2]
color.high <- viridis::viridis(
  3,
  option = "F"
)[1]

# #plot of climate_bio1_average
# p1 <- ggplot2::ggplot() +
#   ggplot2::geom_sf(
#     data = shp,
#     fill = "white"
#   ) +
#   ggplot2::geom_point(
#     data = local.importance,
#     ggplot2::aes(
#       x = y,
#       y = x,
#       color = mrd
#     )
#   ) +
#   ggplot2::scale_color_gradient2(
#     low = color.low,
#     high = color.high
#   ) +
#   ggplot2::theme_bw() +
#   ggplot2::theme(legend.position = "bottom") +
#   ggplot2::ggtitle("mrd") +
#   ggplot2::theme(
#     plot.title = ggplot2::element_text(hjust = 0.5),
#     legend.key.width = ggplot2::unit(1,"cm")
#   ) +
#   ggplot2::labs(color = "Importance") +
#   ggplot2::xlab("Longitude") +
#   ggplot2::ylab("Latitude")
# 
# p2 <- ggplot2::ggplot() +
#   ggplot2::geom_sf(
#     data = shp,
#     fill = "white"
#   ) +
#   ggplot2::geom_point(
#     data = local.importance,
#     ggplot2::aes(
#       x = y,
#       y = x,
#       color = soil
#     )
#   ) +
#   ggplot2::scale_color_gradient2(
#     low = color.low,
#     high = color.high
#   ) +
#   ggplot2::theme_bw() +
#   ggplot2::theme(legend.position = "bottom") +
#   ggplot2::ggtitle("soil") +
#   ggplot2::theme(
#     plot.title = ggplot2::element_text(hjust = 0.5),
#     legend.key.width = ggplot2::unit(1,"cm")
#   ) +
#   ggplot2::labs(color = "Importance") +
#   ggplot2::xlab("Longitude") +
#   ggplot2::ylab("Latitude")
# 
# p1 + p2
# # In these maps, values lower than 0 indicate that for a given record, the permuted version of the variable led to an accuracy score even higher than the one of the non-permuted variable, so again these negative values can be interpreted as “worse than chance”.
# spatialRF::plot_response_curves(
#   model.non.spatial,
#   quantiles = 0.5,
#   ncol = 3
# )
local.importance

# Model performance
spatialRF::print_performance(model.non.spatial)

# Response curves
spatialRF::plot_response_curves(
  model.non.spatial,
  quantiles = 0.5,
  ncol = 3
)

# Spatial cross-validation
# The function rf_evaluate() overcomes the limitations of the performance scores
# explained above by providing honest performance based on spatial
# cross-validation. The function separates the data into a number of spatially
# independent training and testing folds. Then, it fits a model on each training
# fold, predicts over each testing fold, and computes statistics of performance
# measures across folds. Let’s see how it works.
model.non.spatial <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy,                  #data coordinates
  repetitions = 30,         #number of spatial folds
  training.fraction = 0.75, #training data fraction on each fold
  metrics = "r.squared",
  seed = random.seed,
  verbose = FALSE
)

# spatial folds
spatialRF::plot_evaluation(model.non.spatial)
spatialRF::print_evaluation(model.non.spatial)




## The spatial RF!
spatialRF::plot_moran(
  model.non.spatial,
  verbose = FALSE
)

model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed
)
#The plot below shows the Moran’s I of the residuals of the spatial model, and
#indicates that the residuals are not autocorrelated at any distance.

spatialRF::plot_moran(
  model.spatial,
  verbose = FALSE
)
# much better!

p1 <- spatialRF::plot_importance(
  model.non.spatial,
  verbose = FALSE) +
  ggplot2::ggtitle("Non-spatial model")

p2 <- spatialRF::plot_importance(
  model.spatial,
  verbose = FALSE) +
  ggplot2::ggtitle("Spatial model")

p1 | p2

kableExtra::kbl(
  head(model.spatial$importance$per.variable, n = 15),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)
# spatial_predictor_X_Y, where X is the neighborhood distance at which the
# predictor has been generated, and Y is the index of the predictor.
# so I got one named 2000000_6, i.e. distance threshold 2000000 with variable 6: pre_mean

# look at spatial predictors
spatial.predictors <- spatialRF::get_spatial_predictors(model.spatial)
pr <- data.frame(spatial.predictors, dat.ns[,c("x", "y")])

p1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = dat, fill = "white") +
  ggplot2::geom_point(
    data = pr,
    ggplot2::aes(
      x = y,
      y = x,
      color = spatial_predictor_2000000_6
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "Eigenvalue") +
  ggplot2::ggtitle("Variable: spatial_predictor_2000000_6") +
  ggplot2::theme(legend.position = "bottom")+
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude")

p2 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = dat, fill = "white") +
  ggplot2::geom_point(
    data = pr,
    ggplot2::aes(
      x = y,
      y = x,
      color = spatial_predictor_1000000_11,
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "Eigenvalue") +
  ggplot2::ggtitle("Variable: spatial_predictor_1000000_11") +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("")

p1 | p2



p <- spatialRF::plot_optimization(model.spatial)

model.spatial.repeat <- spatialRF::rf_repeat(
  model = model.spatial,
  repetitions = 30,
  seed = random.seed,
  verbose = FALSE
)
#The importance scores of a model fitted with rf_repeat() are plotted as a
#violin plot, with the distribution of the importance scores of each predictor
#across repetitions.

spatialRF::plot_importance(
  model.spatial.repeat,
  verbose = FALSE
)

spatialRF::plot_response_curves(
  model.spatial.repeat,
  quantiles = 0.5,
  ncol = 3
)

spatialRF::print_performance(model.spatial.repeat)

# Full model
# model.full <- rf_spatial(
#   data = st_drop_geometry(shp2),
#   dependent.variable.name = dependent.variable.name,
#   predictor.variable.names = predictor.variable.names,
#   distance.matrix = dm,
#   distance.thresholds = distance.thresholds,
#   xy = xy
# ) %>%
#   rf_tuning() %>%
#   rf_evaluate() %>%
#   rf_repeat()


comparison <- spatialRF::rf_compare(
  models = list(
    `Non-spatial` = model.non.spatial,
    `Spatial` = model.spatial
  ),
  xy = xy,
  repetitions = 30,
  training.fraction = 0.8,
  metrics = "r.squared",
  seed = random.seed
)

x <- comparison$comparison.df %>%
  dplyr::group_by(model, metric) %>%
  dplyr::summarise(value = round(median(value), 3)) %>%
  dplyr::arrange(metric) %>%
  as.data.frame()
colnames(x) <- c("Model", "Metric", "Median")
kableExtra::kbl(
  x,
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)






# Spatial GBM -------------------

shp <- readRDS("fin_shp.rds")

# remove not needed data
shp <- shp[,-grep("obs_p|obs_rank|reps|_cw$|LEVEL2|LEVEL1|LEVEL_3_CO|LEVEL_NAME|ID|\\.3|_rw|CONTI|REGION", names(shp))]
names(shp)<- gsub("\\.1", "_mean", names(shp))
names(shp)<- gsub("\\.2", "_sd", names(shp))
dat <- shp[,c(1,6,13,16:44)]
names(dat)
dat <- na.omit(dat)
dim(dat)


#coordinates of the cases
dat$centroids <- st_centroid(dat) %>% 
  st_coordinates()
dat$y <- dat$centroids[,1] # x=lng
dat$x <- dat$centroids[,2]
xy <- st_drop_geometry(dat[, c("x", "y")])
dat.ns <- st_drop_geometry(dat)

#distance matrix
dm <- dist(xy, method = "euclidean", diag = TRUE, upper = TRUE)
dm <- as.matrix(dm)

#distance thresholds (same units as distance_matrix)
distance.thresholds <- 100000*c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)

## get spatial predictors ----
mems <- spatialRF::mem_multithreshold(
  distance.matrix = dm,
  distance.thresholds = distance.thresholds
)
kableExtra::kbl(
  head(mems[, 1:4], n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

mem.rank <- spatialRF::rank_spatial_predictors(
  distance.matrix = dm,
  spatial.predictors.df = mems,
  ranking.method = "moran"
)
mems <- mems[, mem.rank$ranking]
dat.ns <- cbind(dat.ns, mems[,1:10])


## tuning parameter with caret package ----
gbmControl <- trainControl(method = "repeatedcv", number = 10,
                           repeats = 3, savePredictions = "final",
                           returnResamp ="final")# allowParallel=F has become necessary recently, not sure why.
gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3),
                        n.trees = (1:10)*50,
                        shrinkage = c(0.1, 0.01, 0.001),
                        n.minobsinnode = c(5, 10, 15))
# define gbm function
run.gbm <- function(i, seeds, data, trControl=gbmControl,
                    tuneGrid=gbmGrid, f){
  set.seed(s[i])
  if(!i%%1)cat(i,"\r")
  temp <- train(f, data, method = "gbm",
                trControl = gbmControl, verbose=FALSE,
                tuneGrid = gbmGrid)
  return(temp)
}

## define gbm formula. choose here spatial predictors if desired ----
n_spatial <- 2

n <- names(dat.ns[,-grep("LEVEL|centroids|x|y", names(dat.ns))])
spat <- grep("spatial", n) 
non_spat <- seq(n)[-spat]
n <- n[c(non_spat, spat[seq(n_spatial)])]
f.pd <- as.formula(paste("SES.PD ~", paste(n[!n %in% c("SES.PD", "SES.PE")], collapse = " + ")))
f.pe <- as.formula(paste("SES.PE ~", paste(n[!n %in% c("SES.PD", "SES.PE")], collapse = " + ")))

## RUN GBMs (sequential run) ----
foreach::registerDoSEQ()
s <- seq(500,509,1)
res <- list()
for(i in 1:10){
  res[[i]] <- run.gbm(i=i, seeds=s[i], data=dat.ns, f=f.pd)
  beepr::beep(2)
}
saveRDS(res, paste0("idata/res10GBM_PD", n_spatial, "spat_predictors.rds"))

foreach::registerDoSEQ()
s <- seq(500,509,1)
res <- list()
for(i in 1:10){
  res[[i]] <- run.gbm(i=i, seeds=s[i], data=dat.ns, f=f.pe)
  beepr::beep(2)
}
saveRDS(res, paste0("idata/res10GBM_PE", n_spatial, "spat_predictors.rds"))


# RESULTS ----
# 
# ## analyse PD-----
# shp <- readRDS("fin_shp.rds")
# shp <- shp[,-grep("obs_p|obs_rank|reps|_cw$|LEVEL2|LEVEL1|LEVEL_3_CO|LEVEL_NAME|ID|\\.3|_rw|CONTI|REGION", names(shp))]
# names(shp)<- gsub("\\.1", "_mean", names(shp))
# names(shp)<- gsub("\\.2", "_sd", names(shp))
# dat <- shp[,c(1,6,13,16:44)]
# dat <- na.omit(dat)
# dat.ns <- st_drop_geometry(dat)
# 
# #coordinates of the cases
# dat$centroids <- st_centroid(dat) %>% 
#   st_coordinates()
# dat$y <- dat$centroids[,1] # x=lng
# dat$x <- dat$centroids[,2]
# xy <- st_drop_geometry(dat[, c("x", "y")])
# 
# n_spatial <- 2
# res <- readRDS(paste0("idata/res10GBM_PD", n_spatial, "spat_predictors.rds"))
# 
# ## morans I for residuals ----
# # get model averages
# resids_list <- lapply(res, residuals)
# resids_df <- as.data.frame(resids_list, col.names=c(1:length(res)))
# dat.ns$PD_residual_median <- apply(resids_df, 1, median)
# resids_df$SES.PD <- dat.ns$SES.PD
# 
# resids_df_l <- tidyr::pivot_longer(resids_df, cols=1:length(res))
# ggplot(resids_df_l, aes(x=SES.PD, y=value, group=SES.PD))+
#   geom_boxplot()+
#   theme(axis.text.x=element_text(angle=45))+
#   geom_hline(yintercept=0)+
#   ylab(paste0("SES.PD GBM model residuals ", n_spatial, "spat pred, 10 runs"))
# 
# # distance matrix
# dm <- dist(xy, method = "euclidean", diag = TRUE, upper = TRUE)
# dm <- as.matrix(dm)
# 
# #distance thresholds (same units as distance_matrix = METERS)
# distance.thresholds <- 100000*c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)
# 
# # distance bands
# moran <- data.frame(dist.class = distance.thresholds,
#                     moransI = NA,
#                     moransp = NA)
# 
# coo <- cbind(dat$y, dat$x)
# for(i in 1:length(moran$dist.class)){
#   S.dist  <-  spdep::dnearneigh(coo, 0, moran$dist.class[i], longlat = FALSE)
#   lw <- spdep::nb2listw(S.dist, style="W",zero.policy=T) 
#   
#   MI <- spdep::moran.mc(dat.ns$PD_residual_median, lw, nsim=999,zero.policy=T) 
#   moran$moransI[i] <- MI$statistic
#   moran$moransp[i] <- MI$p.value
#   if(!i%%1)cat(i,"\r")
# }
# temp <- tidyr::pivot_longer(moran, cols = c("moransI"))
# ggplot(temp, aes(x=dist.class, y=value)) +
#   geom_line()+
#   geom_point(aes(size=moransp<0.05))+
#   scale_x_continuous("Distance class (m)")+
#   scale_color_discrete("GBM residuals")+
#   ylab("Moran's I for SES.PD GBM residuals")+
#   theme(legend.position=c(.85,.85))
# ggsave(paste0("figures/moransI_SES_PD_gbm_", n_spatial, "spatialpredictor.png"), width=5, height=4, units = "in", dpi = 300)
# 
# 
# ## model results ----
# pdl <- lapply(res, summary)
# 
# pdm <- sapply(pdl, "[", "rel.inf")
# pdv <- sapply(pdl, "[", "var")
# pddf <- data.frame(variable = unlist(pdv), rel.inf = unlist(pdm))
# 
# # sort factors
# su <- tapply(pddf$rel.inf, pddf$variable, median)
# pddf$variable <- factor(pddf$variable)
# pddf$variable <- factor(pddf$variable, levels=names(sort(su, decreasing = F)))
# 
# ggplot(pddf, aes(x=variable, y=rel.inf))+
#   geom_boxplot()+
#   ylab("Relative influence")+
#   theme(axis.title.y = element_blank(), axis.text.y = element_text(size=9))+
#   coord_flip()
# ggsave(paste0("figures/varImp_SES_PD_gbm10_runs_", n_spatial, "spatial.png"), width=6, height=4, units = "in", dpi = 600)
# 
# 
# 
# 
## analyse SES.PD ----
dat <- na.omit(dat)
dat.ns <- st_drop_geometry(dat)

n_spatial <- c(0,1,2)
res0 <- readRDS(paste0("idata/res10GBM_PD", n_spatial[1], "spat_predictors.rds"))
res1 <- readRDS(paste0("idata/res10GBM_PD", n_spatial[2], "spat_predictors.rds"))
res2 <- readRDS(paste0("idata/res10GBM_PD", n_spatial[3], "spat_predictors.rds"))
res <- list(res0, res1, res2)

## morans I for residuals ----
resids_list <- lapply(res, lapply, residuals)
resids_df <- lapply(resids_list, as.data.frame, col.names=c(1:10))
names(resids_df) <- c("spat0", "spat1", "spat2")
PD_residual_median <- lapply(resids_df, apply, 1, median)
resids_df <- as.data.frame(resids_df)
resids_df$SES.PD <- dat.ns$SES.PD

resids_df_l <- tidyr::pivot_longer(resids_df, cols=grep("spat", names(resids_df)))
resids_df_l$name <- gsub("\\.X[0-9]{1,2}", "", resids_df_l$name)
ggplot(resids_df_l, aes(x=SES.PD, y=value, group=SES.PD))+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=45))+
  geom_hline(yintercept=0)+
  ylab(paste0("SES.PD GBM model residuals ", n_spatial, "spat pred, 10 runs"))

# distance matrix
xy <- cbind(dat.ns$x, dat.ns$y)
dm <- dist(xy, method = "euclidean", diag = TRUE, upper = TRUE)
dm <- as.matrix(dm)
#distance thresholds (same units as distance_matrix = METERS)
distance.thresholds <- 100000*c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)
# distance bands
moran <- data.frame(dist.class = distance.thresholds,
                    moransI0 = NA,
                    moransp0 = NA,
                    moransI1 = NA,
                    moransp1 = NA,
                    moransI2 = NA,
                    moransp2 = NA)

coo <- cbind(dat.ns$y, dat.ns$x)
# add residual means to dat.ns dataframe
dat.ns <- cbind(dat.ns, as.data.frame(PD_residual_median))
dat.ns2 <- tidyr::pivot_longer(dat.ns[,c("SES.PD","spat0", "spat1", "spat2")], cols=c(2,3,4), values_to="residual")


for(i in 1:length(moran$dist.class)){
  S.dist  <-  spdep::dnearneigh(coo, 0, moran$dist.class[i], longlat = FALSE)
  lw <- spdep::nb2listw(S.dist, style="W",zero.policy=T) 
  MI0 <- spdep::moran.mc(dat.ns$spat0, lw, nsim=999,zero.policy=T) 
  MI1 <- spdep::moran.mc(dat.ns$spat1, lw, nsim=999,zero.policy=T) 
  MI2 <- spdep::moran.mc(dat.ns$spat2, lw, nsim=999,zero.policy=T) 
  moran$moransI0[i] <- MI0$statistic
  moran$moransp0[i] <- MI0$p.value
  moran$moransI1[i] <- MI1$statistic
  moran$moransp1[i] <- MI1$p.value
  moran$moransI2[i] <- MI2$statistic
  moran$moransp2[i] <- MI2$p.value
  if(!i%%1)cat(i,"\r")
}

moran$obs_id <- seq(1:nrow(moran))
setDT(moran)
dat.mlt <- melt(moran, id.vars = c("obs_id", "dist.class"))
# Ericks comment: now extract group number into new column. I know you can do this more elegantly
# with your grepl magic. I prolly could too, but you get the point.
dat.mlt[grepl("0", variable), group := 0]
dat.mlt[grepl("1", variable), group := 1]
dat.mlt[grepl("2", variable), group := 2]

dat.mlt[, variable := gsub("[0-9]", "", variable)]
dat.mlt.cst <- dcast(data = dat.mlt,
                     obs_id + dist.class + group ~ variable,
                     value.var = "value")

moran_PD <- ggplot(dat.mlt.cst, aes(x=dist.class, y=moransI, group=group, col=factor(group))) +
  geom_line()+
  geom_point(aes(size=moransp<0.05), alpha=0.7)+
  scale_x_continuous("Distance class (m)")+
  scale_color_scico_d("spat. predictors", palette="batlow", end=0.7)+
  scale_size_discrete("p<0.05")+
  ylab("Moran's I GBM residuals")+
  theme(legend.position=c(.8,.8), legend.spacing=unit(1, "mm"), 
        legend.key.height=unit(4, "mm"), legend.box="horizontal", legend.margin=margin(c(1,1,1,1),unit="mm"))
#ggsave(paste0("figures/moransI_SES_PD_gbm_all_spatialpredictor.png"), width=4, height=4, units = "in", dpi = 300)



## analyse model results ----
pdl <- lapply(res, lapply, summary)
pdm <- lapply(pdl, sapply, "[", "rel.inf")
pdv <- lapply(pdl, sapply, "[", "var")
pddf <- data.frame(variable = unlist(lapply(pdv, unlist)), 
                   rel.inf = unlist(lapply(pdm, unlist)),
                   spat_predictors=c(rep(0, length(unlist(pdv[[1]]))),
                                     rep(1, length(unlist(pdv[[2]]))),
                                     rep(2, length(unlist(pdv[[3]])))))
# sort factors
su <- tapply(pddf$rel.inf, pddf$variable, median)
pddf$variable <- factor(pddf$variable)
pddf$variable <- factor(pddf$variable, levels=names(sort(su, decreasing = F)))
varimp_PD <- ggplot(pddf, aes(x=variable, y=rel.inf, fill=factor(spat_predictors)))+
  geom_boxplot()+
  ylab("Relative influence")+
  scale_fill_scico_d(palette="batlow", end=0.7, alpha=0.7)+
  theme(axis.title.y = element_blank())+
  coord_flip()+
  facet_wrap(~spat_predictors)+ 
  theme(strip.background=element_blank(), legend.position="none", panel.grid.major.y=element_line(size=0.05, color="grey50"))
ggsave(paste0("figures/varImp_SES_PD_gbm10_runs_all_spatial.png"), width=6, height=4, units = "in", dpi = 600)


# alt layout cumulative importance plot
ggplot(pddf)+
  geom_boxplot(aes(x=variable, y=rel.inf, fill=factor(spat_predictors)))+
  ylab("Relative influence")+
  scale_fill_scico_d(palette="batlow", end=0.7, alpha=0.7)+
  theme(axis.title.y = element_blank())+
  coord_flip()+
#  facet_wrap(~spat_predictors)+
  theme(strip.background=element_blank(), legend.position=c(.8,.2), panel.grid.major.y=element_line(size=0.05, color="grey50"))

# pdm <- lapply(pdl, sapply, "[", "rel.inf")
# # switch to cumulative importance
# pdv <- lapply(pdl, sapply, "[", "var")
# pddf <- data.frame(variable = unlist(lapply(pdv, unlist)), 
#                    rel.inf = unlist(lapply(pdm, unlist)),
#                    spat_predictors=c(rep(0, length(unlist(pdv[[1]]))),
#                                      rep(1, length(unlist(pdv[[2]]))),
#                                      rep(2, length(unlist(pdv[[3]])))))
dfdf.alt <- aggregate(pddf$rel.inf, by=list(pddf$spat_predictors, pddf$variable), median)
ggplot(dfdf.alt)+
  geom_smooth(aes(y=x, x=Group.2, col=factor(Group.1), fill=factor(Group.1), group=Group.1), alpha=0.1)+
  ylab("Relative influence")+
  scale_color_scico_d(palette="batlow", end=0.7, alpha=0.7)+
  scale_fill_scico_d(palette="batlow", end=0.7, alpha=0.7)+
  theme(axis.text.x=element_text(angle=45, hjust=1))


### performance ----
get_fin_gbm <- function(x){x$results[as.numeric(rownames(x$bestTune)),]}
res.gbm <- lapply(res, lapply, get_fin_gbm)
mat <- matrix(ncol=10, nrow=10)
r2 <- lapply(res.gbm, sapply, "[[", "Rsquared")
rmse <- lapply(res.gbm, sapply, "[[", "RMSE")
PDrf <- data.frame(spat_pred = c(rep(0,10,),rep(1,10,),rep(2,10,)),
                   r2= unlist(r2),
                   rmse = unlist(rmse))
PDrf.l <- data.table::melt(PDrf, id="spat_pred")

performance_PD <- ggplot(PDrf.l, aes(x=factor(spat_pred), y=value, fill=factor(spat_pred), group=spat_pred))+
  geom_boxplot()+
  scale_y_continuous("performance metric")+
  scale_x_discrete("Spatial predictors")+
  scale_fill_scico_d(palette="batlow", end=0.7, alpha=0.7)+
  facet_wrap(~variable, scales="free_y", labeller=as_labeller(c("r2"="Rsquared", "rmse"="RMSE")))+
  theme(strip.background=element_blank(), legend.position="none")
ggsave("figures/GBM_SES.PD_PDrformance.png", width=4, height=3, unit="in", dpi=300)


plot_grid(plot_grid(moran_PD, performance_PD, ncol=2, rel_widths=c(2,1), labels=c("A", "B"), label_fontface=1, label_size=12),
          varimp_PD, ncol=1, labels=c("", "C"), label_fontface=1, label_size=12, rel_heights=c(1,1.5))
ggsave("figures/GBM_SES.PD_combiplot.png", width=20, height=15, unit="cm", dpi=300)




## analyse SES.PE ----
n_spatial <- c(0,1,2)
res0 <- readRDS(paste0("idata/res10GBM_PE", n_spatial[1], "spat_predictors.rds"))
res1 <- readRDS(paste0("idata/res10GBM_PE", n_spatial[2], "spat_predictors.rds"))
res2 <- readRDS(paste0("idata/res10GBM_PE", n_spatial[3], "spat_predictors.rds"))
res <- list(res0, res1, res2)

## morans I for residuals ----
# get model averages
resids_list <- lapply(res, lapply, residuals)
resids_df <- lapply(resids_list, as.data.frame, col.names=c(1:10))
names(resids_df) <- c("spat0", "spat1", "spat2")
PE_residual_median <- lapply(resids_df, apply, 1, median)
resids_df <- as.data.frame(resids_df)
resids_df$SES.PE <- dat.ns$SES.PE

resids_df_l <- tidyr::pivot_longer(resids_df, cols=grep("spat", names(resids_df)))
resids_df_l$name <- gsub("\\.X[0-9]{1,2}", "", resids_df_l$name)
ggplot(resids_df_l, aes(x=SES.PE, y=value, group=SES.PE))+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=45))+
  geom_hline(yintercept=0)+
  ylab(paste0("SES.PE GBM model residuals ", n_spatial, "spat pred, 10 runs"))

# distance matrix
# xy <- cbind(dat$x, dat$y)
# dm <- dist(xy, method = "euclidean", diag = TRUE, upper = TRUE)
# dm <- as.matrix(dm)
# 
# #distance thresholds (same units as distance_matrix = METERS)
# distance.thresholds <- 100000*c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)

# distance bands
moran <- data.frame(dist.class = distance.thresholds,
                    moransI0 = NA,
                    moransp0 = NA,
                    moransI1 = NA,
                    moransp1 = NA,
                    moransI2 = NA,
                    moransp2 = NA)

coo <- cbind(dat$y, dat$x)
# add residual means to dat.ns dataframe
dat.ns <- cbind(dat.ns, as.data.frame(PE_residual_median))
dat.ns2 <- tidyr::pivot_longer(dat.ns[,c("SES.PD","spat0", "spat1", "spat2")], cols=c(2,3,4), values_to="residual")


for(i in 1:length(moran$dist.class)){
  S.dist  <-  spdep::dnearneigh(coo, 0, moran$dist.class[i], longlat = FALSE)
  lw <- spdep::nb2listw(S.dist, style="W",zero.policy=T) 
  MI0 <- spdep::moran.mc(dat.ns$spat0, lw, nsim=999,zero.policy=T) 
  MI1 <- spdep::moran.mc(dat.ns$spat1, lw, nsim=999,zero.policy=T) 
  MI2 <- spdep::moran.mc(dat.ns$spat2, lw, nsim=999,zero.policy=T) 
  moran$moransI0[i] <- MI0$statistic
  moran$moransp0[i] <- MI0$p.value
  moran$moransI1[i] <- MI1$statistic
  moran$moransp1[i] <- MI1$p.value
  moran$moransI2[i] <- MI2$statistic
  moran$moransp2[i] <- MI2$p.value
  if(!i%%1)cat(i,"\r")
}

moran$obs_id <- seq(1:nrow(moran))
setDT(moran)
dat.mlt <- melt(moran, id.vars = c("obs_id", "dist.class"))

# Ericks comment: now extract group number into new column. I know you can do this more elegantly
# with your grepl magic. I prolly could too, but you get the point.
dat.mlt[grepl("0", variable), group := 0]
dat.mlt[grepl("1", variable), group := 1]
dat.mlt[grepl("2", variable), group := 2]

dat.mlt[, variable := gsub("[0-9]", "", variable)]
dat.mlt.cst <- dcast(data = dat.mlt,
                     obs_id + dist.class + group ~ variable,
                     value.var = "value")

moran_pe <- ggplot(dat.mlt.cst, aes(x=dist.class, y=moransI, group=group, col=factor(group))) +
  geom_line()+
  geom_point(aes(size=moransp<0.05), alpha=0.7)+
  scale_x_continuous("Distance class (m)")+
  scale_color_scico_d("spat. predictors", palette="batlow", end=0.7)+
  scale_size_discrete("p<0.05")+
  ylab("Moran's I GBM residuals")+
  theme(legend.position=c(.8,.8), legend.spacing=unit(1, "mm"), 
        legend.key.height=unit(4, "mm"), legend.box="horizontal", legend.margin=margin(c(1,1,1,1),unit="mm"))
#ggsave(paste0("figures/moransI_SES_PE_gbm_all_spatialpredictor.png"), width=4, height=4, units = "in", dpi = 300)



## analyse model results ----
pdl <- lapply(res, lapply, summary)
pdm <- lapply(pdl, sapply, "[", "rel.inf")
pdv <- lapply(pdl, sapply, "[", "var")
pddf <- data.frame(variable = unlist(lapply(pdv, unlist)), 
                   rel.inf = unlist(lapply(pdm, unlist)),
                   spat_predictors=c(rep(0, length(unlist(pdv[[1]]))),
                                     rep(1, length(unlist(pdv[[2]]))),
                                     rep(2, length(unlist(pdv[[3]])))))
# sort factors
su <- tapply(pddf$rel.inf, pddf$variable, median)
pddf$variable <- factor(pddf$variable)
pddf$variable <- factor(pddf$variable, levels=names(sort(su, decreasing = F)))
varimp_pe <- ggplot(pddf, aes(x=variable, y=rel.inf, fill=factor(spat_predictors)))+
  geom_boxplot()+
  ylab("Relative influence")+
  scale_fill_scico_d(palette="batlow", end=0.7, alpha=0.7)+
  theme(axis.title.y = element_blank())+
  coord_flip()+
  facet_wrap(~spat_predictors)+ 
  theme(strip.background=element_blank(), legend.position="none", panel.grid.major.y=element_line(size=0.05, color="grey50"))
ggsave(paste0("figures/varImp_SES_PE_gbm10_runs_all_spatial.png"), width=6, height=4, units = "in", dpi = 600)

# alt layout: 



### performance ----
get_fin_gbm <- function(x){x$results[as.numeric(rownames(x$bestTune)),]}
res.gbm <- lapply(res, lapply, get_fin_gbm)
mat <- matrix(ncol=10, nrow=10)
r2 <- lapply(res.gbm, sapply, "[[", "Rsquared")
rmse <- lapply(res.gbm, sapply, "[[", "RMSE")
perf <- data.frame(spat_pred = c(rep(0,10,),rep(1,10,),rep(2,10,)),
                   r2= unlist(r2),
                   rmse = unlist(rmse))
perf.l <- melt(perf, id="spat_pred")

performance_pe <- ggplot(perf.l, aes(x=factor(spat_pred), y=value, fill=factor(spat_pred), group=spat_pred))+
  geom_boxplot()+
  scale_y_continuous("Performance metric")+
  scale_x_discrete("Spatial predictors")+
  scale_fill_scico_d(palette="batlow", end=0.7, alpha=0.7)+
  facet_wrap(~variable, scales="free_y", labeller=as_labeller(c("r2"="Rsquared", "rmse"="RMSE")))+
  theme(strip.background=element_blank(), legend.position="none")
ggsave("figures/GBM_SES.PE_performance.png", width=4, height=3, unit="in", dpi=300)


plot_grid(plot_grid(moran_pe, performance_pe, ncol=2, rel_widths=c(2,1), labels=c("A", "B"), label_fontface=1, label_size=12),
          varimp_pe, ncol=1, labels=c("", "C"), label_fontface=1, label_size=12, rel_heights=c(1,1.5))
ggsave("figures/GBM_SES.PE_combiplot.png", width=20, height=15, unit="cm", dpi=300)


# Figure for both 2 spatial predictor model ----
respd <- readRDS(paste0("idata/res10GBM_PD2spat_predictors.rds"))
respe <- readRDS(paste0("idata/res10GBM_PE2spat_predictors.rds"))

pdl <- lapply(respd, summary)
pel <- lapply(respe, summary)

pdm <- sapply(pdl, "[", "rel.inf")
pdv <- sapply(pdl, "[", "var")
pddf <- data.frame(variable = unlist(pdv), rel.inf = unlist(pdm))

pem <- sapply(pel, "[", "rel.inf")
pev <- sapply(pel, "[", "var")
pedf <- data.frame(variable = unlist(pev), rel.inf = unlist(pem))

df <- aggregate(pddf$rel.inf, by=list(variable=pddf$variable), FUN=median)
df2 <- aggregate(pedf$rel.inf, by=list(variable=pedf$variable), FUN=median)

df <- merge(df, df2, by="variable", all=T)
names(df) <- c("variable", "SES.PD", "SES.PE")
df.mlt <- melt(df, id.vars="variable", variable.name="group", value.name="rel.inf")
# sort factors
su <- tapply(df.mlt$rel.inf, df.mlt$variable, median)
df.mlt$variable <- factor(df.mlt$variable)
df.mlt$variable <- factor(df.mlt$variable, levels=names(sort(su, decreasing = F)))

varimp_bar <- ggplot(df.mlt, aes(x=factor(variable), y=rel.inf, fill=group))+
  geom_bar(stat="identity", position="dodge", width=0.7)+
  ylab("Relative influence")+
  scale_fill_scico_d(palette="batlow", end=0.7, alpha=0.7)+
  theme(axis.title.y = element_blank())+
  coord_flip()+
  theme(legend.position=c(0.9,0.8), legend.title=element_blank())
ggsave(varimp_bar, paste0("figures/varImp_gbm_best_models_summary.png"), width=4, height=3, units = "in", dpi = 300)

varimp_ranks <- ggplot(df, aes(x=rank(SES.PD), y=rank(SES.PE), label=variable))+
  #geom_point()+
  geom_label(label.padding=unit(0,"mm"), label.size=0, size=1.7)+
  scale_x_continuous("SES.PD")+
  scale_y_continuous("SES.PE")+
  theme(text=element_text(size=6), plot.background=element_blank(), panel.background=element_blank())+
  ggtitle("Variable model ranks")
ggsave(varimp_ranks, paste0("figures/var_ranks_summary.png"), width=4, height=4, units = "in", dpi = 300)

ggdraw() + draw_plot(varimp_bar, 0, 0, 1, 1) + draw_plot(varimp_ranks, 0.51, 0.1, 0.45, 0.65)
ggsave(paste0("figures/varimp_combi.png"), width=6, height=4, units = "in", dpi = 300)


# marginal effect plots ----

plot_grid(ncol=4, 
plot.gbm(res[[1]][[1]]$finalModel, "tra_mean", ylab="SES.PE"),
plot.gbm(res[[1]][[1]]$finalModel, "elev_range"),
plot.gbm(res[[1]][[1]]$finalModel, "sub_trop_mbf"),
plot.gbm(res[[1]][[1]]$finalModel, "mat_mean"),
plot.gbm(res[[1]][[1]]$finalModel, "pre_sd"),
plot.gbm(res[[1]][[1]]$finalModel, "area"),
plot.gbm(res[[1]][[1]]$finalModel, "soil"))









# Classification tree for Hotspots via environment ------------------------

load("workspace_maps.RData")

# clean unneeded data
shp2 <- shp2[,-grep("obs_p|obs_rank|reps|_cw$|LEVEL2|LEVEL1|LEVEL_3_CO|LEVEL_NAME|ID|\\.3|_rw|CONTI|REGION|mpd|TTD|AvTD|^bio_|^FC|^PC", names(shp2))]

#coordinates of the cases
shp2$centroids <- st_centroid(shp2) %>% 
  st_coordinates()
shp2$y <- shp2$centroids[,1] # x=lng
shp2$x <- shp2$centroids[,2]
xy <- st_drop_geometry(shp2[, c("x", "y")])
shp2.ns <- st_drop_geometry(shp2)

#distance matrix
dm <- dist(xy, method = "euclidean", diag = TRUE, upper = TRUE)
dm <- as.matrix(dm)

#distance thresholds (same units as distance_matrix)
distance.thresholds <- 100000*c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)

## get spatial predictors ----
mems <- spatialRF::mem_multithreshold(
  distance.matrix = dm,
  distance.thresholds = distance.thresholds
)
kableExtra::kbl(
  head(mems[, 1:4], n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

mem.rank <- spatialRF::rank_spatial_predictors(
  distance.matrix = dm,
  spatial.predictors.df = mems,
  ranking.method = "moran"
)
mems <- mems[, mem.rank$ranking]
shp2.ns <- cbind(shp2.ns, mems[,1:10])
shp2.ns$hotspot_type[is.na(shp2.ns$hotspot_type)] <- "none"
shp2.ns$hotspot_type <- as.factor(shp2.ns$hotspot_type)


## tuning parameter with caret package ----
gbmControl <- trainControl(method = "repeatedcv", number = 10,
                           repeats = 3, savePredictions = "final",
                           returnResamp ="final")# allowParallel=F has become necessary recently, not sure why.
gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3),
                        n.trees = (1:10)*50,
                        shrinkage = c(0.1, 0.01, 0.001),
                        n.minobsinnode = c(5, 10, 15))
# define gbm function
run.gbm <- function(i, seeds, data, trControl=gbmControl,
                    tuneGrid=gbmGrid, f){
  set.seed(s[i])
  if(!i%%1)cat(i,"\r")
  temp <- train(f, data, method = "gbm",
                trControl = gbmControl, verbose=FALSE,
                tuneGrid = gbmGrid,
                na.action=na.omit)
  return(temp)
}

## define gbm formula. choose here spatial predictors if desired ----
vars <- paste(names(shp2.ns)[c(12:40,78,79)], collapse="+")
f.hs <- as.formula(paste0("hotspot_type ~ ", vars))

## RUN GBMs (sequential run) ----
foreach::registerDoSEQ()
s <- seq(500,509,1)
res.class.tree <- list()
for(i in 1:1){
  res.class.tree[[i]] <- run.gbm(i=i, seeds=s[i], data=shp2.ns, f=f.hs)
  beepr::beep(2)
}
tmp <- res.class.tree[[1]]
dev.off()
summary(tmp)
tmp
tmp$bestTune
tmp$results[184,] # 3% correctly identified when taking random chance into account. lol


# Accuracy is the percentage of correctly classified instances out of all
# instances. It is more useful on a binary classification than multi-class
# classification problems because it can be less clear exactly how the accuracy
# breaks down across those classes (e.g. you need to go deeper with a confusion
# matrix).
#
# Kappa is like classification accuracy, except that it is
# normalized at the baseline of random chance on your dataset. It is a more
# useful measure to use on problems that have an imbalance in the classes (e.g.
# 70-30 split for classes 0 and 1 and you can achieve 70% accuracy by predicting
# all instances are for class 0).

## Confusion matrix
expected <- tmp$pred$obs
predicted <- tmp$pred$pred
results <- confusionMatrix(data=predicted, reference=expected)
print(results)
table(expected)


# Simple binomial model
shp2.ns$hs <- "no"
shp2.ns$hs[shp2.ns$hotspot_type %in% c("PE","PD","both")] <- "yes"
shp2.ns$hs <- as.factor(shp2.ns$hs)

nums <- c(12:40,78,79)
vars <- paste(names(shp2.ns)[nums[1:30]], collapse="+")
f.hs <- as.formula(paste0("hs ~ ", vars))

dat <- shp2.ns[,names(shp2.ns) %in% unlist(strsplit(as.character(f.hs), split=" "))]
dat <- na.omit(dat)
b1 <- glm(f.hs, data=shp2.ns, family="binomial")
summary(b1)
b1.steps <- MASS::stepAIC(b1, steps=20)
b1.steps

b2 <- glm(hs ~ soil + mat_mean + mat_sd + tra_sd + prs_sd + 
            mio_mat_ano_mean + tri + sub_trop_mbf + temp_bmf + flooded_gs, 
          family = "binomial", data = shp2.ns)
summary(b2)

b3 <- glm(hs ~ soil + sqrt(tra_sd) + mat_mean + sub_trop_mbf , dat, family="binomial")
summary(b3)

# problem: not normal
