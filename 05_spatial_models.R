# old school spatial models 


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
library(raster)
library(sf)
library(ggplot2)
theme_set(theme_bw())
library(cowplot)
library(rgdal)
library(mgcv)
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

tmp <- shp[,-grep("rand_|hfp|deforest|bio|change", names(shp))]
# tmp$centroids <- st_centroid(tmp) %>% 
#   st_coordinates()
# tmp$y <- tmp$centroids[,1] # x=lng
# tmp$x <- tmp$centroids[,2]
#tmp <- st_drop_geometry(tmp)
tmp <- tmp[,-grep("centroids|PD_obs|PE_sd|PE|WE|mdr", names(tmp))]




# SpatialRF ---------------------------------------------------------------

library(spatialRF)
#loading training data and distance matrix from the package
#data(plant_richness_df)
#data(distance_matrix)

# omit NAs
shp2 <- shp[,c(1:48)]
shp2 <- na.omit(shp2)

#names of the response variable and the predictors
dependent.variable.name <- "SES.PD"
predictor.variable.names <- colnames(shp2)[16:48]

#coordinates of the cases
shp2$centroids <- st_centroid(shp2) %>% 
   st_coordinates()
shp2$y <- shp2$centroids[,1] # x=lng
shp2$x <- shp2$centroids[,2]
xy <- st_drop_geometry(shp2[, c("x", "y")])

#distance matrix
#dist(x, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)
dm <- dist(xy, method = "euclidean", diag = TRUE, upper = TRUE)
dm <- as.matrix(dm)
#dm <- st_distance(shp, which="Euclidean") # this is extremely slow

# Weighted distance matrix
distMat <- as.matrix(dist(cbind(shp2$x, shp2$y)))
diag(distMat) <- 0


#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 100000, 500000, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000)

#random seed for reproducibility
random.seed <- 1

spatialRF::plot_training_df(
  data = shp2,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  ncol = 5,
  point.color = viridis::viridis(100, option = "F"),
  line.color = "gray30"
)

## assess autocorrelation

spatialRF::plot_training_df_moran(
  data = shp2,
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

# Reducing multicollinearity in the predictors

preference.order <- c(
  "mat_mean",
    "sub_trop_mbf",
  "tra_mean",
  "pre_mean",
  "soil",
  "area",
  "pre_lgm_ano_mean")

predictor.variable.names2 <- spatialRF::auto_cor(
  x = st_drop_geometry(shp2[, predictor.variable.names]),
  cor.threshold = 0.85,
  preference.order = preference.order
) %>%
  spatialRF::auto_vif(
    vif.threshold = 5,
    preference.order = preference.order
  )

predictor.variable.names <- predictor.variable.names2

names(predictor.variable.names)

# Finding promising variable interactions
interactions <- spatialRF::the_feature_engineer(
  data = st_drop_geometry(shp2),
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  xy = xy,
  importance.threshold = 0.50, #uses 50% best predictors
  cor.threshold = 0.60, #max corr between interactions and predictors
  seed = random.seed,
  repetitions = 100,
  verbose = TRUE
)


## non spatial RF
model.non.spatial <- spatialRF::rf(
  data = st_drop_geometry(shp2),
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

#plot of climate_bio1_average
p1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = shp,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = local.importance,
    ggplot2::aes(
      x = y,
      y = x,
      color = mrd
    )
  ) +
  ggplot2::scale_color_gradient2(
    low = color.low,
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ggtitle("mrd") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) +
  ggplot2::labs(color = "Importance") +
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude")

p2 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = shp,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = local.importance,
    ggplot2::aes(
      x = y,
      y = x,
      color = soil
    )
  ) +
  ggplot2::scale_color_gradient2(
    low = color.low,
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ggtitle("soil") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) +
  ggplot2::labs(color = "Importance") +
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude")

p1 + p2
# In these maps, values lower than 0 indicate that for a given record, the permuted version of the variable led to an accuracy score even higher than the one of the non-permuted variable, so again these negative values can be interpreted as “worse than chance”.
spatialRF::plot_response_curves(
  model.non.spatial,
  quantiles = 0.5,
  ncol = 3
)

# Model performance
spatialRF::print_performance(model.non.spatial)

# Spatial cross-validation
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
#The plot below shows the Moran’s I of the residuals of the spatial model, and indicates that the residuals are not autocorrelated at any distance.

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
  head(model.spatial$importance$per.variable, n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# look at spatial predictors
spatial.predictors <- spatialRF::get_spatial_predictors(model.spatial)
pr <- data.frame(spatial.predictors, st_drop_geometry(shp2[, c("x", "y")]))

p1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = shp2, fill = "white") +
  ggplot2::geom_point(
    data = pr,
    ggplot2::aes(
      x = y,
      y = x,
      color = spatial_predictor_1000000_11
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "Eigenvalue") +
  ggplot2::ggtitle("Variable: spatial_predictor_1000000_11") +
  ggplot2::theme(legend.position = "bottom")+
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude")

p2 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = shp2, fill = "white") +
  ggplot2::geom_point(
    data = pr,
    ggplot2::aes(
      x = y,
      y = x,
      color = spatial_predictor_2000000_5,
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "Eigenvalue") +
  ggplot2::ggtitle("Variable: spatial_predictor_2000000_5") +
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
model.full <- rf_spatial(
  data = st_drop_geometry(shp2),
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = dm,
  distance.thresholds = distance.thresholds,
  xy = xy
) %>%
  rf_tuning() %>%
  rf_evaluate() %>%
  rf_repeat()


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

# GAM ---------------------------------------------------------------------

library(spdep)
# create neighbor list
tmp <- na.omit(tmp)
nb <- spdep::poly2nb(tmp, row.names = tmp$LEVEL3_COD, )
names(nb) <- tmp$LEVEL3_COD
#names(nb) <- attr(nb, "region.id")

# example
data(oldcol)
col.W <- nb2listw(COL.nb, style="W")
crime <- COL.OLD$CRIME
lee.test(crime, crime, col.W, zero.policy=TRUE)
#A positive Lee’s L indicates that clusters match for the two variables. A
#negative value indicates that the clusters have an opposite spatial
#distribution . A value around zero indicates that the spatialstructures of the
#two variables do not match. The significance of thevalues of Lee’s L was
#evaluated using a Monte Carlo test with 999 randomizations.
col.W <- nb2listw(nb, style="W", zero.policy = TRUE)
lee.test(tmp$SES.PD_ts, tmp$SES.PD_ts, col.W, zero.policy = TRUE)
# SES.PD_ts shows significant spatial autocorrelation 
lee.test(tmp$SES.PD_ts, tmp$SES.PD_ts, col.W, zero.policy = TRUE)


tmp$LEVEL3_COD <- as.factor(tmp$LEVEL3_COD)

# model 
ctrl <- gam.control(nthreads = 6) # parallel threads

gam_mrf <- gam(SES.PD_ts ~ s(LEVEL3_COD, bs = 'mrf', xt = list(nb = nb)) + # define MRF smooth (markov random fields)
                 mrd +soil + mat_mean + tra_mean + pre_mean + prs_mean + area + elev_range,
               data = st_drop_geometry(tmp),
               method = 'REML', 
               family = betar,  # fit a beta regression
               control = ctrl) 
summary(gam_mrf)


l1 <- lm(data=shp2, log(richness)~SES.PD*SES.PE)
summary(l1)
b <- mgcv::gam(data=shp2, log(richness)~s(SES.PD)+s(SES.PE))
summary(b)
plot(b,pages=1,residuals=TRUE)  ## show partial residuals
mgcv::gam.check(b)
mgcv::vis.gam(b,theta=30,phi=30,ticktype="detailed")