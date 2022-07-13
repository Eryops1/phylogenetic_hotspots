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


# compare spatial and non spatial models
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


# Get spatial predictors for other models--------------------------------
#single distance (0km by default)
mems <- spatialRF::mem(distance.matrix = dm)

#several distances
mems <- spatialRF::mem_multithreshold(
  distance.matrix = dm,
  distance.thresholds = distance.thresholds
)

# But not all MEMs are made equal, and you will need to rank them by their
# Moran’s I. The function rank_spatial_predictors() will help you do so.

mem.rank <- spatialRF::rank_spatial_predictors(
  distance.matrix = dm,
  spatial.predictors.df = mems,
  ranking.method = "moran"
)
#The output of rank_spatial_predictors() is a list with three slots: “method”, a
#character string with the name of the ranking method; “criteria”, an ordered
#data frame with the criteria used to rank the spatial predictors; and
#“ranking”, a character vector with the names of the spatial predictors in the
#order of their ranking (it is just the first column of the “criteria” data
#frame). We can use this “ranking” object to reorder or mems data frame.

mems <- mems[, mem.rank$ranking]

# From here, spatial predictors can be included in any model one by one, in the
# order of the ranking, until the spatial autocorrelation of the residuals
# becomes neutral, if possible. A little example with a linear model follows.

#model definition
predictor.variable.names <- colnames(shp2)[16:48]
predictors <- predictor.variable.names

model.formula <- as.formula(
  paste(
    dependent.variable.name,
    " ~ ",
    paste(
      predictors,
      collapse = " + "
    )
  )
)

#scaling the data
#model.data <- scale(st_drop_geometry(shp2)) %>%
#  as.data.frame()

#fitting the model
m <- lm(model.formula, data = shp2)

#Moran's I test of the residuals
moran.test <- spatialRF::moran(
  x = residuals(m),
  distance.matrix = dm,
  verbose = FALSE
)
print(moran.test$plot)


# GBM -------------------------------------------------------------------
library(caret)
library(gbm)
library(parallel)
# parameter tuning using caret
gbmControl <- trainControl(method = "repeatedcv", number = 10,
                           repeats = 3, savePredictions = "final",
                           returnResamp ="final")

gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3),
                        n.trees = (1:10)*50,
                        shrinkage = c(0.1, 0.01, 0.001),
                        n.minobsinnode = c(5, 10, 15))


# define function
## formula
dat_no.na <- st_drop_geometry(shp2)
n <- names(dat_no.na[,-grep("LEVEL", names(dat_no.na))])
do_not_include <- c("SES.PD", "centroids", "y", "x")
f_ses.pd <- as.formula(paste("SES.PD ~", paste(n[!n %in% do_not_include], collapse = " + ")))

run.gbm <- function(i, seeds, data, trControl=gbmControl,
                    tuneGrid=gbmGrid, form=form){
  set.seed(s[i])
  if(!i%%1)cat(i,"\r")
  temp <- train(form, data, method = "gbm",
                trControl = gbmControl, verbose=FALSE,
                tuneGrid = gbmGrid)
  return(temp)
}

n_cores=1
s <- seq(500,507,1) # models are rather slow with that many variables, ca. 5 minutes per model
system.time(
  PD_list <- mclapply(1:length(s), run.gbm, seeds = s,
                      data = dat_no.na, mc.cores=n_cores, form=f_ses.pd)
)

# this throws errors currently, no fucking idea why.... single rungs outside parallel are fine

#pd_list <- unlist(PD_list, recursive = F)
pdl <- lapply(PD_list, varImp)
pdm <- as.data.frame(sapply(pdl, "[", "importance"))
su <- rowSums(pdm)
pdm$pred = row.names(pdm) 
pddf <- reshape2::melt(pdm, id.name= "pred", value.name = "varImp")
# sort factors
pddf$pred <- factor(pddf$pred)
pddf$pred <- factor(pddf$pred, levels=names(sort(su, decreasing = F)))

ggplot(pddf, aes(x=pred, y=varImp))+
  geom_boxplot()+
  ylab("Variable importance")+
  #  facet_wrap(~Var2, scales = "free")+
  #  theme(axis.text.x = element_text(size=6, angle = 45, hjust=1),
  #        strip.background = element_blank())+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size=9))+
  coord_flip()
ggsave("figures/varImp_SES_PD_gbm10_runs.png", width=7, height=10, units = "in", dpi = 600)



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
