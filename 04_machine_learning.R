# machine learning instead of linear assumptions



# Neural net --------------------------------------------------------------

library(neuralnet)

# Binary classification
nn <- neuralnet(Species == "setosa" ~ Petal.Length + Petal.Width, iris, linear.output = FALSE)
print(nn)
plot(nn)

# Multiclass classification
nn <- neuralnet(Species ~ Petal.Length + Petal.Width, iris, linear.output = FALSE)
print(nn)
plot(nn)

# Custom activation function
softplus <- function(x) log(1 + exp(x))
nn <- neuralnet((Species == "setosa") ~ Petal.Length + Petal.Width, iris, 
                linear.output = FALSE, hidden = c(3, 2), act.fct = softplus)
print(nn)
plot(nn)

confidence.interval(nn)
gwplot(nn)

# that seems to be a classification thing for factorial data
nn2 <- neuralnet(SES.PD_ts ~ mrd + soil + richness + x + y, tmp, linear.output = F, hidden = c(5))
print(nn2)
plot(nn2)
summary(nn2)


# GPBoost -----------------------------------------------------------------
# https://towardsdatascience.com/tree-boosting-for-spatial-data-789145d6d97d


library(gpboost)
data(GPBoost_data, package = "gpboost")

#--------------------Grouped random effects model: single-level random effect----------------
gp_model <- GPModel(group_data = group_data[,1], likelihood="gaussian")

#--------------------Gaussian process model----------------
gp_model <- GPModel(gp_coords = coords, cov_function = "exponential",
                    likelihood="gaussian")

#--------------------Combine Gaussian process with grouped random effects----------------
gp_model <- GPModel(group_data = group_data,
                    gp_coords = coords, cov_function = "exponential",
                    likelihood="gaussian")


# gp_model <- GPModel(gp_coords = coords, 
#                     cov_function = "exponential")
# Training
bst <- gpboost(data = X, label = y, gp_model = gp_model,
               nrounds = 247, learning_rate = 0.01,
               max_depth = 3, min_data_in_leaf = 10, 
               num_leaves = 2^10,
               objective = "regression_l2", verbose = 0)
summary(gp_model) # Estimated covariance parameters
# Make predictions: latent variables and response variable
pred <- predict(bst, data = X_test, gp_coords_pred = coords_test, 
                predict_var = TRUE, pred_latent = TRUE)
# pred[["fixed_effect"]]: predictions from the tree-ensemble.
# pred[["random_effect_mean"]]: predicted means of the gp_model.
# pred["random_effect_cov"]]: predicted (co-)variances 
# of the gp_model
pred_resp <- predict(bst, data = X_test, 
                     gp_coords_pred = coords_test, 
                     pred_latent = FALSE)
y_pred <- pred_resp[["response_mean"]] # predicted response mean
# Calculate mean square error
MSE_GPBoost <- mean((y_pred-y_test)^2)



library("SHAPforxgboost")
shap.plot.summary.wrap1(bst, X = X)
shap_long <- shap.prep(bst, X_train = X)
shap.plot.dependence(data_long = shap_long, x = "Covariate_1",
                     color_feature = "Covariate_2", smooth = FALSE)



# Create random effects model and datasets
gp_model <- GPModel(gp_coords = coords, 
                    cov_function = "exponential") 
dtrain <- gpb.Dataset(data = X, label = y)
# Candidate parameter grid
param_grid = list("learning_rate" = c(1,0.1,0.01),
                  "min_data_in_leaf" = c(1,10,100),
                  "max_depth" = c(1,3,5,10))
# Other parameters not contained in the grid of tuning parameters
params <- list(objective = "regression_l2", verbose = 0, 
               "num_leaves" = 2^10)
# Use random grid search and cross-validation. Set 'num_try_random=NULL' to use deterministic grid search
set.seed(1)
opt_params <- gpb.grid.search.tune.parameters(
  param_grid = param_grid,
  params = params,
  num_try_random = 20,
  nfold = 4,
  data = dtrain,
  gp_model = gp_model,
  verbose_eval = 1,
  nrounds = 1000,
  early_stopping_rounds = 5,
  eval = "l2")
# Found the following parameters:
# ***** New best score (0.397766105827059) found for the following parameter combination: learning_rate: 0.01, min_data_in_leaf: 10, max_depth: 3, nrounds:247



# Spatial GBM -------------------------------------------------------------

library(spm)
data(sponge)
set.seed(1234)
rvi1 <- rvi(sponge[, -c(3)], sponge[, 3], family = "poisson", n.cores=2)

gbmokcv1 <- gbmokcv(sponge[, c(1,2)], sponge[,-c(3)], sponge[, 3],
                    cv.fold = 10, family = "poisson", n.cores=2, predacc = "ALL")
gbmokcv1

n <- 20 # number of iterations, 60 to 100 is recommended.
VEcv <- NULL
for (i in 1:n) {
  gbmokcv1 <- gbmokcv(sponge[, c(1,2)], sponge[, -c(3)], sponge[, 3],
                      cv.fold = 10, family = "poisson", n.cores=4, predacc = "ALL")
  VEcv [i] <- gbmokcv1
}
plot(VEcv ~ c(1:n), xlab = "Iteration for gbmok", ylab = "VEcv (%)")
points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
abline(h = mean(VEcv), col = 'blue', lwd = 2)

rvi(sponge[, -c(3)], sponge[, 3], family = "poisson", n.cores=2)

# my data
library(sf)
shp <- readRDS("fin_shp.rds")

# remove not needed data
shp <- cbind(PD_obs=shp$PD_obs_ts, shp)
shp <- shp[,-grep("obs_ts|obs_rw|obs_cw|obs_p|obs_rank|reps|_cw$|LEVEL2|LEVEL1|LEVEL_3_CO|LEVEL_NAME|ID|\\.3|_rw|CONTI|REGION", names(shp))]
names(shp)<- gsub("\\.1", "_mean", names(shp))
names(shp)<- gsub("\\.2", "_sd", names(shp))
# remove BOU that has no data
shp <- shp[!shp$LEVEL3_COD=="BOU",]

tmp <- shp[,-grep("LEVEL3_COD|rand_|hfp|deforest|bio|change", names(shp))]
tmp$centroids <- st_centroid(tmp) %>% 
  st_coordinates()
tmp$y <- tmp$centroids[,1] # x=lng
tmp$x <- tmp$centroids[,2]
tmp <- st_drop_geometry(tmp)
tmp <- tmp[,-grep("centroids|PD_obs|PE_sd|PE|WE|mdr", names(tmp))]
# regular GBM
tmp <- na.omit(tmp)
rvi2 <- rvi(tmp[, -2], tmp[, 2], family = "gaussian", n.cores=4)

gbmokcv2 <- gbmokcv(tmp[, c("x", "y")], tmp[, -2], tmp[, 2], family = "gaussian", n.cores=4,
                    cv.fold = 10, predacc = "ALL")
gbmokpred1 <- gbmokpred(tmp[, c("x", "y")], tmp[, c("x", "y", 6:9)], petrel[, 3],
                        petrel.grid[, c(1,2)], petrel.grid, family = "gaussian", n.cores=6,
                        nmax = 12, vgm.args = ("Sph"))
