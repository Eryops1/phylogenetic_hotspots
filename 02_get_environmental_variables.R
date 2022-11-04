# variables to compare to phylogeny hotspots


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
library(data.table) # fast csv reading
library(castor) # fast tree reading
library(phyloregion) # PD calculations
library(raster)
library(sf)
library(ggplot2)
theme_set(theme_bw())
library(cowplot)
library(beepr)
library(ncdf4)


# *** Human footprint ------------------------------------------------------

# source: https://www.nature.com/articles/s41597-022-01284-8#Sec12
hfp <- raster("data/hfp2018.tif")


# *** Deforestation -----------------------------------------------------------

# https://data.globalforestwatch.org/documents/tree-cover-loss/explore

# tx <- readLines("data/deforestation/lossyear.txt")
# for(i in 1:length(tx)){
#   url <- tx[i]
#   destfile <- paste0("data/deforestation/", sub("^.*?GFC-2020-v1\\.8/", "", tx[i]))
#   if(!file.exists(destfile)){
#   download.file(url, destfile)}
# }
# 
# tx <- readLines("data/deforestation/treecover2000.txt")
# for(i in 1:length(tx)){
#   url <- tx[i]
#   destfile <- paste0("data/treecover2000/", sub("^.*?GFC-2021-v1\\.9/", "", tx[i]))
#   if(!file.exists(destfile)){
#     download.file(url, destfile)}
# }

# combine in some GIS and then save as one tif:
treecover <- raster("data/treecover300dpi.tif")
treecover <- raster::crop(treecover, extent(-180, 180, -59, 78))
deforest <- raster("data/deforest_combined.tif")



# convert to binary tree cover yes/no for clipping 
x <- values(treecover)
hist(x)
x[x<10] <- 0
x[x>=10] <- 1
values(treecover) <- x
plot(treecover)

plot(deforest) # different values are different years!!! its ALL deforestation
extent(deforest)
deforest <- raster::crop(deforest, extent(-17355314, 17328826, -6353170, 7174430))
x <- values(deforest)
x[x<10] <- 0
x[x>=10] <- 1
values(deforest) <- x
plot(deforest)

# when calculating the deforestation per bot country, crop bot country layer
# with with tree cover layer first to get only the values for area actually
# covered with forest


writeRaster(deforest, "data/deforest_fin.tif", overwrite=T)
writeRaster(treecover, "data/treecover_fin.tif", overwrite=T)

# *** Climate change ----------------------------------------------------------

# https://www.worldclim.org/data/cmip6/cmip6_clim30s.html, CMIP6
# https://www.carbonbrief.org/cmip6-the-next-generation-of-climate-models-explained
# 
# These scenarios – RCP2.6, RCP4.5, RCP6.0, and RCP8.5 – have new versions in
# CMIP6. These updated scenarios are called SSP1-2.6, SSP2-4.5, SSP4-6.0, and
# SSP5-8.5. 
# --> worst case (SSP5-8.5), middle of the road (SSP3-7.0) and more
# optimistic (SSP4-6.0) Using IPSL-CM6A-LR ssp370 as relatively conservative in
# regards of ECS, + available via Chelsa (years 2071-2100)

# BIO1 = Annual Mean Temperature
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO12 = Annual Precipitation
# BIO 15 = precipitation seasonality (!)

urls <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2071-2100/IPSL-CM6A-LR/ssp370/bio/CHELSA_bio", c(1,4,5,6,7,12,15), "_2071-2100_ipsl-cm6a-lr_ssp370_V.2.1.tif")
destfiles <- paste0("data/climate_change/", sub("^.*?ssp370/bio/", "", urls))
for(i in 1:length(urls)){
  if(!file.exists(destfiles[i])){
  download.file(urls[i], destfiles[i], method="wget", quiet=TRUE, options(timeout = max(600, getOption("timeout"))))
  # use method="wget" + quiet=TRUE to use increased timeout option (bio12 is >500mb)
  }
  print(i)
}



# # save data for cluster
# save("shp", "hfp", "deforest", "CHELSA_bio1", "CHELSA_bio5",
#     "CHELSA_bio6", "CHELSA_bio7", "CHELSA_bio12", "CHELSA_bio15",
#      "pes2", "fes2",
#      file="environment_vars.RData")




## OUTSOURCING to cluster for time efficiency --------------

# run: environment_vars_job_array.sh > env.sh > env_var.R





# Different method for deforestation --------------------
shp <- readRDS("data/fin_shape.rds")
deforest <- terra::rast("data/deforest_fin.tif")
treecover <- terra::rast("data/treecover_fin.tif")
shp2 <- terra::vect(shp)

# set up variables
var="deforest"
vars_stat <- c('proportion')
combs <- nrow(expand.grid(var, vars_stat))
m <- matrix(seq(1:combs), ncol=3, byrow = TRUE)
num_list <- split(m, rep(1:nrow(m)))
res <- matrix(nrow=nrow(shp), ncol=combs)
rownames(res) <- shp$LEVEL3_COD
disag_id <- c()
upsale_count <- c()


for(i in 1:nrow(shp2)){
  # loop over botanical countries
    shape_sub <- shp2[i]
  # kill switch for ANT
    if(shape_sub$LEVEL3_COD=="ANT"){next}
    treecover_crop <- terra::crop(treecover, shape_sub, mask=TRUE)
    tree_poly <- terra::as.polygons(treecover_crop)
  # kill switch for no treecover:
    if(all(dim(tree_poly)==c(0,0))){next}
    
    # get the area affected by deforestation
    deforest_crop <- terra::crop(deforest, tree_poly, mask=TRUE)
    x <- as.numeric(na.omit(values(deforest_crop)))
    affected <- length(which(x==1))
    unaffected <- length(which(x==0))
    total <- length(x)
    
    # get mean, sd and sample size
    res[i,1] <- affected/total
    # res[i,2] <- sd(rest$value)
    # res[i,3] <- length(rest$value)
    print(i)
}

res <- rbind(res, NA)
row.names(res)[369] <- "BOU"
res <- res[order(row.names(res)),]
res <- as.matrix(res)

#saveRDS(res, file=paste0('data/deforestation2.rds'))





# Conservation hotpot coverages for botanical countries ---------------------
h <- st_read("data/hotspots_fixed.gpkg")
h <- st_wrap_dateline(h, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))

shp <- readRDS("data/fin_shape.rds")
s <- as(st_geometry(shp), "Spatial")
m3 <- as(st_geometry(h), "Spatial")


library(terra)
a <- c()
for(i in 1:length(s)){
  tmp <- s[i,]
  res <- raster::intersect(tmp, m3) # zero width buffering happens automatically now
  if(class(res)!="NULL"){
    hot_area <- sum(area(res))
    tmp_area <- area(tmp)
    a <- c(a, hot_area / tmp_area)
  }else{a <- c(a, 0)}
  if(!i%%1)cat(i,"\r")
}

shp$hotspot_coverage <- a
saveRDS(shp, file='data/hotspot_coverage.rds')





# Merge + save ------------------------------------------------------------



shp <- readRDS("fin_shape.rds")

## Load data from Global Plant Diversity Drivers --------------------------
global_drivers <- readRDS("../Global_Drivers_2022/processed_data/shp_object_fin_analysis.RDS")
global_drivers <- global_drivers[,-grep("ID|_NAM|CONTINENT|sr|pet_|_n|_abs|mangroves", names(global_drivers))]
global_drivers <- sf::st_drop_geometry(global_drivers)

# merge into shape object
shp <- merge(shp, global_drivers, all.x=TRUE)


## Load Environment data ---------------------------------------------------


vars <- c('hfp','deforestation2','bio1','bio5','bio6','bio7','bio12','bio15')

var.list <- lapply(paste0("data/", vars, ".rds"), readRDS)
names(var.list) <- vars
env.df <- do.call(cbind.data.frame, var.list)
env.df$LEVEL3_COD <- row.names(env.df)

# merge into shape object
shp <- merge(shp, env.df, all.x=TRUE)
names(shp)

# BIO1 = Annual Mean Temperature
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO12 = Annual Precipitation
# Bio15 = Precipitation Seasonality

## Add hotspot coverage --------------------------------------------------
hc <- readRDS('data/hotspot_coverage.rds')
shp <- merge(shp, st_drop_geometry(hc), all.x=T)


## Get future MAT + PRE change -------------------------------------------

shp$mat_change <- shp$bio1.1 - shp$mat_mean
shp$pre_change <- shp$bio12.1 - shp$pre_mean # future - current = increase if positive



# transform to Behrmann projection
#behr <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
#tmp <- st_transform(shp, behr) # Behrmann

saveRDS(tmp, "data/fin_shp.rds")
