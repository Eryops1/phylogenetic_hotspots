# variables to compare to phylogeny hotspots



wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
#library(arrow) # fast parquet file reading
library(data.table) # fast csv reading
library(castor) # fast tree reading
library(phyloregion) # PD calculations
library(raster)
library(sf)
library(ggplot2)
theme_set(theme_bw())
library(cowplot)
library(rgdal)
library(beepr)
library(ncdf4)


# Measures of long-term climate stability (e.g. Miocene climate anomalies),
# present environment (climate, soils), topography, geographic isolation etc.
# will be included as potential drivers of phylogenetic diversity/endemism.
# Measures of threat, such as past and projected future land use change,
# deforestation, and climate change will be included to identify phylogenetic
# diversity hotspots.


# Prep for cluster runs  ---------------------------------------------------



# *** Human footprint ------------------------------------------------------

# https://www.nature.com/articles/s41597-022-01284-8#Sec12
hfp <- raster("../DATA/PDiv/hfp2018.tif")
behr <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
hfp <- projectRaster(hfp, crs=behr)
writeRaster(hfp, "../DATA/PDiv/hfp.tif")


# *** Deforestation -----------------------------------------------------------

# https://data.globalforestwatch.org/documents/tree-cover-loss/explore
tx <- readLines("../DATA/PDiv/lossyear.txt")
for(i in 1:length(tx)){
  url <- tx[i]
  destfile <- paste0("../DATA/PDiv/deforestation/", sub("^.*?GFC-2020-v1\\.8/", "", tx[i]))
  if(!file.exists(destfile)){
  download.file(url, destfile)}
}

deforest <- raster("../DATA/PDiv/deforestation/deforest_combined.tif")
deforest[deforest == 255] <- NA
deforest <- projectRaster(deforest, crs=behr)
writeRaster(deforest, "../DATA/PDiv/deforest.tif")

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
destfiles <- paste0("../DATA/PDiv/climate_change/", sub("^.*?ssp370/bio/", "", urls))
for(i in 1:length(urls)){
  if(!file.exists(destfiles[i])){
  download.file(urls[i], destfiles[i], method="wget", quiet=TRUE, options(timeout = max(600, getOption("timeout"))))
  # use method="wget" + quiet=TRUE to use increased timeout option (bio12 is >500mb)
  }
  print(i)
}



# *** Land use change past + future -------------------------------------------
## land use change 1960 - 2019: urban, cropland, rangeland/pastures
## https://doi.org/10.1594/PANGAEA.921846

# use now: 
# https://luh.umd.edu/data.shtml
# primf: forested primary land
# primn: non-forested primary land
# secdf: potentially forested secondary land
# secdn: potentially non-forested secondary land
# pastr: managed pasture
# range: rangeland
# urban: urban land
# c3ann: C3 annual crops
# c3per: C3 perennial crops
# c4ann: C4 annual crops
# c4per: C4 perennial crops
# c3nfx: C3 nitrogen-fixing crops
# secma: secondary mean age (units: years)
# secmb: secondary mean biomass density (units: kg C/m^2)

# pasT: 
tmp <- nc_open("../DATA/PDiv/land_use/states_hist.nc")
ncatt_get(tmp,"time","units")
vars <- names(tmp$var)[1:14]

pes <- list()
for(i in 1:length(vars)){pes[[i]] <- brick("../DATA/PDiv/land_use/states_hist.nc", varname=vars[i])}
names(pes) <- vars

funclist <- function(x, years){x[[1166]]-x[[1166-years]]}
pes2 <- lapply(pes, funclist, years=100)
rm(pes)
for(i in 1:length(vars)){writeRaster(pes2[[i]], paste0("../DATA/PDiv/land_use/PC_", vars[i], ".tif"))}

# futurE:
tmp <- nc_open("../DATA/PDiv/land_use/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-AIM-ssp370-2-1-f_gn_2015-2100.nc")
ncatt_get(tmp,"time","units")
vars <- names(tmp$var)[1:14]

fes <- list()
for(i in 1:length(vars)){fes[[i]] <- brick(
  "../DATA/PDiv/land_use/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-AIM-ssp370-2-1-f_gn_2015-2100.nc",
  varname=vars[i])}
names(fes) <- vars

futlist <- function(x, years){x[[1]]-x[[1+years]]}
fes2 <- lapply(fes, futlist, years=85)
rm(fes)
for(i in 1:length(vars)){writeRaster(fes2[[i]], paste0("../DATA/PDiv/land_use/FC_", vars[i], ".tif"))}






# save data for cluster
# save("shp", "hfp", "deforest", "CHELSA_bio1", "CHELSA_bio5",
#     "CHELSA_bio6", "CHELSA_bio7", "CHELSA_bio12", "CHELSA_bio15",
#      "pes2", "fes2",
#      file="environment_vars.RData")


# *** Build scripts -----------------------------------------------------------
# read .tif files

# build R scripts
Rscript <- "
# get environmental variable layer
args <- commandArgs()
print(args)
var <- args[6]

library(raster)
library(sf)
library(rgdal)
library(exactextractr)

# load data
shp <- readRDS('fin_shape.rds')
lay <- raster(paste0(var, '.tif')
  
# set up variables
vars_stat <- c('mean', 'sd', 'n')
combs <- nrow(expand.grid(var, vars_stat))
m <- matrix(seq(1:combs), ncol=3, byrow = TRUE)
num_list <- split(m, rep(1:nrow(m)))

res <- matrix(nrow=nrow(shp@data), ncol=combs)
rownames(res) <- shp@data$LEVEL3_COD
disag_id <- c()
upsale_count <- c()

for(i in 1:nrow(shp@data)){
  # loop over botanical countries
  shape_sub <- subset(shp, shp$LEVEL3_COD==shp$LEVEL3_COD[[i]])
    rest <- exact_extract(lay, shape_sub)
    #rest <- raster::extract(lay, shape_sub)
    rest <- na.omit(rest[[1]])
    # increase resolution necessary?
    if(all(is.na(rest))==TRUE){
      print('disaggregate to increase resolution')
      upsale_count <- c(upsale_count, 1)
      disag_id <- c(disag_id, paste(i))
      # to avoid huge raster files crop to extent of shapefile sub * 20
      newExtent <- extent(bbox(shape_sub))
      lay2 <- crop(lay, newExtent*20)
      lay2 <- disaggregate(lay2, 10)
      rest <- exact_extract(lay2, shape_sub)
      rest <- na.omit(rest[[1]])
      
    }else{}
    
    # get mean, sd and sample size
    res[i,1] <- mean(rest$value)
    res[i,2] <- sd(rest$value)
    res[i,3] <- length(rest$value)
  print(i)}
  
saveRDS(res, file=paste0(var,'.rds'))
"
cat(Rscript, file="env_var.R")

### build subordination bash script
bashscript2 <- '#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=env_var_child
#SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=8gb
#SBATCH --cpus-per-task 1
#SBATCH --time 02:00:00
##SBATCH --output=NULL

source ~/miniconda3/bin/activate pdiv
Rscript env_var.R $this_var > logfile.txt
'
cat(bashscript2, file="env.sh")

### build bash array script
bashscript <- "#!/bin/bash
# submit_array.sh

#SBATCH --account PDiv
#SBATCH --job-name=env_var_boss
#SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=1gb
#SBATCH --cpus-per-task 1
#SBATCH --time 00:05:00

vars=('hfp' 'deforest' 'bio1' 'bio5' 'bio6' 'bio7' 'bio12')

# pass on variables to child scripts
for ((i=0; i<=6; i++)) do
  this_var=${vars[${i}]}
  export this_var
  sbatch env.sh
done
"
cat(bashscript, file="environment_vars_job_array.sh")

system({
  "git add .
git commit -m 'env model script changes'
git push"
})

# stuff on GenomeDK: git pull etc (doesnot work yet....)
# system2({
#   "ssh -T 'mtietje@login.genome.au.dk'"
# })
# #system("ssh -T 'mtietje@login.genome.au.dk'")








# Merge + save ------------------------------------------------------------



shp <- readRDS("fin_shape.rds")
shp <- st_as_sf(shp)
# transform to Behrmann projection
# shp <- st_transform(shp, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") # Behrmann


## Load data from Global Plant Diversity Drivers --------------------------
global_drivers <- readRDS("../Global_Drivers_2022/processed_data/shp_object_fin_analysis.RDS")
global_drivers <- global_drivers[,-grep("ID|_NAM|CONTINENT|sr|pet_|_n|_abs|mangroves", names(global_drivers))]
global_drivers <- sf::st_drop_geometry(global_drivers)

# merge into shape object
shp <- merge(shp, global_drivers, by.x="LEVEL3_COD", by.y="LEVEL_3_CO", all.x=TRUE)


## Load Environment data ---------------------------------------------------

vars <- c('hfp','deforest','bio1','bio5','bio6','bio7','bio12','bio15','PC_primf','PC_primn','PC_secdf','PC_secdn','PC_urban','PC_c3ann','PC_c4ann','PC_c3per','PC_c4per','PC_c3nfx','PC_pastr','PC_range','PC_secmb','PC_secma','FC_primf','FC_primn','FC_secdf','FC_secdn','FC_urban','FC_c3ann','FC_c4ann','FC_c3per','FC_c4per','FC_c3nfx','FC_pastr','FC_range','FC_secmb','FC_secma')

var.list <- lapply(paste0("environment/", vars, ".rds"), readRDS)
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


## Get future MAT + PRE change -------------------------------------------

shp$mat_change <- shp$bio1.1 - shp$mat_mean
shp$pre_change <- shp$bio12.1 - shp$pre_mean # future - current = increase if positive



# SAVE --------------------------------------------------------------------
# saveRDS(shp, "shp.rds")
# 
# 
# some shapefile wrangling since new TDWG shapefile is corrupt....-> GIT##
# rm(list = ls())
# so <- readOGR("../DATA/shapefile_bot_countries/level3.shp")
# s <- readRDS("shp.rds")
# so <- merge(so@data, s, all.x=TRUE, by.x="LEVEL_3_CO", by.y="LEVEL3_COD")
# 
# # transform to Behrmann projection
# behr <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
# strans <- spTransform(so, CRS(behr))
# 
# shp <- st_as_sf(so)
# 
# #m = st_buffer(shp, 0)  ## fixes some issues
# shp <- st_make_valid(shp)
# st_is_valid(shp)
# 
# 
# shp <- st_crop(shp, st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 90)))
# transform to Behrmann projection
behr <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
tmp <- st_transform(shp, behr) # Behrmann

saveRDS(tmp, "fin_shp.rds")
