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


# Measures of long-term climate stability (e.g. Miocene climate anomalies),
# present environment (climate, soils), topography, geographic isolation etc.
# will be included as potential drivers of phylogenetic diversity/endemism.
# Measures of threat, such as past and projected future land use change,
# deforestation, and climate change will be included to identify phylogenetic
# diversity hotspots.


# Load data from Global Plant Diversity Drivers ---------------------------

global_drivers <- readRDS("../Global_Drivers_2022/processed_data/shp_object_fin_analysis.RDS")
global_drivers <- global_drivers[,-grep("ID|_NAM|CONTINENT|sr|pet_|_n|_abs|mangroves", names(global_drivers))]
global_drivers <- sf::st_drop_geometry(global_drivers)



# Human footprint ------------------------------------------------------

# https://www.nature.com/articles/s41597-022-01284-8#Sec12
hfp <- raster("../DATA/PDiv/hfp2018.tif")
behr <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
hfp <- projectRaster(hfp, crs=behr)
writeRaster(hfp, "../DATA/PDiv/hfp.tif")


# Deforestation -----------------------------------------------------------

# https://data.globalforestwatch.org/documents/tree-cover-loss/explore
# tx <- readLines("../DATA/PDiv/lossyear.txt")
# for(i in 1:length(tx)){
#   url <- tx[i]
#   destfile <- paste0("../DATA/PDiv/deforestation/", sub("^.*?GFC-2020-v1\\.8/", "", tx[i]))
#   download.file(url, destfile)
# }

#deforest <- raster("../DATA/PDiv/deforestation/Hansen_GFC-2020-v1.8_lossyear_00N_050W.tif")
deforest <- raster("../DATA/PDiv/deforestation/deforest_combined.tif")
deforest[deforest == 255] <- NA
deforest <- projectRaster(deforest, crs=behr)
writeRaster(deforest, "../DATA/PDiv/deforest.tif")

# Climate change ----------------------------------------------------------

# https://www.worldclim.org/data/cmip6/cmip6_clim30s.html, CMIP6
# https://www.carbonbrief.org/cmip6-the-next-generation-of-climate-models-explained
# 
# These scenarios – RCP2.6, RCP4.5, RCP6.0, and RCP8.5 – have new versions in
# CMIP6. These updated scenarios are called SSP1-2.6, SSP2-4.5, SSP4-6.0, and
# SSP5-8.5. 
# --> worst case (SSP5-8.5), middle of the road (SSP3-7.0) and more
# optimistic (SSP4-6.0) Using IPSL-CM6A-LR ssp370 as relatively conservative in
# regards of ECS, + avaialaible via Chelsa (years 2071-2100)

# BIO1 = Annual Mean Temperature
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO12 = Annual Precipitation

urls <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2071-2100/IPSL-CM6A-LR/ssp370/bio/CHELSA_bio", c(1,4,5,6,7,12), "_2071-2100_ipsl-cm6a-lr_ssp370_V.2.1.tif")
destfiles <- paste0("../DATA/PDiv/climate_change/", sub("^.*?ssp370/bio/", "", urls))
# for(i in 1:6){
#   download.file(urls[i], destfiles[i], method="wget", quiet=TRUE, options(timeout = max(600, getOption("timeout"))))
#   # use method="wget" + quiet=TRUE to use increased timeout option (bio12 is >500mb)
#   print(i)
# }

# future <- lapply(destfiles, raster)
# names(future) <- regmatches(destfiles, regexpr("CHELSA_bio[1-9]{1,2}", destfiles))
# #list2env(future, envir=.GlobalEnv)
# for(i in 1:length(future)){
#   #future[[i]] <- projectRaster(future[[i]], crs=behr)
#   writeRaster(future[[i]], paste0("../DATA/PDiv/", names(future[i]), ".grd"))
#   print(i)
# }
# dont do this, too big. just use the original tifs on the cluster

# past and projected future land use change

# STILL MISSING #
#
#
#
#


# Average per LEVEL3 unit -------------------------------------------------

shp <- readRDS("fin_shape.rds")

# save data for cluster
# save("shp", "hfp", "deforest", "CHELSA_bio1", "CHELSA_bio5",
#               "CHELSA_bio6", "CHELSA_bio7", "CHELSA_bio12",
#      file="environment_vars.RData")

### build R scripts
Rscript <- "
# get environmental variable layer
args <- commandArgs()
print(args)
var <- args[6]

library(raster)
library(sf)
library(rgdal)

# load data
shp <- readRDS('fin_shape.rds'')
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
    rest <- raster::extract(lay, shape_sub)
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
      rest <- raster::extract(lay2, shape_sub)
      rest <- na.omit(rest[[1]])
      
    }else{}
    
    # get mean, sd and sample size
    res[i,1] <- mean(rest)
    res[i,2] <- sd(rest)
    res[i,3] <- length(rest)
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
shp <- st_as_sf(s)
# transform to Behrmann projection
shp <- st_transform(shp, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") # Behrmann


# merge into PD file
shp <- merge(shp, global_drivers, by.x="LEVEL3_COD", by.y="LEVEL_3_CO", all.x=TRUE)


