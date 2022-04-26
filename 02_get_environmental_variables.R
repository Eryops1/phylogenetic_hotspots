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

# merge into
shp <- merge(shp, global_drivers, by.x="LEVEL3_COD", by.y="LEVEL_3_CO", all.x=TRUE)


# Human footprint ------------------------------------------------------

# https://www.nature.com/articles/s41597-022-01284-8#Sec12
hfp <- raster("../DATA/PDiv/hfp2018.tif")
behr <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
hfp <- projectRaster(hfp, crs=behr)


# Deforestation -----------------------------------------------------------

tx <- readLines("../DATA/PDiv/lossyear.txt")
for(i in 1:length(tx)){
  url <- tx[i]
  destfile <- paste0("../DATA/PDiv/deforestation/", sub("^.*?GFC-2020-v1\\.8/", "", tx[i]))
  download.file(url, destfile)
}

deforest <- raster("../DATA/PDiv/deforestation/deforest_combined.tif")
deforest[deforest == 255] <- NA


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

future <- lapply(destfiles, raster)
future <- stack(future)
names(future) <- regmatches(destfiles, regexpr("CHELSA_bio[1-9]{1,2}", destfiles))


# past and projected future land use change



