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
