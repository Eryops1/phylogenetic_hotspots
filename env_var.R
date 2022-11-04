
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
lay <- raster(paste0(var, '.tif'))
  
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
#    rest <- raster::extract(lay, shape_sub)
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
