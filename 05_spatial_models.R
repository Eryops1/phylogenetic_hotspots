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
shp <- cbind(PD_obs=shp$PD_obs_ts, shp)
shp <- shp[,-grep("obs_ts|obs_rw|obs_cw|obs_p|obs_rank|reps|_cw$|LEVEL2|LEVEL1|LEVEL_3_CO|LEVEL_NAME|ID|\\.3|_rw|CONTI|REGION", names(shp))]
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
