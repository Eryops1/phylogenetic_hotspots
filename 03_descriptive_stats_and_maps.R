
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
library(phyloregion) # PD calculations
library(raster)
library(sf)
library(ggplot2)
theme_set(theme_bw())
library(cowplot)
library(rgdal)
library(beepr)
if(!dir.exists("figures"))dir.create("figures")



# Load data ---------------------------------------------------------------

shp <- readRDS("fin_shp.rds")

# remove not needed data
shp <- cbind(PD_obs=shp$PD_obs_ts, shp)

shp <- shp[,-grep("obs_ts|obs_rw|obs_cw|obs_p|obs_rank|reps|_cw$|LEVEL2|LEVEL1|LEVEL_3_CO|LEVEL_NAME|ID|\\.3|_rw|CONTI|REGION", names(shp))]
names(shp)<- gsub("\\.1", "_mean", names(shp))
names(shp)<- gsub("\\.2", "_sd", names(shp))

# remove BOU that has no data
shp <- shp[!shp$LEVEL3_COD=="BOU",]


# Explorative plots -------------------------------------------------------

plot(st_drop_geometry(shp)[,c(1,4:12,13:23)])
plot(st_drop_geometry(shp)[,c(1,4:12,24:35)])
plot(st_drop_geometry(shp)[,c(1,4:12,36:57)])




#### COMPLETENESS and CORRELATION ########################################

cor.test(shp$richness, shp$SES.PD_ts)
cor.test(shp$richness, shp$PD_obs)

dat <- st_drop_geometry(shp)

# missing data
library(naniar)
library(UpSetR)
vis_miss(dat)

# how many missings?
n_var_miss(dat)
gg_miss_upset(dat, nsets=15)

dat.sub <- dat[,-grep("_sd$", names(dat))]
n_var_miss(dat.sub)
gg_miss_upset(dat.sub, nsets=9)

plot_grid(ggplot(dat.sub, aes(y=area, group=is.na(hfp_mean)))+
  geom_boxplot(varwidth = T)+
  scale_y_log10()+
  xlab("HFP == NA (positive=T: n=27)"),
ggplot(dat.sub, aes(y=PD_obs, group=is.na(hfp_mean)))+
  geom_boxplot(varwidth = T)+
  xlab("HFP == NA (positive=T: n=27)"),
ggplot(dat.sub, aes(y=SES.PD_ts, group=is.na(hfp_mean)))+
  geom_boxplot(varwidth = T)+
  xlab("HFP == NA (positive=T: n=27)"),
ggplot(dat.sub, aes(y=PE, group=is.na(hfp_mean)))+
  geom_boxplot(varwidth = T)+
  scale_y_log10()+
  xlab("HFP == NA (positive=T: n=27)"))
ggsave("figures/missing_hfp_influence.png", width=5, height=4, 
       units = "in", dpi = 300, bg = "white")


dat_no.na <- na.omit(dat.sub)
dim(dat_no.na)

source("99_functions.R")
library(rstatix)
cmat <- cor_mat(dat_no.na[,-grep("LEVEL3|hfp|deforest|bio|change",names(dat_no.na))], method = "s")
cpmat <- cor_pmat(dat_no.na[,-grep("LEVEL3|hfp|deforest|bio|change",names(dat_no.na))], method = "s")
#cor.dat_no.na<- cor(dat_no.na[,-grep("LEVEL3|hfp|deforest|bio|change",names(dat_no.na))], method = "s")
#p.dat_no.na <- cor_pmat(dat_no.na[,-grep("LEVEL3|hfp|deforest|bio|change",names(dat_no.na))], method = "s")

# make names more readable
fn <- colnames(cpmat)
fn <- gsub("mdr", "DR", fn)
fn <- gsub("mrd", "MRD", fn)
fn <- gsub("mio_pre_ano_mean", "mio_pre_ano", fn)
fn <- gsub("mio_mat_ano_mean", "mio_mat_ano", fn)
fn <- gsub("mat_lgm_ano_mean", "lgm_mat_ano", fn)
fn <- gsub("pre_lgm_ano_mean", "lgm_pre_ano", fn)
fn <- gsub("sub_trop_", "(sub)trop ", fn)
colnames(cmat) <- fn
row.names(cmat) <- fn


my.corrplot(cmat, lab=T, p.mat = cpmat, insig = "blank",
            tl.cex = 8, digits = 1, type = "lower", hc.order = F,
            lab_size = 4, highlight = TRUE, method="square", lab_col = "grey20",
            colors = c("#6874DE", "white", "#D80000"), 
            legend.title = "Spearman rho")+
  scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text(margin=margin(0,0,0,0)),  # Order: top, right, bottom, left
        axis.text.y = element_text(margin=margin(0,0,0,0)),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.grid = element_blank(),
        legend.position = c(.15, .8), legend.text.align = 1)
ggsave("figures/correlation.png", width=7, height=6, units = "in", dpi = 600, bg = "white")



# Spatial autocorrelation -------------------------------------------------


# get coordinates 
shp$centroids <- st_centroid(shp) %>% 
  st_coordinates()
shp$y <- shp$centroids[,1] # x=lng
shp$x <- shp$centroids[,2]

# subset sf object to complete cases??? maybe, check how the sarlm works later.
# use all for now
# shp <- shp[shp$LEVEL3_COD %in% dat_no.na$level3,]
# shp <- merge(shp, dat_no.na[,c("sr_trans", "level3")], by.x="LEVEL_3_CO", by.y="level3")

# check variable distributions
hist(shp)


# Get model residuals
pd <- 
  

zinc.lm <- lm(log(zinc) ~ elev + sqrt(dist), data = meuse)
summary(zinc.lm)
  
fm.p <- sem(fm.npn, data = dat_no.na, estimator="MLM", meanstructure=TRUE)
resids <- as.data.frame(my.residuals_lavaan(fm.p)) 
rawdata <- data.frame(inspect(fm.p, "data")) # get raw data from the model
names(resids) <- paste0(names(resids),"_residuals")
plot(rawdata$sr_trans, dat_no.na$sr_trans) # check order for attaching level3 ID
resids$level3 <- dat_no.na$level3

# add residuals to spatial object
shp <- merge(shp, resids, by.x="LEVEL_3_CO", by.y="level3", all.x=TRUE) 

# add fitted SR values to spatial object
fitted_vals <- my.fitted_lavaan(fm.p)
names(fitted_vals) <- paste0(names(fitted_vals),"_fitted")
fitted_vals$level3 <- dat_no.na$level3 
shp <- merge(shp, fitted_vals[,c("sr_trans_fitted", "level3")], by.x="LEVEL_3_CO", by.y="level3", all.x=TRUE)

# plot original vs fitted and residuals
plot_grid(labels = c("A", "B", "C"), label_fontface = "plain", ncol = 2,
          ggplot(shp, aes(x=sr_trans, y=sr_trans_fitted))+
            geom_point()+
            xlab("species richness")+
            ylab("species richness fitted")+
            geom_abline(slope=1, intercept=0)
          ,
          ggplot(shp, aes(x=sr_trans, y=sr_trans_residuals))+
            geom_point()+
            xlab("species richness")+
            ylab("species richness SEM residuals")+
            geom_abline(slope=0, intercept=0)
          ,
          ggplot(shp, aes(x=abs(y), y=sr_trans_residuals, col=sr_trans))+
            geom_point()+
            ylab("species richness SEM residuals")+
            xlab("absolute latitude")+
            scale_color_continuous("SR")+
            geom_smooth(method="lm")
)
#ggsave(file="../figures/SEM_residuals_SR.png",
#       width=7, height=7, units = "in", dpi = 600)



## Residuals spatial autocorrelation ###################################################

# Weighted distance matrix
distMat <- as.matrix(dist(cbind(shp$x, shp$y)))
distsInv <- 1/distMat # invert matrix for weights
diag(distsInv) <- 0

## GLOBAL autocorrelation
MI.sr <- moran.mc(shp$sr_trans_residuals[,1], mat2listw(distsInv, style = "B"), nsim=599, zero.policy=T)
MI.mrd  <-  moran.mc(shp$mrd_residuals[,1], mat2listw(distsInv, style = "B"), nsim=599, zero.policy=T)
MI.sr$statistic
MI.mrd$statistic

## DISTANCE bands
# distances are in kilometers
moran <- data.frame(dist.class = seq(100, 10000, 100),
                    sr.moransI = NA,
                    sr.moransp = NA,
                    mrd.moransI = NA,
                    mrd.moransp = NA)
coo <- cbind(shp$y, shp$x)
for(i in 1:length(moran$dist.class)){
  S.dist  <-  dnearneigh(coo, 0, moran$dist.class[i], longlat = TRUE)
  lw <- nb2listw(S.dist, style="W",zero.policy=T) 
  
  MI <- moran.mc(shp$sr_trans_residuals[,1], lw, nsim=599,zero.policy=T) 
  moran$sr.moransI[i] <- MI$statistic
  moran$sr.moransp[i] <- MI$p.value
  
  MI.mrd <- moran.mc(shp$mrd_residuals[,1], lw, nsim=599,zero.policy=T) 
  moran$mrd.moransI[i] <- MI.mrd$statistic
  moran$mrd.moransp[i] <- MI.mrd$p.value
  #moran$avg.n[i] <- median(card(lw$neighbours))
  if(!i%%1)cat(i,"\r")
}

temp <- pivot_longer(moran, cols = c("sr.moransI", "mrd.moransI"))
ggplot(temp, aes(x=dist.class, y=value, col=name)) +
  geom_line()+
  scale_x_continuous("Distance class (km)")+
  scale_color_discrete("SEM residuals", labels=c("MRD", "SR"))+
  ylab("Moran's I")
#ggsave("../figures/SAC_sem_residuals.png", width=5, height=4, units = "in", dpi = 600)



# Multicollinearity ########################################################
# check potential full regressions for multicollinearity
names(dat_no.na) <- sub("_mean$", "_m", names(dat_no.na)) # shorten mean to _m

temp <- dat_no.na[,-grep("LEVEL3|hfp|deforest|bio|change|pd_rand|PD_obs|WE|PE",
                         names(dat_no.na))]
lm_pd <- lm(SES.PD_ts ~ . , data=temp)
summary(lm_pd)
sort(car::vif(lm_pd))
# mat, subtropmbf are problematic

# test changes
lm_pd1 <- lm(SES.PD_ts ~ . -sub_trop_mbf, data=temp)
sort(car::vif(lm_pd1))
lm_pd2 <- lm(SES.PD_ts ~ . -mat_m, data=temp)
sort(car::vif(lm_pd2))
lm_pd3 <- lm(SES.PD_ts ~ . -mat_m -sub_trop_mbf, data=temp)
sort(car::vif(lm_pd3))
lm_pd4 <- lm(SES.PD_ts ~ . -mat_m -sub_trop_mbf -pre_lgm_ano_m, data=temp)
sort(car::vif(lm_pd4))





plot_grid(ncol=2,
          ggplot(data=shp, aes(fill=PD_obs))+
            geom_sf(lwd=0)
          ,
          ggplot(data=shp, aes(fill=SES.PD_ts))+
            geom_sf(lwd=0)
          ,
          ggplot(data=shp, aes(x=richness, y=SES.PD_ts))+
            geom_text(aes(label=LEVEL3_COD))
          ,
          ggplot(data=s@data, aes(x=pd_rand_mean_ts, y=PD_obs_ts, col=log(richness)))+
            geom_point()+
            geom_abline(x=1)+
            scale_x_continuous(trans="sqrt")+
            scale_y_continuous(trans="sqrt")
          # observed PD is usually lower than for a random species sample for the phylogey that matches the SR
)


ggplot(data=s@data, aes(x=richness, y=PD_obs_ts))+
  geom_point()+
  geom_point(aes(y=pd_rand_mean_ts), col="salmon")+
  scale_x_continuous(trans="sqrt")+
  scale_y_continuous(trans="sqrt")

hist(shp$pd_rand_mean_ts)



# PD_obs: observed PD in community
# pd_rand_mean: mean PD in null communities
# pd_obs_rank: Rank of observed PD vs. null communities
# pd_obs_z: Standardized effect size of PD vs. null communities = (PD_obs - pd_rand_mean) / pd_rand_sd
# pd_obs_p: P-value (quantile) of observed PD vs. null communities = mpd_obs_rank / iter + 1








# Maps --------------------------------------------------------------------

# format to sf object for plotting
shp <- st_as_sf(s)

# transform to Behrmann projection
shp <- st_transform(shp, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") # Behrmann

(pd_map <- ggplot(shp) + 
    geom_sf(aes(fill=PD_obs_ts),lwd=0, col=NA) + 
    #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
    #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
    #                        begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("PD", option = "plasma")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          #panel.border = element_blank(),
          text = element_text(size = 10),
          # axis.ticks.x = element_line(color="white"),
          # axis.text.x = element_text(color="white"),
          # plot.margin = margin(0, 0, 0, 0.5, "cm")
    )+
    #coord_sf(expand = F, label_axes = "-N--", ylim=c(-6200000, 8200000))+
    xlab(" ")
)
(pd_zscore_map <- ggplot(shp) + 
    geom_sf(aes(fill=zscore_ts),lwd=0, col=NA) + 
    #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
    #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
    #                        begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("PD zscore", option = "plasma")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")
)
(pd_rand_map <- ggplot(shp) + 
    geom_sf(aes(fill=pd_rand_mean_ts),lwd=0, col=NA) + 
    #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
    #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
    #                        begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("PD null distribution", option = "plasma")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")
)

# (pd_ses_map <- ggplot(shp) + 
#     geom_sf(aes(fill=PD_ses),lwd=0, col=NA) + 
#     #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
#     #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
#     #                        begin = lcol, end = sqrt(ucol))+
#     scale_fill_viridis_c("PD_ses", option = "plasma")+ #, 
#     theme(legend.position = c(0.18, 0.3),
#           legend.key.height = unit(6,"mm"),
#           legend.background = element_blank(),
#           legend.key = element_blank(),
#           panel.background = element_blank(),
#           #panel.border = element_blank(),
#           text = element_text(size = 10),
#           # axis.ticks.x = element_line(color="white"),
#           # axis.text.x = element_text(color="white"),
#           # plot.margin = margin(0, 0, 0, 0.5, "cm")
#     )+
#     #coord_sf(expand = F, label_axes = "-N--", ylim=c(-6200000, 8200000))+
#     xlab(" ")
# )

(pe_map <- ggplot(shp) + 
    geom_sf(aes(fill=PE),lwd=0, col=NA) + 
    #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
    #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
    #                        begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("PE", option = "plasma")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          #panel.border = element_blank(),
          text = element_text(size = 10),
          # axis.ticks.x = element_line(color="white"),
          # axis.text.x = element_text(color="white"),
          # plot.margin = margin(0, 0, 0, 0.5, "cm")
    )+
    #coord_sf(expand = F, label_axes = "-N--", ylim=c(-6200000, 8200000))+
    xlab(" ")
)

(sr_map <- ggplot(shp) + 
    geom_sf(aes(fill=richness),lwd=0, col=NA) + 
    #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
    #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
    #                        begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("SR", option = "plasma")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          #panel.border = element_blank(),
          text = element_text(size = 10),
          # axis.ticks.x = element_line(color="white"),
          # axis.text.x = element_text(color="white"),
          # plot.margin = margin(0, 0, 0, 0.5, "cm")
    )+
    #coord_sf(expand = F, label_axes = "-N--", ylim=c(-6200000, 8200000))+
    xlab(" ")
)

plot_grid(pd_map, pd_rand_map, pd_zscore_map, sr_map, nrow = 2)


plot(shp$PD_obs_ts, shp$pd_rand_mean_ts)
plot(shp$PD_obs_ts, shp$PE)

# lm1 <- lm(data=shp, PD_obs_ts~mrd + soil + mat_mean + tra_mean + pre_mean + prs_mean + elev_range+
#      tri + area + hfp.1 + deforest.1 + bio1.1 + bio12.1)




# High threat regions -----------------------------------------------------

# Definition?# regions where deforestation is high, anticipated future climate change is strong, and hfp is strong. 
# Top x percent of each factor?


# Get top countries for each variable -------------------------------------

co <- 0.025 # top 9 countries

# top X % PD_obs

orders <- apply(shp[,c(5:(ncol(shp)-1))], 2, order, decreasing=T)
res <- matrix(nrow=nrow(orders), ncol=ncol(orders))
for(i in 1:ncol(orders)){
  res[,i] <- shp$LEVEL3_COD[orders[,i]]
}

colnames(res) <- colnames(shp[,c(5:(ncol(shp)-1))])

# Top 2.5% in richness, PD observed, PD zscore + PE:
res[1:9, c("richness", "PD_obs_ts", "zscore_ts", "PE")]



# Hotspot definition ------------------------------------------------------





# Endemism hotspots ##
C <- coldspots(shp$PE) # coldspots
H <- hotspots(shp$PE) # hotspots

## Merge endemism values to shapefile of grid cells.
DF <- data.frame(LEVEL3_COD=shp$LEVEL3_COD, cold=C, hot=H)
DF$PE_spot <- NA
DF$PE_spot[DF$hot==1] <- 1
DF$PE_spot[DF$cold==1] <- 2
shp <- merge(shp, DF, by = "LEVEL3_COD", all = TRUE)
shp$PE_spot <- as.factor(shp$PE_spot)

sdf <- sf::st_drop_geometry(shp) # for quick object access


# add area to define thick lines for small countries
shp$area <- st_area(shp)

min.area <- 1.5e+9
class(min.area)
library(units)
units(min.area) <- as_units("m^2")
thicc_lines <- shp[which(shp$area<min.area),]
thicc_lines <- thicc_lines[!is.na(thicc_lines$PE_spot),]

ggplot(shp)+
  geom_sf(aes(fill=PE_spot),lwd=0) + 
  geom_sf(data=thicc_lines, col="#00BFC4", show.legend=F, lwd=2)+
  scale_fill_discrete(na.value="grey80")+
  theme(panel.border = element_blank())+
  ggtitle("Phylogenetic Endemism Hotspots and Coldspots")
ggsave("figures/PE_hot_cold.png", width=7, height=4, units = "in", dpi = 600, bg = "white")

# Diversity hotspots ##
C <- coldspots(shp$PD_obs_ts) # coldspots
H <- hotspots(shp$PD_obs_ts) # hotspots
DF <- data.frame(LEVEL3_COD=shp$LEVEL3_COD, cold=C, hot=H)
DF$PD_spot <- NA
DF$PD_spot[DF$hot==1] <- 1
DF$PD_spot[DF$cold==1] <- 2
shp <- merge(shp, DF, by = "LEVEL3_COD", all = TRUE)
shp$PD_spot <- as.factor(shp$PD_spot)

thicc_lines <- shp[which(shp$area<min.area),]
thicc_lines <- thicc_lines[!is.na(thicc_lines$PD_spot),]

ggplot(shp)+
  geom_sf(aes(fill=PD_spot),lwd=0) + 
  geom_sf(data=thicc_lines, col="#00BFC4", show.legend=F, lwd=2)+
  scale_fill_discrete(na.value="grey80")+
  theme(panel.border = element_blank())+
  ggtitle("Observed PD Hotspots and Coldspots")
ggsave("figures/PD_hot_cold.png", width=7, height=4, units = "in", dpi = 600, bg = "white")

# Standardized PD hotspots
C <- coldspots(shp$zscore_ts) # coldspots
H <- hotspots(shp$zscore_ts) # hotspots
DF <- data.frame(LEVEL3_COD=shp$LEVEL3_COD, cold=C, hot=H)
DF$PD_zspot <- NA
DF$PD_zspot[DF$hot==1] <- 1
DF$PD_zspot[DF$cold==1] <- 2
shp <- merge(shp, DF, by = "LEVEL3_COD", all = TRUE)
shp$PD_zspot <- as.factor(shp$PD_zspot)

thicc_lines <- shp[which(shp$area<min.area),]
thicc_lines <- thicc_lines[!is.na(thicc_lines$PD_zspot),]

ggplot(shp)+
  geom_sf(aes(fill=PD_zspot),lwd=0) + 
  geom_sf(data=thicc_lines, aes(col=PD_zspot), show.legend=F, lwd=2)+
  scale_fill_discrete(na.value="grey80")+
  theme(panel.border = element_blank())+
  ggtitle("Standardized effect size of PD vs. null communities", subtitle = "Hotspots and Coldspots")
ggsave("figures/SES.PD_hot_cold.png", width=7, height=4, units = "in", dpi = 600, bg = "white")


co <- 0.025 # top 9 countries
orders <- apply(shp[,c("mat_change", "pre_change", "deforest.1", "hfp.1")], 2, order, decreasing=T)
res2 <- matrix(nrow=nrow(orders), ncol=ncol(orders))
for(i in 1:ncol(orders)){
  res2[,i] <- shp$LEVEL3_COD[orders[,i]]
}
colnames(res2) <- colnames(shp[,c("mat_change", "pre_change", "deforest.1", "hfp.1")])

# Top 2.5% in richness, PD observed, PD zscore + PE:
res[1:9, c("richness", "PD_obs_ts", "zscore_ts", "PE")]

# top s.5% in temp change, pre change, deforestation and hfp:
res2[1:9,c("mat_change", "pre_change", "deforest.1", "hfp.1")]

# get a cumulative threat index:
#reshape::rescaler(abs(shp$mat_change))
shp$threat <- abs(reshape::rescaler(shp$mat_change))+
  abs(reshape::rescaler(shp$pre_change))+
  abs(reshape::rescaler(shp$deforest.1))+
  reshape::rescaler(shp$hfp.1)


ggplot(shp) + 
  geom_sf(aes(fill=threat),lwd=0, col=NA) + 
  #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
  #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
  #                        begin = lcol, end = sqrt(ucol))+
  scale_fill_viridis_c("threat", option = "plasma")+ #, 
  theme(legend.position = c(0.18, 0.3),
        legend.key.height = unit(6,"mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        #panel.border = element_blank(),
        text = element_text(size = 10),
        # axis.ticks.x = element_line(color="white"),
        # axis.text.x = element_text(color="white"),
        # plot.margin = margin(0, 0, 0, 0.5, "cm")
  )+
  #coord_sf(expand = F, label_axes = "-N--", ylim=c(-6200000, 8200000))+
  xlab(" ")



plot(shp$mat_change~shp$zscore_ts)
plot(shp$pre_change~shp$zscore_ts)
plot(shp$deforest.1~shp$zscore_ts)
plot(shp$hfp.1~shp$zscore_ts)
cor.test(shp$hfp.1,shp$zscore_ts)

plot(shp$pre_change~shp$PD_obs_ts)
plot(shp$deforest.1~shp$PD_obs_ts)
plot(shp$hfp.1~shp$PD_obs_ts, pch=shp$LEVEL3_COD)

ggplot(shp) + 
  geom_sf(aes(fill=bio12.1),lwd=0, col=NA) + 
  #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
  #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
  #                        begin = lcol, end = sqrt(ucol))+
  scale_fill_viridis_c("MAT", option = "plasma")+ #, 
  theme(legend.position = c(0.18, 0.3),
        legend.key.height = unit(6,"mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        #panel.border = element_blank(),
        text = element_text(size = 10),
        # axis.ticks.x = element_line(color="white"),
        # axis.text.x = element_text(color="white"),
        # plot.margin = margin(0, 0, 0, 0.5, "cm")
  )+
  #coord_sf(expand = F, label_axes = "-N--", ylim=c(-6200000, 8200000))+
  xlab(" ")
