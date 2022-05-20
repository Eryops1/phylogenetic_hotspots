
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




# COMPLETENESS  ########################################

cor.test(shp$richness, shp$SES.PD_ts)
cor.test(shp$richness, shp$PD_obs)

dat <- st_drop_geometry(shp)

# missing data
library(naniar)
library(UpSetR)
vis_miss(dat)

# how many missings?
n_var_miss(dat)
gg_miss_upset(dat, nsets=20)

# get rid of SD variables
dat.sub <- dat[,-grep("_sd", names(dat))]
n_var_miss(dat.sub)
gg_miss_upset(dat.sub, nsets=37)

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

ggplot(dat.sub, aes(y=area, group=is.na(PC_primf_mean)))+
            geom_boxplot(varwidth = T)+
            scale_y_log10()+
            xlab("PRIMF_past change == NA (positive=T)")


dat_no.na <- na.omit(dat.sub)
dim(dat_no.na)

# Correlation ##########################
source("99_functions.R")
library(rstatix)
cmat <- cor_mat(dat_no.na[,-grep("LEVEL3|hfp|deforest|bio|change|PC|FC",names(dat_no.na))], method = "s")
cpmat <- cor_pmat(dat_no.na[,-grep("LEVEL3|hfp|deforest|bio|change|PC|FC",names(dat_no.na))], method = "s")
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

# subset spatial object to dat_no.na variables
shp <- shp[,which(names(shp) %in% names(dat_no.na))]
library(spdep)
# create neighbor list
tmp <- na.omit(shp)
nb <- spdep::poly2nb(tmp, row.names = shp$LEVEL3_COD)
names(nb) <- tmp$LEVEL3_COD

# distance object
col.W <- nb2listw(nb, style="W", zero.policy = TRUE)


# *** Moran's I -----------------------------------------------------------

#moran.test(tmp$richness, col.W, zero.policy = TRUE)
tmpm <- st_drop_geometry(tmp)
tmpm <- tmpm[,-1]
res <- apply(tmpm, 2, moran.test, listw=col.W, zero.policy = TRUE)
res.df <- sapply(res, "[[", "estimate")
res.df <- rbind(res.df, sapply(res, "[[", "p.value"))
row.names(res.df) <- c("moran", "expect", "var", "pvalue")
res.df <- as.data.frame(t(res.df))
res.df$twoSD <- 2*sqrt(res.df$var)
ggplot(res.df, aes(y=expect, x=row.names(res.df), shape=factor(pvalue<0.05)))+
  geom_pointrange(aes(ymin=expect-twoSD, ymax=expect+twoSD))+
  geom_point(aes(y=moran), x=row.names(res.df))+
  xlab("Moran's I")+
  coord_flip()





# *** Lee's spatial correlation -------------------------------------------

#A positive Lee’s L indicates that clusters match for the two variables. A
#negative value indicates that the clusters have an opposite spatial
#distribution . A value around zero indicates that the spatialstructures of the
#two variables do not match. The significance of the values of Lee’s L was
#evaluated using a Monte Carlo test with 999 randomizations.

# use this to check where PD hotspots and threats overlap!

nb <- spdep::poly2nb(tmp, row.names = tmp$LEVEL3_COD)
col.W <- nb2listw(nb, style="W", zero.policy = TRUE)

les <- apply(tmpm[,-3], 2, lee.test, y=tmpm$SES.PD_ts, 
             listw=col.W, zero.policy = TRUE, alternative="two.sided")
les.df <- sapply(les, "[[", "estimate")
les.df <- rbind(les.df, sapply(les, "[[", "p.value"))
row.names(les.df) <- c("Lee", "expect", "var", "pvalue")
les.df <- as.data.frame(t(les.df))
les.df$twoSD <- 2*sqrt(les.df$var)

ggplot(les.df, aes(y=expect, x=row.names(les.df)))+
  geom_pointrange(aes(ymin=expect-twoSD, ymax=expect+twoSD))+
  geom_point(aes(y=Lee, x=row.names(les.df), shape=factor(pvalue<0.05)), col="red", size=3)+
  ylab("Lee's L with SES.PD")+
  scale_shape("p<0.05")+
  xlab("")+
  coord_flip()
ggsave("figures/LeesL.png", width=5, height=8, units = "in", dpi = 300)
#A positive Lee’s L indicates that clusters match for the two variables. A
#negative value indicates that the clusters have an opposite spatial
#distribution . A value around zero indicates that the spatialstructures of the
#two variables do not match. The significance of the values of Lee’s L was
#evaluated using a Monte Carlo test with 999 randomizations.





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


ggplot(data=shp, aes(x=richness, y=PD_obs))+
  geom_point()+
  geom_point(aes(y=pd_rand_mean_ts), col="salmon")+
  scale_x_continuous(trans="log")+
  scale_y_continuous(trans="log")

hist(shp$pd_rand_mean_ts)



# PD_obs: observed PD in community
# pd_rand_mean: mean PD in null communities
# pd_obs_rank: Rank of observed PD vs. null communities
# pd_obs_z: Standardized effect size of PD vs. null communities = (PD_obs - pd_rand_mean) / pd_rand_sd
# pd_obs_p: P-value (quantile) of observed PD vs. null communities = mpd_obs_rank / iter + 1



# Variable importance -----------------------------------------------------

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
n <- names(dat_no.na[,-grep("LEVEL3_COD", names(dat_no.na))])
f <- as.formula(paste("SES.PD_ts ~", paste(n[!n %in% "SES.PD_ts"], collapse = " + ")))

run.gbm <- function(i, seeds, data, trControl=gbmControl,
                    tuneGrid=gbmGrid){
  set.seed(s[i])
  if(!i%%1)cat(i,"\r")
  gbm_temp <- train(f, data, method = "gbm",
                  trControl = gbmControl, verbose=FALSE,
                  tuneGrid = gbmGrid)
  temp <- list(summary(gbm_temp)$var, s[i])
  return(temp)
}

n_cores=4
s <- seq(500,501,1)
system.time(
  PD_list <- mclapply(1:length(s), seeds = s, run.gbm,
                      data = dat_no.na, mc.cores=n_cores)
)
#saveRDS(sr_list, "../processed_data/sr_list100_lgm.rds")
PD_list[1]



# Distribution, transformations, standardization ---------------------

dat <- dat[,-grep("_sd$", names(dat))]
# Distributions
tmp <- apply(dat[,-grep("LEVEL3_COD", names(dat))], 2, shapiro.test)
any(as.numeric(sapply(tmp, "[[", "p.value"))>0.05)
#names(tmp)[which(as.numeric(sapply(tmp, "[[", "p.value"))>0.05)] # nothings normally distributed

library(fitdistrplus)
par(mfrow=c(5,5))
apply(dat[,c(1,3:25)], 2, hist, breaks=20, main="bla") 
apply(dat[,c(26:51)], 2, hist, breaks=20)
apply(dat[,c(52:73)], 2, hist, breaks=20)


# Standardization of PD
## compare ses.pd to pd_obs/richness or area
hist(shp$PD_obs/shp$richness)
hist(shp$PD_obs)
ggplot(shp, aes(x=PD_obs, y=SES.PD_ts))+
  geom_point()
ggplot(shp, aes(x=PD_obs, y=richness))+
  geom_point()
ggplot(shp, aes(x=PD_obs, y=PD_obs/richness))+
  geom_point()

# Maps --------------------------------------------------------------------

# format to sf object for plotting
shp <- st_as_sf(s)

# transform to Behrmann projection
shp <- st_transform(shp, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") # Behrmann

min.area <- 1.5e+9
class(min.area)
library(units)
units(min.area) <- as_units("m^2")
thicc_lines <- shp[which(shp$area<min.area),]
lcol <- min(thicc_lines$PD_obs)/max(shp$PD_obs)
ucol <- max(thicc_lines$PD_obs)/max(shp$PD_obs)
(pd_map <- ggplot(shp) + 
    geom_sf(aes(fill=PD_obs),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=PD_obs), show.legend=F)+
    scale_colour_viridis_c("PD", option = "plasma", trans = "sqrt", 
                            begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("PD", option = "plasma", trans="sqrt")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10))+
    xlab(" ")+
    ggtitle("Raw PD Faith")
)

# +100 tweak.....
lcol <- min(thicc_lines$SES.PD_ts+100)/max(shp$SES.PD_ts+100)
ucol <- max(thicc_lines$SES.PD_ts+100)/max(shp$SES.PD_ts+100)
(pd_ses_map <- ggplot(shp) + 
    geom_sf(aes(fill=SES.PD_ts+100),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=SES.PD_ts+100), show.legend=F)+
    scale_colour_viridis_c("SES.PD", option = "plasma", trans="sqrt",
                            begin = lcol, end = ucol)+
    scale_fill_viridis_c("SES.PD+100", option = "plasma",trans="sqrt")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")+
    ggtitle("SES PD Faith")
)
# (pd_rand_map <- ggplot(shp) + 
#     geom_sf(aes(fill=pd_rand_mean_ts),lwd=0, col=NA) + 
#     #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
#     #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
#     #                        begin = lcol, end = sqrt(ucol))+
#     scale_fill_viridis_c("PD null distribution", option = "plasma")+ #, 
#     theme(legend.position = c(0.18, 0.3),
#           legend.key.height = unit(6,"mm"),
#           legend.background = element_blank(),
#           legend.key = element_blank(),
#           panel.background = element_blank(),
#           text = element_text(size = 10),
#     )+
#     xlab("")
#)


lcol <- min(thicc_lines$PE)/max(shp$PE)
ucol <- max(thicc_lines$PE)/max(shp$PE)
(pe_map <- ggplot(shp) + 
    geom_sf(aes(fill=PE),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=PE), show.legend=F)+
    scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
                            begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("PE", option = "plasma", trans="sqrt")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )+
    xlab(" ")+
    ggtitle("Phylogenetic endemism")
)

lcol <- min(thicc_lines$richness)/max(shp$richness)
ucol <- max(thicc_lines$richness)/max(shp$richness)
(sr_map <- ggplot(shp) + 
    geom_sf(aes(fill=richness),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=richness), show.legend=F)+
    scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
                            begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("SR", option = "plasma", trans = "sqrt")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )+
    xlab(" ")
)

plot_grid(pd_map, pd_ses_map, pe_map, sr_map, nrow = 4)
ggsave("figures/maps.png", width=8, height=16, units = "in", dpi = 600, bg = "white")

plot(shp$PD_obs_ts, shp$pd_rand_mean_ts)
plot(shp$PD_obs_ts, shp$PE)

# lm1 <- lm(data=shp, PD_obs_ts~mrd + soil + mat_mean + tra_mean + pre_mean + prs_mean + elev_range+
#      tri + area + hfp.1 + deforest.1 + bio1.1 + bio12.1)




# Hotspot definition ------------------------------------------------------

s@data$LEVEL3_COD <- s@data$LEVEL_3_CO
shp <- merge(shp, s@data[,c("LEVEL3_COD", "LEVEL_NAME")])

# Endemism hotspots ##
C <- coldspots(shp$PE) # coldspots
H <- hotspots(shp$PE) # hotspots

## Merge endemism values to shapefile of grid cells.
DF <- data.frame(LEVEL3_COD=shp$LEVEL3_COD, cold=C, hot=H)
DF$PE_spot <- 0
DF$PE_spot[DF$hot==1] <- 1
DF$PE_spot[DF$cold==1] <- -1
shp <- merge(shp, DF[,c("LEVEL3_COD", "PE_spot")], by = "LEVEL3_COD", all = TRUE)
shp$PE_spot


# add area to define thick lines for small countries
shp$area <- st_area(shp)

min.area <- 1.5e+9
class(min.area)
library(units)
units(min.area) <- as_units("m^2")
thicc_lines <- shp[which(shp$area<min.area),]
thicc_lines <- thicc_lines[thicc_lines$PE_spot!=0,]

# 
# (sr_map2 <- ggplot(test) + 
#     geom_sf(aes(fill=sr),lwd=0, col=NA) + 
#     geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
#     scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
#                            begin = lcol, end = sqrt(ucol))+
#     scale_fill_viridis_c("SR", option = "plasma", trans = "sqrt")+ #, 

ggplot(shp)+
  geom_sf(aes(fill=factor(PE_spot)),lwd=0, col=NA) + 
  geom_sf(data=thicc_lines, aes(col=factor(PE_spot)), show.legend=F, lwd=2)+
  scale_fill_manual(values = c("blue", "grey80", "red"))+
  scale_color_manual(values = c("blue", "red"))+
  theme(panel.border = element_blank())+
  ggtitle("Phylogenetic Endemism Hotspots and Coldspots")
ggsave("figures/PE_hot_cold.png", width=7, height=4, units = "in", dpi = 600, bg = "white")

# Names islands:
View(st_drop_geometry(thicc_lines[,c("PE_spot.x", "LEVEL3_COD", "LEVEL_NAME")]))
thicc_lines$LEVEL_NAME


# Diversity hotspots ##
C <- coldspots(shp$PD_obs) # coldspots
H <- hotspots(shp$PD_obs) # hotspots
DF <- data.frame(LEVEL3_COD=shp$LEVEL3_COD, cold=C, hot=H)
DF$PD_spot <- 0
DF$PD_spot[DF$hot==1] <- 1
DF$PD_spot[DF$cold==1] <- -1
shp <- merge(shp, DF[,c("LEVEL3_COD", "PD_spot")], by = "LEVEL3_COD", all = TRUE)


thicc_lines <- shp[which(shp$area<min.area),]
thicc_lines <- thicc_lines[thicc_lines$PD_spot!=0,]

ggplot(shp)+
  geom_sf(aes(fill=factor(PD_spot)),lwd=0, col=NA) + 
  geom_sf(data=thicc_lines, aes(col=factor(PD_spot)), show.legend=F, lwd=2)+
  scale_fill_manual(values = c("blue", "grey80", "red"))+
  scale_color_manual(values = c("blue", "grey80", "red"))+
  theme(panel.border = element_blank())+
  ggtitle("Observed PD Hotspots and Coldspots")
ggsave("figures/PD_hot_cold.png", width=7, height=4, units = "in", dpi = 600, bg = "white")

# Names islands:
View(st_drop_geometry(thicc_lines[,c("PD_spot.x", "LEVEL3_COD", "LEVEL_NAME")]))
thicc_lines$LEVEL_NAME


# Standardized PD hotspots
C <- coldspots(shp$SES.PD_ts) # coldspots
H <- hotspots(shp$SES.PD_ts) # hotspots
DF <- data.frame(LEVEL3_COD=shp$LEVEL3_COD, cold=C, hot=H)
DF$SES.PD_spot <- 0
DF$SES.PD_spot[DF$hot==1] <- 1
DF$SES.PD_spot[DF$cold==1] <- -1
shp <- merge(shp, DF[,c("LEVEL3_COD", "SES.PD_spot")], by = "LEVEL3_COD", all = TRUE)


thicc_lines <- shp[which(shp$area<min.area),]
thicc_lines <- thicc_lines[thicc_lines$SES.PD_spot!=0,]

ggplot(shp)+
  geom_sf(aes(fill=factor(SES.PD_spot)),lwd=0, col=NA) + 
  geom_sf(data=thicc_lines, aes(col=factor(SES.PD_spot)), show.legend=F, lwd=2)+
  scale_fill_manual("spot", values = c("blue", "grey80", "red"))+
  scale_color_manual(values = c("red"))+ # MANUALLY CHECK if all are positive
  theme(panel.border = element_blank())+
  ggtitle("SES.PD Hotspots and Coldspots")
ggsave("figures/SES.PD_hot_cold.png", width=7, height=4, units = "in", dpi = 600, bg = "white")

# Names islands:
View(st_drop_geometry(thicc_lines[,c("SES.PD_spot", "LEVEL3_COD", "LEVEL_NAME")]))
thicc_lines$LEVEL_NAME

# 
# co <- 0.025 # top 9 countries
# orders <- apply(shp[,c("mat_change", "pre_change", "deforest.1", "hfp.1")], 2, order, decreasing=T)
# res2 <- matrix(nrow=nrow(orders), ncol=ncol(orders))
# for(i in 1:ncol(orders)){
#   res2[,i] <- shp$LEVEL3_COD[orders[,i]]
# }
# colnames(res2) <- colnames(shp[,c("mat_change", "pre_change", "deforest.1", "hfp.1")])
# 
# # Top 2.5% in richness, PD observed, PD zscore + PE:
# res[1:9, c("richness", "PD_obs_ts", "zscore_ts", "PE")]
# 
# # top s.5% in temp change, pre change, deforestation and hfp:
# res2[1:9,c("mat_change", "pre_change", "deforest.1", "hfp.1")]
# 
# # get a cumulative threat index:
# #reshape::rescaler(abs(shp$mat_change))
# shp$threat <- abs(reshape::rescaler(shp$mat_change))+
#   abs(reshape::rescaler(shp$pre_change))+
#   abs(reshape::rescaler(shp$deforest.1))+
#   reshape::rescaler(shp$hfp.1)
# 
# 
# ggplot(shp) + 
#   geom_sf(aes(fill=threat),lwd=0, col=NA) + 
#   #geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
#   #scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
#   #                        begin = lcol, end = sqrt(ucol))+
#   scale_fill_viridis_c("threat", option = "plasma")+ #, 
#   theme(legend.position = c(0.18, 0.3),
#         legend.key.height = unit(6,"mm"),
#         legend.background = element_blank(),
#         legend.key = element_blank(),
#         panel.background = element_blank(),
#         #panel.border = element_blank(),
#         text = element_text(size = 10),
#         # axis.ticks.x = element_line(color="white"),
#         # axis.text.x = element_text(color="white"),
#         # plot.margin = margin(0, 0, 0, 0.5, "cm")
#   )+
#   #coord_sf(expand = F, label_axes = "-N--", ylim=c(-6200000, 8200000))+
#   xlab(" ")


# Diversity hotspots ##
C <- coldspots(shp$richness) # coldspots
H <- hotspots(shp$richness) # hotspots
DF <- data.frame(LEVEL3_COD=shp$LEVEL3_COD, cold=C, hot=H)
DF$sr_spot <- 0
DF$sr_spot[DF$hot==1] <- 1
DF$sr_spot[DF$cold==1] <- -1
shp <- merge(shp, DF[,c("LEVEL3_COD", "sr_spot")], by = "LEVEL3_COD", all = TRUE)


thicc_lines <- shp[which(shp$area<min.area),]
thicc_lines <- thicc_lines[thicc_lines$sr_spot!=0,]

ggplot(shp)+
  geom_sf(aes(fill=factor(sr_spot)),lwd=0, col=NA) + 
  geom_sf(data=thicc_lines, aes(col=factor(sr_spot)), show.legend=F, lwd=2)+
  scale_fill_manual(values = c("blue", "grey80", "red"))+
  scale_color_manual(values = c("blue", "grey80", "red"))+
  theme(panel.border = element_blank())+
  ggtitle("Observed SR Hotspots and Coldspots")
ggsave("figures/SR_hot_cold.png", width=7, height=4, units = "in", dpi = 600, bg = "white")

View(st_drop_geometry(thicc_lines[,c("sr_spot", "LEVEL3_COD", "LEVEL_NAME")]))
thicc_lines$LEVEL_NAME



# High threat regions -----------------------------------------------------


nb <- spdep::poly2nb(shp, row.names = shp$LEVEL3_COD)
col.W <- nb2listw(nb, style="W", zero.policy = TRUE)
thr <- shp[,grep("LEVEL|SES|_spot|zspot|PE|hfp|deforest|_change|FC|PC", names(shp))]
thr <- st_drop_geometry(thr)

tes <- apply(thr[,grep("hfp|deforest|_change|FC|PC", names(thr))], 2, lee.test, y=thr$SES.PD_ts, 
             listw=col.W, zero.policy = TRUE, alternative="two.sided", na.action=na.omit)
tes.df <- sapply(tes, "[[", "estimate")
tes.df <- rbind(tes.df, sapply(tes, "[[", "p.value"))
row.names(tes.df) <- c("Lee", "expect", "var", "pvalue")
tes.df <- as.data.frame(t(tes.df))
tes.df$twoSD <- 2*sqrt(tes.df$var)

ggplot(tes.df, aes(y=expect, x=gsub("_mean", "", row.names(tes.df))))+
  geom_pointrange(aes(ymin=expect-twoSD, ymax=expect+twoSD))+
  geom_point(aes(y=Lee, x=gsub("_mean", "", row.names(tes.df)), shape=factor(pvalue<0.05)), col="red", size=3)+
  ylab("Lee's L with SES.PD")+
  scale_shape("p<0.05")+
  xlab("")+
  coord_flip()
ggsave("figures/LeesL_PD_threat.png", width=5, height=4, units = "in", dpi = 300)
#A positive Lee’s L indicates that clusters match for the two variables. A
#negative value indicates that the clusters have an opposite spatial
#distribution . A value around zero indicates that the spatialstructures of the
#two variables do not match. The significance of the values of Lee’s L was
#evaluated using a Monte Carlo test with 999 randomizations.

tes <- apply(thr[,grep("hfp|deforest|_change|FC|PC", names(thr))], 2, lee.test, y=thr$SES.PD_spot, 
             listw=col.W, zero.policy = TRUE, alternative="two.sided", na.action=na.omit)
tes.df <- sapply(tes, "[[", "estimate")
tes.df <- rbind(tes.df, sapply(tes, "[[", "p.value"))
row.names(tes.df) <- c("Lee", "expect", "var", "pvalue")
tes.df <- as.data.frame(t(tes.df))
tes.df$twoSD <- 2*sqrt(tes.df$var)

ggplot(tes.df, aes(y=expect, x=gsub("_mean", "", row.names(tes.df))))+
  geom_pointrange(aes(ymin=expect-twoSD, ymax=expect+twoSD))+
  geom_point(aes(y=Lee, x=gsub("_mean", "", row.names(tes.df)), shape=factor(pvalue<0.05)), col="red", size=3)+
  ylab("Lee's L with SES.PD hotspots")+
  scale_shape("p<0.05")+
  xlab("")+
  coord_flip()
ggsave("figures/LeesL_SES.PD_hotspots_threat.png", width=5, height=4, units = "in", dpi = 300)


# Threat hotspots ##
C <- coldspots(shp$hfp_mean) # coldspots
H <- hotspots(shp$hfp_mean) # hotspots
DF <- data.frame(LEVEL3_COD=shp$LEVEL3_COD, cold=C, hot=H)
DF$hfp_spot <- 0
DF$hfp_spot[DF$hot==1] <- 1
DF$hfp_spot[DF$cold==1] <- -1
shp <- merge(shp, DF[,c("LEVEL3_COD", "hfp_spot")], by = "LEVEL3_COD", all = TRUE)


thicc_lines <- shp[which(shp$area<min.area),]
thicc_lines <- thicc_lines[thicc_lines$hfp_spot!=0,]

ggplot(shp)+
  geom_sf(aes(fill=factor(hfp_spot)),lwd=0, col=NA) + 
  geom_sf(data=thicc_lines, aes(col=factor(hfp_spot)), show.legend=F, lwd=2)+
  scale_fill_manual(values = c("blue", "grey80", "red"))+
  scale_color_manual(values = c("blue", "grey80", "red"))+
  theme(panel.border = element_blank())+
  ggtitle("Observed HFP Hotspots and Coldspots")
ggsave("figures/PD_hot_cold.png", width=7, height=4, units = "in", dpi = 600, bg = "white")


