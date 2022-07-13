
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
gc()
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
shp <- shp[,-grep("obs_p|obs_rank|reps|LEVEL2|LEVEL1|LEVEL_3_CO|ID|\\.3|_rw|CONTI|REGION", names(shp))]
names(shp)<- gsub("\\.1", "_mean", names(shp))
names(shp)<- gsub("\\.2", "_sd", names(shp))

# remove BOU that has no data
shp <- shp[!shp$LEVEL3_COD=="BOU",]


# Explorative plots -------------------------------------------------------

# plot(st_drop_geometry(shp)[,c(1,4:12,13:23)])
# plot(st_drop_geometry(shp)[,c(1,4:12,24:35)])
# plot(st_drop_geometry(shp)[,c(1,4:12,36:57)])




# COMPLETENESS  ########################################

cor.test(shp$richness, shp$SES.PD)
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
gg_miss_upset(dat.sub, nsets=40)

# remove SES.PD NAs
dat.sub <- dat.sub[-which(is.na(dat.sub$SES.PD)),]
n_var_miss(dat.sub)
gg_miss_upset(dat.sub, nsets=37)


plot_grid(
ggplot(dat.sub, aes(y=area, group=is.na(hfp_mean)))+
  geom_boxplot(varwidth = T)+
  scale_y_log10()+
  xlab("HFP == NA (positive=T: n=27)"),
ggplot(dat.sub, aes(y=PD_obs, group=is.na(hfp_mean)))+
  geom_boxplot(varwidth = T)+
  xlab("HFP == NA (positive=T: n=27)"),
ggplot(dat.sub, aes(y=SES.PD, group=is.na(hfp_mean)))+
  geom_boxplot(varwidth = T)+
  xlab("HFP == NA (positive=T: n=27)"),
ggplot(dat.sub, aes(y=SES.PE, group=is.na(hfp_mean)))+
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

# 316 bot countries got data!

# Correlation ##########################
source("99_functions.R")
library(rstatix)
cmat <- cor_mat(dat_no.na[,-grep("LEVEL|hfp|deforest|bio|change|PC|FC",names(dat_no.na))], method = "s")
cpmat <- cor_pmat(dat_no.na[,-grep("LEVEL|hfp|deforest|bio|change|PC|FC",names(dat_no.na))], method = "s")
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
#row.names(cmat) <- fn


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
tmpm <- tmpm[,-c(1,2)] # exclude level names
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

les <- apply(tmpm[,-3], 2, lee.test, y=tmpm$SES.PD, 
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

temp <- dat_no.na[,-grep("LEVEL|hfp|deforest|bio|change|pd_rand|PD_obs|WE|PE",
                         names(dat_no.na))]
lm_pd <- lm(SES.PD ~ . , data=temp)
summary(lm_pd)
sort(car::vif(lm_pd))
# mat, subtropmbf are problematic

# test changes
lm_pd1 <- lm(SES.PD ~ . -sub_trop_mbf, data=temp)
sort(car::vif(lm_pd1))
lm_pd2 <- lm(SES.PD ~ . -mat_m, data=temp)
sort(car::vif(lm_pd2))
lm_pd3 <- lm(SES.PD ~ . -mat_m -sub_trop_mbf, data=temp)
sort(car::vif(lm_pd3))
lm_pd4 <- lm(SES.PD ~ . -mat_m -sub_trop_mbf -pre_lgm_ano_m, data=temp)
sort(car::vif(lm_pd4))





plot_grid(ncol=2,
          
          ggplot(data=shp, aes(fill=PD_obs))+
            geom_sf(lwd=0)
          ,
          ggplot(data=shp, aes(fill=SES.PD))+
            geom_sf(lwd=0)
          ,
          ggplot(data=shp, aes(x=richness, y=SES.PD))+
            geom_text(aes(label=LEVEL3_COD))
          ,
          ggplot(data=shp, aes(x=pd_rand_mean, y=PD_obs, col=log(richness)))+
            geom_point()+
            geom_abline(x=1)+
            scale_x_continuous(trans="sqrt")+
            scale_y_continuous(trans="sqrt")
          # observed PD is usually lower than for a random species sample for the phylogey that matches the SR
)


ggplot(data=shp, aes(x=richness, y=PD_obs))+
  geom_point()+
  geom_point(aes(y=pd_rand_mean), col="salmon")+
  scale_x_continuous(trans="log")+
  scale_y_continuous(trans="log")

hist(shp$pd_rand_mean)



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
n <- names(dat_no.na[,-grep("LEVEL", names(dat_no.na))])
f <- as.formula(paste("SES.PD ~", paste(n[!n %in% "SES.PD"], collapse = " + ")))

run.gbm <- function(i, seeds, data, trControl=gbmControl,
                    tuneGrid=gbmGrid, f=f){
  set.seed(s[i])
  if(!i%%1)cat(i,"\r")
  temp <- train(f, data, method = "gbm",
                  trControl = gbmControl, verbose=FALSE,
                  tuneGrid = gbmGrid)
  return(temp)
}

n_cores=4
s <- seq(500,507,1) # models are rather slow with that many variables, ca. 5 minutes per model
system.time(
  PD_list <- mclapply(i=1:length(s), seeds = s, run.gbm,
                      data = dat_no.na, mc.cores=n_cores, f=f)
)

# this throws errors currently, no fucking idea why.... single rungs outside parallel are fine

#pd_list <- unlist(PD_list, recursive = F)
pdl <- lapply(PD_list, varImp)
pdm <- as.data.frame(sapply(pdl, "[", "importance"))
su <- rowSums(pdm)
pdm$pred = row.names(pdm) 
pddf <- reshape2::melt(pdm, id.name= "pred", value.name = "varImp")
# sort factors
pddf$pred <- factor(pddf$pred)
pddf$pred <- factor(pddf$pred, levels=names(sort(su, decreasing = F)))

ggplot(pddf, aes(x=pred, y=varImp))+
  geom_boxplot()+
  ylab("Variable importance")+
#  facet_wrap(~Var2, scales = "free")+
#  theme(axis.text.x = element_text(size=6, angle = 45, hjust=1),
#        strip.background = element_blank())+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size=9))+
  coord_flip()
ggsave("figures/varImp_SES_PD_gbm10_runs.png", width=7, height=10, units = "in", dpi = 600)



## formula
f2 <- as.formula(paste("SES.PE ~", paste(n[!n %in% "SES.PE"], collapse = " + ")))

n_cores=4
s <- seq(500,507,1) # models are rather slow with that many variables, ca. 5 minutes per model
system.time(
  PE_list <- mclapply(1:length(s), seeds = s, run.gbm, f=f2,
                      data = dat_no.na, mc.cores=n_cores)
)

pel <- lapply(PE_list, varImp)
pem <- as.data.frame(sapply(pel, "[", "importance"))
su <- rowSums(pem)
pem$pred = row.names(pem) 
pedf <- reshape2::melt(pem, id.name= "pred", value.name = "varImp")
# sort factors
pedf$pred <- factor(pedf$pred)
pedf$pred <- factor(pedf$pred, levels=names(sort(su, decreasing = F)))

ggplot(pedf, aes(x=pred, y=varImp))+
  geom_boxplot()+
  ylab("Variable importance")+
  #  facet_wrap(~Var2, scales = "free")+
  #  theme(axis.text.x = element_text(size=6, angle = 45, hjust=1),
  #        strip.background = element_blank())+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size=9))+
  coord_flip()
ggsave("figures/varImp_SES_PE_gbm10_runs.png", width=7, height=14, units = "in", dpi = 600)





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
ggplot(shp, aes(x=PD_obs, y=richness))+
  geom_point()
cor.test(shp$PD_obs, shp$richness, method="s") # strong correlation!!

ggplot(shp, aes(x=PD_obs/richness, y=richness))+
  geom_point()
cor.test(shp$PD_obs/shp$richness, shp$richness, method="s") # strong correlation!!

ggplot(shp, aes(x=PD_obs, y=area))+
  geom_point()
cor.test(shp$PD_obs, shp$area, method="s") # 0.6

ggplot(shp, aes(x=SES.PD, y=area))+
  geom_point()
cor.test(shp$SES.PD, shp$area, method="s") # -.7
# SES.PD is stronger correlated with area than PD_obs is

# all standardized for area
cor.test(shp$SES.PD/shp$area, shp$richness/shp$area, method="s")


saveRDS(shp, "for_plotting.rds")







# Maps --------------------------------------------------------------------

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
gc()
library(phyloregion) # PD calculations
library(raster)
library(sf)
library(ggplot2)
theme_set(theme_bw())
library(cowplot)
library(rgdal)
library(beepr)
if(!dir.exists("figures"))dir.create("figures")


# format to sf object for plotting
shp <- readRDS("for_plotting.rds")
sfixed <- read_sf("../DATA/shapefile_bot_countries/level3_fixed.gpkg")
sfixed <- sfixed[!sfixed$LEVEL_3_CO=="BOU",]
# move data over to fixed shapefile, can be removed later when fix is completely
# established in earlier pipeline
shp <- merge(sfixed, st_drop_geometry(shp))


# transform to Behrmann projection
shp <- st_transform(shp, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") # Behrmann

min.area <- 1.5e+9
class(min.area)
# library(units)
# units(min.area) <- as_units("m^2")
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


(pd_ses_map <- ggplot(na.omit(shp)) + 
    geom_sf(aes(fill=SES.PD),lwd=0, col=NA) + 
    geom_sf(data=na.omit(thicc_lines), lwd=1.5, aes(col=SES.PD), show.legend=F)+
    scale_colour_viridis_c("SES.PD", option = "plasma")+
    scale_fill_viridis_c("SES.PD", option = "plasma")+ #, 
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

# Just dividing PD by richness instead:
# lcol <- min(thicc_lines$PD_obs/thicc_lines$richness, na.rm = T)/max(thicc_lines$PD_obs/thicc_lines$richness, na.rm = T)
# ucol <- max(thicc_lines$PD_obs/thicc_lines$richness, na.rm = T)/max(thicc_lines$PD_obs/thicc_lines$richness, na.rm = T)
shp2 <- shp[!is.na(shp$SES.PD),]
thicc_lines <- shp2[which(shp2$area<min.area),]
lcol <- min(thicc_lines$PD_obs/thicc_lines$richness)/max(shp2$PD_obs/shp2$richness)
ucol <- max(thicc_lines$PD_obs/thicc_lines$richness)/max(shp2$PD_obs/shp2$richness)
(pd_richness_map <- ggplot(shp2) + 
    geom_sf(aes(fill=PD_obs/richness),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=thicc_lines$PD_obs/thicc_lines$richness), show.legend=F)+
    scale_colour_viridis_c("PD/SR", option = "plasma", trans="sqrt",
                           begin = lcol, end = ucol)+
    scale_fill_viridis_c("PD/SR", option = "plasma",trans="sqrt")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")+
    ggtitle("PD Faith / richness")
  )


thicc_lines <- shp2[which(shp2$area<min.area),]
lcol <- min(thicc_lines$PE_obs)/max(shp2$PE_obs)
ucol <- max(thicc_lines$PE_obs)/max(shp2$PE_obs)
(pe_map <- ggplot(shp2) + 
    geom_sf(aes(fill=PE_obs),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=thicc_lines$PE_obs), show.legend=F)+
    scale_colour_viridis_c("PE", option = "plasma", trans="sqrt",
                           begin = lcol, end = ucol)+
    scale_fill_viridis_c("PE", option = "plasma",trans="sqrt")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")+
    ggtitle("PE")
)

thicc_lines <- shp2[which(shp2$area<min.area),]
# negative numbers fix. range of number spans negative to positive, so the color range needs to be adjusted
lcol <- min(thicc_lines$SES.PE+abs(min(shp2$SES.PE)))/diff(range(shp2$SES.PE)) 
ucol <- max(thicc_lines$SES.PE+abs(min(shp2$SES.PE)))/diff(range(shp2$SES.PE))
(ses.pe_map <- ggplot(shp2) + 
    geom_sf(aes(fill=SES.PE),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=thicc_lines$SES.PE), show.legend=F)+
    scale_colour_viridis_c("PE", option = "plasma",
                           begin = lcol, end = ucol)+
    scale_fill_viridis_c("SES.PE", option = "plasma")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")+
    ggtitle("SES.PE")
)

# Simple SR
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

plot_grid(pd_map, pd_ses_map, pe_map, ses.pe_map, sr_map, nrow = 3)
ggsave("figures/maps.png", width=10, height=10, units = "in", dpi = 600, bg = "white")


# Standardization effects ------------------------------------------------
normalized = function(x)(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
shp2$SES.PD.norm <- normalized(shp2$SES.PD)
shp2$PD_obs.norm <- normalized(shp2$PD_obs)

# get plots and linear model
par(mfrow=c(2,2))
hist(shp2$PD_obs)
hist(shp2$SES.PD)
plot(shp2$PD_obs.norm, shp2$SES.PD.norm)
abline(a=1, b=-1)
#plot(order(shp2$PD_obs), order(shp2$SES.PD))
cor.test(shp2$PD_obs.norm, shp2$SES.PD.norm, method = "s")

df <- st_drop_geometry(shp2)
df$PD_rank <- rank(df$PD_obs)
df$SES.PD_rank <- rank(df$SES.PD)
df$PE_rank <- rank(df$PE_obs)
df$SES.PE_rank <- rank(df$SES.PE)

library(tidyr); library(ggbump); library(dplyr)
df <- df[,c("LEVEL_NAME", "PD_rank", "SES.PD_rank", "PE_rank", "SES.PE_rank", "richness", "area")]
df2 <- pivot_longer(df, cols = c("PD_rank", "SES.PD_rank", "PE_rank", "SES.PE_rank"))

alevel = 0.4
sm <- 4
bonbon_SR <-
  ggplot(df2[grep("PD", df2$name),], aes(name, value, color = richness, group=LEVEL_NAME))+
  geom_bump(size = 1, position = "identity", smooth=sm) +
  scale_color_viridis_c("Species richness", trans="sqrt", alpha = alevel, option = "plasma", guide = "none")+
  scale_x_discrete("", expand = c(0.05,0.05))+
  scale_y_discrete("rank")+
#  geom_point(size = 0, aes(fill=richness)) +
#  scale_fill_viridis_c("Species richness", trans="sqrt", alpha = 1, option = "plasma")+
  theme_minimal()


# bonbon_area <- ggplot(df2[grep("PD", df2$name),], aes(name, value, color = area, group=LEVEL_NAME))+
#   geom_bump(size = 1, position = "identity", smooth=sm) +
#   scale_color_viridis_c(trans="sqrt", alpha = alevel, option = "plasma", guide = "none")+
#   scale_x_discrete("", expand = c(0.05,0.05))+
#   scale_y_discrete("rank")+
#   geom_point(size = 0, aes(fill=area)) +
#   scale_fill_viridis_c("Area", trans="sqrt", alpha = 1, option = "plasma")+
#   theme_minimal()

bonbon_SR_PE <- ggplot(df2[grep("PE", df2$name),], aes(name, value, color = richness, group=LEVEL_NAME))+
  geom_bump(size = 1, position = "identity", smooth=sm) +
  scale_color_viridis_c("Species richness", trans="sqrt", alpha = alevel, option = "plasma", guide = "none")+
  scale_x_discrete("", expand = c(0.05,0.05))+
  scale_y_discrete("rank")+
  geom_point(size = 0, aes(fill=richness)) +
  scale_fill_viridis_c("Species richness", trans="sqrt", alpha = 1, option = "plasma")+
  theme_minimal()

# bonbon_area_PE <- ggplot(df2[grep("PE", df2$name),], aes(name, value, color = area, group=LEVEL_NAME))+
#   geom_bump(size = 2, position = "identity", smooth=sm) +
#   scale_color_viridis_c(trans="sqrt", alpha = alevel, option = "plasma", guide = "none")+
#   scale_x_discrete("", expand = c(0.05,0.05))+
#   scale_y_discrete("rank")+
#   geom_point(size = 0, aes(fill=area)) +
#   scale_fill_viridis_c("Area", trans="sqrt", alpha = 1, option = "plasma")+
#   theme_minimal()

plot_grid(bonbon_SR, bonbon_SR_PE, ncol=2, rel_widths = c(.41,.59))  
ggsave("figures/bonbon_plots.png", width=8, height=4, units = "in", dpi = 600, bg = "white")

#plot_grid(bonbon_SR+coord_flip(), bonbon_SR_PE+coord_flip(), ncol=1)  


       
# SES.PD vs PD_obs colored for SR
ggplot(shp2, aes(x=PD_obs, y=SES.PD, color=richness))+
  geom_point(aes(size=3))+
  scale_color_viridis_c(trans="sqrt")

source("99_functions.R")
shp2$SES.PE.norm <- normalized(shp2$SES.PE)
shp2$PE_obs.norm <- normalized(shp2$PE_obs)

thicc_lines <- shp2[which(shp2$area<min.area),]
lcol <- min(thicc_lines$SES.PE.norm)/max(shp2$SES.PE.norm)
ucol <- max(thicc_lines$SES.PE.norm)/max(shp2$SES.PE.norm)
(pe_standardization_effects_map <- ggplot(shp2) + 
    geom_sf(aes(fill=PE_obs.norm-SES.PE.norm)) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=PE_obs.norm-SES.PE.norm), show.legend=F)+
    scale_color_gradient2("stand. differences",low = "red", mid = "white", high = "blue")+
    scale_fill_gradient2("stand. differences",low = "red", mid = "white", high = "blue")+
    ggtitle("PE[0,1] - SES.PE[0,1]. Blue=SES is smaller, red=SES is bigger than raw values")
)
plot_grid(pd_standardization_effects_map, pe_standardization_effects_map, nrow=2)



# lm1 <- lm(data=shp, PD_obs~mrd + soil + mat_mean + tra_mean + pre_mean + prs_mean + elev_range+
#      tri + area + hfp.1 + deforest.1 + bio1.1 + bio12.1)




# Choropleth maps -------------------------------------------------
## get overlap of PD and PE, and both with SR
library(biscale)

shp2 <- bi_class(shp2, x = PD_obs, y = PE_obs, style = "jenks", dim = 3)
thicc_lines <- shp2[which(shp2$area<min.area),]

# test colors
# "Bluegill", "BlueGold", "BlueOr", "BlueYl", "Brown"/"Brown2", "DkBlue"/"DkBlue2", "DkCyan"/"DkCyan2", "DkViolet"/"DkViolet2", "GrPink"/"GrPink2", "PinkGrn", "PurpleGrn", or "PurpleOr"
# bi_pal(pal = "BlueGold", dim = 3, preview = TRUE)

(chloropl1 <- ggplot() +
  geom_sf(na.omit(shp2), mapping = aes(fill = bi_class), color = NA, size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal = "BlueGold", dim = 3, na.value="white") +
  geom_sf(data=thicc_lines, lwd=1.5, aes(col=bi_class), show.legend=F)+
  bi_scale_color(pal = "BlueGold", dim = 3, na.value="white") +
  bi_theme()+
  #geom_sf(data=m[m$Type=="hotspot area",], fill="green", color=NA, alpha=0.5, show.legend=F)+
  #geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="green", color=NA, alpha=0.5, show.legend=F)+
  ggtitle("PD and PE")+
  theme(title = element_text(size=10)))
legend <- bi_legend(pal = "BlueGold",
                    dim = 3,
                    xlab = "PD",
                    ylab = "PE ",
                    size = 8)
(pd_pe_map <- ggdraw() +
  draw_plot(chloropl1, 0, 0, 1, 1) +
  draw_plot(legend, 0.05, 0.25, 0.2, 0.2))

data <- bi_class(shp2, x = SES.PD, y = SES.PE, style = "jenks", dim = 3)
data$bi_class_ses <- data$bi_class
shp2$bi_class_ses <- data$bi_class_ses
rm(data)

thicc_lines <- data[which(shp2$area<min.area),]
(chloropl2 <- ggplot() +
    geom_sf(shp2, mapping = aes(fill = bi_class_ses), color = NA, size = 0.1, show.legend = FALSE) +
    bi_scale_fill(pal = "BlueGold", dim = 3, na.value="white") +
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=bi_class), show.legend=F)+
    bi_scale_color(pal = "BlueGold", dim = 3, na.value="white") +
    bi_theme()+
    #geom_sf(data=m[m$Type=="hotspot area",], fill="green", color=NA, alpha=0.5, show.legend=F)+
    #geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="green", color=NA, alpha=0.5, show.legend=F)+
    ggtitle("SES.PD and SES.PE")+
  theme(title = element_text(size=10)))
legend <- bi_legend(pal = "BlueGold",
                    dim = 3,
                    xlab = "SES.PD",
                    ylab = "SES.PE ",
                    size = 8)

(ses_pd_pe_map <- ggdraw() +
  draw_plot(chloropl2, 0, 0, 1, 1) +
  draw_plot(legend, 0.05, 0.25, 0.2, 0.2))

plot_grid(pd_pe_map, ses_pd_pe_map, ncol=1)
ggsave("figures/phylo_colorpleth_map2.png", width=10, height=10, units = "in", dpi = 600, bg = "white")


## with Myer hotspots -------------------------------------
m <- st_read("hotspots_fixed.gpkg")

# extract problem polygons
wrld_wrap <- st_wrap_dateline(m, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
m <- st_transform(wrld_wrap, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
plot(m)

#library(ggpattern)
(chloropl3 <- ggplot() +
    geom_sf(shp2, mapping = aes(fill = bi_class_ses), color = NA, size = 0.1, show.legend = FALSE) +
    bi_scale_fill(pal = "BlueGold", dim = 3, na.value="white") +
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=bi_class), show.legend=F)+
    bi_scale_color(pal = "BlueGold", dim = 3, na.value="white") +
    bi_theme()+
    geom_sf(data=m[m$Type=="hotspot area",], fill="red", color=NA, alpha=0.5, show.legend=F)+
    geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="red", color=NA, alpha=0.5, show.legend=F)+
    ggtitle("SES.PD and SES.PE")+
  theme(title = element_text(size=10)))
legend <- bi_legend(pal = "BlueGold", dim = 3, xlab = "SES.PD", ylab = "SES.PE ", size = 8)

(ses_pd_pe_hotspots_map <- ggdraw() +
  draw_plot(chloropl3, 0, 0, 1, 1) +
  draw_plot(legend, 0.05, 0.25, 0.2, 0.2))

plot_grid(ses_pd_pe_hotspots_map)
ggsave("figures/phylo_colorpleth_map2_hotspots.png", width=10, height=5, units = "in", dpi = 600, bg = "white")
ggsave("figures/phylo_colorpleth_map2_hotspots.pdf", width=10, height=5, units = "in", dpi = 600, bg = "white")


# boxplot for richness vs biclass groups -a more quantitative plot
coropleth_palette <- read.csv("figures/coropleth_palette.txt", header = F)$V1 # actual order

bplot <- ggplot(shp2)+
  geom_boxplot(aes(x=bi_class_ses, y=richness, fill=bi_class_ses), varwidth = T,  show.legend=F)+
  scale_fill_manual(values=coropleth_palette)+
  xlab("SES.PD-SES.PE group")
(bplot2 <- ggdraw() +
    draw_plot(bplot, 0, 0, 1, 1) +
    draw_plot(legend, 0.65, 0.65, 0.3, 0.3))
ggsave("figures/boxplot_PD_PE_SR.png", width=5, height=4, units = "in", dpi = 600, bg = "white")

l1 <- lm(data=shp2, log(richness)~SES.PD*SES.PE)
summary(l1)
b <- mgcv::gam(data=shp2, log(richness)~s(SES.PD)+s(SES.PE))
summary(b)
plot(b,pages=1,residuals=TRUE)  ## show partial residuals
mgcv::gam.check(b)
mgcv::vis.gam(b,theta=30,phi=30,ticktype="detailed")

# sort factor to match facet wrap
shp2$bi_class_ses_l <- factor(shp2$bi_class_ses, levels=c("1-3", "2-3", "3-3", "1-2", "2-2", "3-2", "1-1", "2-1" ,"3-1"))
ggplot(shp2)+
  geom_boxplot(aes(x=bi_class_ses_l, y=richness, col=bi_class_ses_l, fill=bi_class_ses_l), varwidth = F)+
  scale_fill_manual(values=coropleth_palette[c(3,6,9,2,5,8,1,4,7)])+
  scale_color_manual(values=coropleth_palette[c(3,6,9,2,5,8,1,4,7)])+
  facet_wrap("bi_class_ses_l", ncol=3, scales="free_x")+
  scale_x_discrete("", expand = c(0.2,0.2))+
  theme_minimal_hgrid()+
  scale_y_continuous(trans="sqrt", expand=c(0,0))+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.spacing.y = unit(0, "lines"),
    axis.text = element_text(size=8),
    legend.position = "none"
    )


# Intersection with hotspot area a la Myer ----------------------------

# calculate quantitative overlay with hotspots per country as proportion?

all(st_is_valid(m))
m2 <- st_make_valid(m)
#write_sf(m2, "hotspots_2016_1/myer_trans.shp")

s <- as(st_geometry(shp2), "Spatial")
m3 <- as(st_geometry(m), "Spatial")

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

shp2$hotspot_coverage <- a
saveRDS(shp2, "including_hotspot_coverage.rds")

plot_grid(chloropl3, 
ggplot(shp2) + 
  geom_sf(aes(fill=hotspot_coverage, col=hotspot_coverage)), nrow = 2 
)

cor.test(shp2$richness, shp2$hotspot_coverage, method="s") # nothing, probs binary if at all
cor.test(shp2$SES.PE, shp2$hotspot_coverage, method="s")
cor.test(shp2$SES.PD, shp2$hotspot_coverage, method="s")
# SES.PD supports PE hotspot definition?

plot_grid(
ggplot(shp2)+
  geom_boxplot(aes(x=hotspot_coverage>0, y=richness)),
ggplot(shp2)+
  geom_boxplot(aes(x=hotspot_coverage>0, y=SES.PD)),
ggplot(shp2)+
  geom_boxplot(aes(x=hotspot_coverage>0, y=SES.PE)), 
nrow=1)
kruskal.test(shp2$richness, shp2$hotspot_coverage>0) #yes
kruskal.test(shp2$SES.PD, shp2$hotspot_coverage>0) # nope
kruskal.test(shp2$SES.PE, shp2$hotspot_coverage>0) # yes






# Hotspot definition ------------------------------------------------------
# Myer2000 uses endemism (at least x number of endemics) and habitat loss(at least 70% habitat loss)

shp2 <- readRDS("including_hotspot_coverage.rds")

# get a reasonable percent:
kmean_calc <- function(df, ...){
  kmeans(df, scaled = ..., nstart = 30)
}
km2 <- kmean_calc(shp2$SES.PE, 2)
km3 <- kmean_calc(shp2$SES.PE, 3)
km4 <- kmeans(shp2$SES.PE, 4)
km5 <- kmeans(shp2$SES.PE, 5)
km6 <- kmeans(shp2$SES.PE, 6)
km7 <- kmeans(shp2$SES.PE, 7)
km8 <- kmeans(shp2$SES.PE, 8)
km9 <- kmeans(shp2$SES.PE, 9)
km_perc <- function(x)x$betweenss/x$totss
plot(unlist(lapply(list(km2,km3,km4,km5,km6,km7,km8,km9),km_perc)))
plot(shp2$SES.PE~shp2$richness, col=km6$cluster) # top 8 look good, matches 2.5% well

# highest 2.5% of everything:
## Endemism hotspots ---------------------
C <- coldspots(shp2$SES.PE) # coldspots
H <- hotspots(shp2$SES.PE) # hotspots
## Merge endemism values to shapefile of grid cells.
DF <- data.frame(LEVEL3_COD=shp2$LEVEL3_COD, cold=C, hot=H)
DF$SES.PE_spot <- 0
DF$SES.PE_spot[DF$hot==1] <- 1
DF$SES.PE_spot[DF$cold==1] <- -1
shp2 <- merge(shp2, DF[,c("LEVEL3_COD", "SES.PE_spot")], by = "LEVEL3_COD", all = TRUE)
shp2$SES.PE_spot

ggplot(shp2, aes(x=factor(SES.PE_spot), y=hotspot_coverage, label=LEVEL3_COD))+
  #geom_point() +
  geom_text(position="jitter")


## PD hotspots --------------------------
C <- coldspots(shp2$SES.PD) # coldspots
H <- hotspots(shp2$SES.PD) # hotspots
DF <- data.frame(LEVEL3_COD=shp2$LEVEL3_COD, cold=C, hot=H)
DF$SES.PD_spot <- 0
DF$SES.PD_spot[DF$hot==1] <- 1
DF$SES.PD_spot[DF$cold==1] <- -1
shp2 <- merge(shp2, DF[,c("LEVEL3_COD", "SES.PD_spot")], by = "LEVEL3_COD", all = TRUE)



## Richness hotspots -----------------------
C <- coldspots(shp2$richness) # coldspots
H <- hotspots(shp2$richness) # hotspots
DF <- data.frame(LEVEL3_COD=shp2$LEVEL3_COD, cold=C, hot=H)
DF$SR_spot <- 0
DF$SR_spot[DF$hot==1] <- 1
DF$SR_spot[DF$cold==1] <- -1
shp2 <- merge(shp2, DF[,c("LEVEL3_COD", "SR_spot")], by = "LEVEL3_COD", all = TRUE)



## Deforestation hotspots -----------------------
H <- hotspots(shp2$deforest_mean) # hotspots
shp2$deforest_spot <- H

# Myers & PD hotspot matches:
shp2$LEVEL_NAME[shp2$hotspot_coverage>0 & shp2$SES.PD_spot==1]
# Myers & PD hotspot no matches:
shp2$LEVEL_NAME[shp2$hotspot_coverage==0 & shp2$SES.PD_spot==1]

# Myers & PE hotspot matches:
shp2$LEVEL_NAME[shp2$hotspot_coverage>0 & shp2$SES.PE_spot==1]
# Myers & PE hotspot no matches:
shp2$LEVEL_NAME[shp2$hotspot_coverage==0 & shp2$SES.PE_spot==1]

# Myers & SR hotspot matches:
shp2$LEVEL_NAME[shp2$hotspot_coverage>0 & shp2$SR_spot==1]
# Myers & SR hotspot no matches:
shp2$LEVEL_NAME[shp2$hotspot_coverage==0 & shp2$SR_spot==1]

# Myers & Deforestation hotspot matches:
shp2$LEVEL_NAME[shp2$hotspot_coverage>0 & shp2$deforest_spot==1]
# Myers & Deforestation hotspot no matches:
shp2$LEVEL_NAME[shp2$hotspot_coverage==0 & shp2$deforest_spot==1]



## Maps hotspots ---------------------------
min.area <- 1.5e+9
thicc_lines <- shp2[which(shp2$area<min.area),]

# extract overlapping parts for another color
hf <- st_read("hotspots_fixed.gpkg")
wrld_wrap <- st_wrap_dateline(hf, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
m <- st_transform(wrld_wrap, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
# fix sf trouble...
#sf::sf_use_s2(FALSE)
#PE.int <- st_intersection(shp2[shp2$SES.PE_spot==1,], m)
# subtract countries from hotspot areas to not cover them
#hs.int <- st_(m, shp2[shp2$SES.PE_spot==1,])


# put myers map on top
thicc_lines$SES.PE_spot[thicc_lines$SES.PE_spot!=0] # no hot/cold zones in small areas
(PE_hotspot_map <- ggplot(shp2)+
  geom_sf(aes(fill=factor(SES.PE_spot==1)),lwd=.1, col=NA, show.legend = F) + 
  geom_sf(data=thicc_lines, aes(col=factor(SES.PE_spot==1)), show.legend=F, lwd=2)+
  scale_fill_manual(values = c("grey90", "#6575B1"))+
  scale_color_manual(values = c("grey90", "#6575B1"))+
  theme_void()+
    # add Myer layer
  geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
  geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)+
    # add hotspot intersection area
  #geom_sf(data=PE.int, fill="#B16575", color=NA, show.legend=F)+  
    # add text
  geom_hline(yintercept=c(-2343636,0,2343636), size=0.1, lty=c("dashed", "solid","dashed"))+
  annotate("text", x= -17067530, y= 2043636, label= "Tropic of Cancer", hjust=0)+
  annotate("text", x= -17067530, y= -330000, label= "Equator", hjust=0)+
  annotate("text", x= -17067530, y= -2643636, label= "Tropic of Capricorn", hjust=0)+
  ggtitle("Top SES.PE 2.5% with Myers biodiv hotspots"))


PD_hotspot_map <- ggplot(shp2)+
  geom_sf(aes(fill=factor(SES.PD_spot==1)),lwd=0, col=NA, show.legend = F) + 
  geom_sf(data=thicc_lines, aes(col=factor(SES.PD_spot==1)), show.legend=F, lwd=2)+
  scale_fill_manual(values = c("grey90", "red"))+
  scale_color_manual(values = c("grey90", "red"))+
  theme(panel.border = element_blank())+
  ggtitle("SES.PD Hotspots")+
  theme_void()+
  # add Myer layer
  geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
  geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)

table(thicc_lines$SR_spot) # only cold zones in small areas, adjust color
SR_hotspot_map <- ggplot(shp2)+
  geom_sf(aes(fill=factor(SR_spot==1)),lwd=0, col=NA, show.legend = F) + 
  geom_sf(data=thicc_lines, aes(col=factor(SR_spot==1)), show.legend=F, lwd=2)+
  scale_fill_manual(values = c("grey90", "red"))+
  scale_color_manual(values = c("grey90"))+
  theme(panel.border = element_blank())+
  ggtitle("SR Hotspots")+
  theme_void()+
  # add Myer layer
  geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
  geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)

plot_grid(PE_hotspot_map, PD_hotspot_map, SR_hotspot_map, ncol=2)
ggsave("figures/single_hotspots.png", width=10, height=6, units = "in", dpi = 600, bg = "white")


table(thicc_lines$deforest_spot)
(deforest_hotspot_map <- ggplot(shp2)+
  geom_sf(aes(fill=factor(deforest_spot)),lwd=0, col=NA) + 
  geom_sf(data=thicc_lines, aes(col=factor(deforest_spot)), show.legend=F, lwd=2)+
  scale_fill_manual(values = c("grey90", "#6575B1"))+
  scale_color_manual(values = c("grey90", "#6575B1"))+
  theme(panel.border = element_blank())+
  ggtitle("Top 2.5% deforestation countries")+
  theme_void()+
  # add Myer layer
  geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
  geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)
)

## Visualize overlap of hotspots -------------------------------
data <- bi_class(shp2, x = deforest_mean, y = SES.PE, style = "jenks", dim = 3)
data$bi_class_ses <- data$bi_class
shp2$bi_class_PE_deforest <- data$bi_class_ses
rm(data)

thicc_lines <- data[which(shp2$area<min.area),]
(PE_deforest <- ggplot() +
    geom_sf(shp2, mapping = aes(fill = bi_class_PE_deforest), color = NA, size = 0.1, show.legend = FALSE) +
    bi_scale_fill(pal = "BlueGold", dim = 3, na.value="white") +
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=bi_class), show.legend=F)+
    bi_scale_color(pal = "BlueGold", dim = 3, na.value="white") +
    bi_theme()+
    #geom_sf(data=m[m$Type=="hotspot area",], fill="green", color=NA, alpha=0.5, show.legend=F)+
    #geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="green", color=NA, alpha=0.5, show.legend=F)+
    ggtitle("SES.PE and deforestation")+
    theme(title = element_text(size=10)))
legend <- bi_legend(pal = "BlueGold",
                    dim = 3,
                    xlab = "deforestation",
                    ylab = "SES.PE ",
                    size = 8)

(ses_pe_deforest_map <- ggdraw() +
    draw_plot(PE_deforest, 0, 0, 1, 1) +
    draw_plot(legend, 0.05, 0.25, 0.2, 0.2))
ggsave("figures/choropleth_PE_deforest.png", width=10, height=5, units = "in", dpi = 600, bg = "white")

data <- bi_class(shp2, x = deforest_mean, y = SES.PD, style = "jenks", dim = 3)
data$bi_class_ses <- data$bi_class
shp2$bi_class_PD_deforest <- data$bi_class_ses
rm(data)
thicc_lines <- shp2[which(shp2$area<min.area),]
(PD_deforest <- ggplot() +
    geom_sf(shp2, mapping = aes(fill = bi_class_PD_deforest), color = NA, size = 0.1, show.legend = FALSE) +
    bi_scale_fill(pal = "BlueGold", dim = 3, na.value="white") +
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=bi_class), show.legend=F)+
    bi_scale_color(pal = "BlueGold", dim = 3, na.value="white") +
    bi_theme()+
    #geom_sf(data=m[m$Type=="hotspot area",], fill="green", color=NA, alpha=0.5, show.legend=F)+
    #geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="green", color=NA, alpha=0.5, show.legend=F)+
    ggtitle("SES.PD and deforestation")+
    theme(title = element_text(size=10)))
legend <- bi_legend(pal = "BlueGold",
                    dim = 3,
                    xlab = "deforestation",
                    ylab = "SES.PD ",
                    size = 8)

(ses_pd_deforest_map <- ggdraw() +
    draw_plot(PD_deforest, 0, 0, 1, 1) +
    draw_plot(legend, 0.05, 0.25, 0.2, 0.2))
ggsave("figures/choropleth_PD_deforest.png", width=10, height=5, units = "in", dpi = 600, bg = "white")


# all hotspots in one map
thicc_lines <- shp2[which(shp2$area<min.area),]
al <- 0.3
ggplot()+
  # world layer
    geom_sf(data=shp2, col=NA, fill="grey95", lwd=0, show.legend = F) +
  # PE layer
    geom_sf(data=shp2[shp2$SES.PE_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
    geom_sf(data=thicc_lines[thicc_lines$SES.PE_spot==1, ], col="#6575B1", show.legend=F, lwd=2) +
  # PD layer
    geom_sf(data=shp2[shp2$SES.PD_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
    geom_sf(data=thicc_lines[thicc_lines$SES.PD_spot==1, ], col=alpha("#6575B1",0.3), show.legend=F, lwd=2) +
  # SR layer
    geom_sf(data=shp2[shp2$SR_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
    geom_sf(data=thicc_lines[thicc_lines$SR_spot==1, ], col="#6575B1", show.legend=F, lwd=2) +
  # style settings
    ggtitle("All hotspots, layered")#+
#  theme_void()+
#  # add Myer layer
#  geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
#  geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)

# get 
x <- st_drop_geometry(shp2[,c("SR_spot","SES.PD_spot", "SES.PE_spot")])
x$SR_spot <- ifelse(x$SR_spot<0, NA,x$SR_spot)
x$SES.PE_spot <- ifelse(x$SES.PE_spot<0, NA,x$SES.PE_spot)
x$SES.PD_spot <- ifelse(x$SES.PD_spot<0, NA,x$SES.PD_spot)
shp2$hotspot_sum <- rowSums(x, na.rm = T)

ggplot()+
  # world layer
  geom_sf(data=shp2, aes(fill=factor(hotspot_sum)),col=NA, lwd=0) +
  geom_sf(data=thicc_lines[which(thicc_lines$LEVEL3_COD %in% shp2$LEVEL3_COD[shp2$hotspot_sum!=0]), ],
          col="#6575B1", show.legend=F, lwd=2) +
  # hotspot sum laye...
  # geom_sf(aes(fill=hotspot_sum), col=NA, lwd=0, show.legend = F, alpha=al) + 
  # geom_sf(data=thicc_lines[thicc_lines$SES.PE_spot==1, ], col="#6575B1", show.legend=F, lwd=2) +
  # # PD layer
  # geom_sf(data=shp2[shp2$SES.PD_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
  # geom_sf(data=thicc_lines[thicc_lines$SES.PD_spot==1, ], col=alpha("#6575B1",0.3), show.legend=F, lwd=2) +
  # # SR layer
  # geom_sf(data=shp2[shp2$SR_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
  # geom_sf(data=thicc_lines[thicc_lines$SR_spot==1, ], col="#6575B1", show.legend=F, lwd=2) +
  # # style settings
   scale_fill_manual(values = c("grey95", "#6575B1"))+
  # scale_color_manual(values = c("grey90", "#6575B1"))+
  ggtitle("Sum all hotspots")+
   # add Myer layer
   geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
   geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)




# # Ternary plots -------------------------------------------------------------
# library(ggtern)
  # shp2$richness.norm <- normalized(shp2$richness)
# 
# # ggtern(shp2, aes(x=PD_obs.norm, y=richness.norm, z=PE_obs.norm))+
# #   geom_confidence_tern()+
# #   tern_limit(T = 1.05, L = 1.05, R = 1.05)+
# #   geom_point()+
# #   theme_bvbw()+
# #   Tlab("SR", labelarrow = "species richness")+
# #   Llab("PD", labelarrow = "phylogenetic diversity")+
# #   Rlab("PE", labelarrow = "phylogenetic endemism")



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
# res[1:9, c("richness", "PD_obs", "zscore", "PE")]
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






# High threat regions -----------------------------------------------------


nb <- spdep::poly2nb(shp, row.names = shp$LEVEL3_COD)
col.W <- nb2listw(nb, style="W", zero.policy = TRUE)
thr <- shp[,grep("LEVEL|SES|_spot|zspot|PE|hfp|deforest|_change|FC|PC", names(shp))]
thr <- st_drop_geometry(thr)

tes <- apply(thr[,grep("hfp|deforest|_change|FC|PC", names(thr))], 2, lee.test, y=thr$SES.PD, 
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
  scale_fill_manual(values = c("blue", "grey90", "red"))+
  scale_color_manual(values = c("blue", "grey90", "red"))+
  theme(panel.border = element_blank())+
  ggtitle("Observed HFP Hotspots and Coldspots")
ggsave("figures/PD_hot_cold.png", width=7, height=4, units = "in", dpi = 600, bg = "white")


