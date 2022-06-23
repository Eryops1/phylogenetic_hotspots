
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
shp <- shp[,-grep("obs_p|obs_rank|reps|LEVEL2|LEVEL1|LEVEL_3_CO|ID|\\.3|_rw|CONTI|REGION", names(shp))]
names(shp)<- gsub("\\.1", "_mean", names(shp))
names(shp)<- gsub("\\.2", "_sd", names(shp))

# remove BOU that has no data
shp <- shp[!shp$LEVEL3_COD=="BOU",]


# Explorative plots -------------------------------------------------------

plot(st_drop_geometry(shp)[,c(1,4:12,13:23)])
plot(st_drop_geometry(shp)[,c(1,4:12,24:35)])
plot(st_drop_geometry(shp)[,c(1,4:12,36:57)])




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


plot_grid(ggplot(dat.sub, aes(y=area, group=is.na(hfp_mean)))+
  geom_boxplot(varwidth = T)+
  scale_y_log10()+
  xlab("HFP == NA (positive=T: n=27)"),
ggplot(dat.sub, aes(y=PD_obs, group=is.na(hfp_mean)))+
  geom_boxplot(varwidth = T)+
  xlab("HFP == NA (positive=T: n=27)"),
ggplot(dat.sub, aes(y=SES.PD, group=is.na(hfp_mean)))+
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
          ggplot(data=s@data, aes(x=pd_rand_mean, y=PD_obs, col=log(richness)))+
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
                    tuneGrid=gbmGrid){
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
  PD_list <- mclapply(1:length(s), seeds = s, run.gbm,
                      data = dat_no.na, mc.cores=n_cores)
)

#pd_list <- unlist(PD_list, recursive = F)
eq <- matrix(sapply(PD_list, "[[", 1), ncol=72, byrow = T)
eq <- reshape2::melt(eq)
ggplot(eq, aes(value))+
  geom_bar()+
  facet_wrap(~Var2, scales = "free")+
  theme(axis.text.x = element_text(size=6, angle = 45, hjust=1),
        strip.background = element_blank())

eqr <-tapply(eq$value, eq$Var2, function(x){names(sort(table(x), decreasing = T))})
fin_sr <- sort(tapply(eq$Var2, eq$value, function(x){sum(x)/length(x)}), decreasing=F)

# select variables that score lowest position relative to their number of
# appearances in different positions (the lower the more important)
fin_sr <- data.frame(fin_sr, variable=names(fin_sr))


fin_sr <- fin_sr[order(fin_sr$fin_sr, decreasing = T),]
fin_sr$variable <- factor(fin_sr$variable, levels = fin_sr$variable)

ggplot(fin_sr, aes(y=variable, x=fin_sr))+
  xlab("Average variable position")+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=9))
ggsave("figures/varImp_SES_PD_gbm10_runs.png", width=7, height=14, units = "in", dpi = 600)






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

# format to sf object for plotting
shp <- readRDS("for_plotting.rds")
shp <- st_as_sf(shp)

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

# +100 tweak.....
#lcol <- min(thicc_lines$SES.PD, na.rm = T)/max(shp$SES.PD, na.rm = T)
#ucol <- max(thicc_lines$SES.PD, na.rm = T)/max(shp$SES.PD, na.rm = T)
(pd_ses_map <- ggplot(shp) + 
    geom_sf(aes(fill=SES.PD),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=SES.PD), show.legend=F)+
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
# (pd_rand_map <- ggplot(shp) + 
#     geom_sf(aes(fill=pd_rand_mean),lwd=0, col=NA) + 
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



#     theme(legend.position = c(0.18, 0.3),
#           legend.key.height = unit(6,"mm"),
#           legend.background = element_blank(),
#           legend.key = element_blank(),
#           panel.background = element_blank(),
#           panel.border = element_blank(),
#           text = element_text(size = 10)
#     )+
#     xlab(" ")+
#     ggtitle("Phylogenetic endemism")
# )

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

plot_grid(pd_map, pd_ses_map, pd_richness_map, sr_map, nrow = 4)
 #ggsave("figures/maps.png", width=8, height=16, units = "in", dpi = 600, bg = "white")

par(mfrow=c(1,2))
plot(shp2$richness, shp2$SES.PD)
plot(shp2$richness, shp2$PD_obs/shp2$richness)

cor.test(shp2$richness, shp2$PD_obs/shp2$richness, method="s")


# lm1 <- lm(data=shp, PD_obs~mrd + soil + mat_mean + tra_mean + pre_mean + prs_mean + elev_range+
#      tri + area + hfp.1 + deforest.1 + bio1.1 + bio12.1)




# Bivariate colorpleth map -------------------------------------------------
# map that shows areas with two variables
library(biscale)
meyer <- readOGR("hotspots_2016_1/hotspots_2016_1.shp")
#sf_use_s2(FALSE)
# library(rgdal)
# meyer <- spTransform(meyer, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"))

# library(terra)
# meyer2 <- vect("hotspots_2016_1/hotspots_2016_1.shp")
# is.valid(meyer2)
# #makeValid(meyer2)
# 
# behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
# p <- terra::project(meyer2, behrmann)
# 
# m <- st_as_sf(meyer)
# st_is_valid(m)
# m = st_buffer(m, 0)  ## fixes some issues
# m <- st_make_valid(m)
# 
# m <- st_transform(m, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") # Behrmann

# extract problem polygons
m <- st_as_sf(meyer)
#tmp <- st_union(m, by_feature = TRUE)
#plot(tmp)
wrld_wrap <- st_wrap_dateline(m, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
m <- st_transform(wrld_wrap, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
plot(m)

# inv <- which(!st_is_valid(m))
# m$NAME[inv]
# m$Type[inv]
# 
# # subset
# msub <- m[-inv,]
# msub <- st_transform(msub, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
# 
# mprob <- m[inv,]
# mprob <- st_shift_longitude(mprob)
# p1 <- st_cast(mprob[1,], "POLYGON")
# p1 <- st_union(p1)
# p1 <- st_cast(p1, "POLYGON")
# 
# p1 <- st_union(p1, is_coverage = T) # makes a multipolygon again....
# 
# mprob <- st_union(mprob)
# 
# mprob <- st_transform(mprob, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs",
#                         )
# 



data <- bi_class(shp, x = PD_obs, y = PE, style = "jenks", dim = 3)
thicc_lines <- data[which(data$area<min.area),]
(map <- ggplot() +
  geom_sf(data = data, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal = "GrPink2", dim = 3, na.value="white") +
  geom_sf(data=thicc_lines, lwd=1.5, aes(col=bi_class), show.legend=F)+
  bi_scale_color(pal = "GrPink2", dim = 3, na.value="white") +
  bi_theme()+
  geom_sf(data=m, aes(fill="green", col="grey", alpha=0.5))+
  ggtitle("PD and WE"))
legend <- bi_legend(pal = "GrPink2",
                    dim = 3,
                    xlab = "PD",
                    ylab = "PE ",
                    size = 8)
ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.05, 0.25, 0.2, 0.2)


#ggsave("figures/phylo_pure_colorpleth_map.png", width=7, height=4, units = "in", dpi = 600, bg = "white")


data <- bi_class(shp, x = SES.PD, y = PE, style = "jenks", dim = 3)
thicc_lines <- data[which(data$area<min.area),]
(map <- ggplot() +
    geom_sf(data = data, mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) +
    bi_scale_fill(pal = "GrPink2", dim = 3, na.value="white") +
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=bi_class), show.legend=F)+
    bi_scale_color(pal = "GrPink2", dim = 3, na.value="white") +
    bi_theme()+
    ggtitle("PD and WE"))
legend <- bi_legend(pal = "GrPink2",
                    dim = 3,
                    xlab = "SES.PD",
                    ylab = "PE ",
                    size = 8)
ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.05, 0.25, 0.2, 0.2)
#ggsave("figures/taxo_colorpleth_map.png", width=7, height=4, units = "in", dpi = 600, bg = "white")

library(ggtern)
source("../Global_Drivers_2022/scripts/0_functions.R")
shp2.sd <- shp2
shp2.sd$PD_obs <- ztrans(shp2.sd$PD_obs)
shp2.sd$richness <- ztrans(shp2.sd$richness)
shp2.sd$PE <- ztrans(shp2.sd$PE)
shp2.sd$SES.PD <- ztrans(shp2.sd$SES.PD)
shp2.sd$SES.PD <- shp2.sd$SES.PD+5

ggtern(shp2.sd, aes(x=PD_obs, y=richness, z=PE))+
  geom_confidence_tern()+
  geom_point(col="grey")+
  theme_bvbw()

ggtern(shp2.sd, aes(x=SES.PD, y=richness, z=PE))+
  geom_confidence_tern()+
  geom_point(col="grey")+
  theme_bvbw()


# classic hotspots a la meyer ---------------------------------------------
meyer <- readOGR("hotspots_2016_1/hotspots_2016_1.shp")
m <- st_as_sf(meyer)

ggplot(m[m$Type=="outer limit",])+
  geom_sf()
ggplot(m[m$Type=="hotspot area",])+
  geom_sf()

ggplot(m)+
  geom_sf(aes(col=NA))+
  geom_sf(data=shp, fill=richness)





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
C <- coldspots(shp$SES.PD) # coldspots
H <- hotspots(shp$SES.PD) # hotspots
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
  scale_fill_manual(values = c("blue", "grey80", "red"))+
  scale_color_manual(values = c("blue", "grey80", "red"))+
  theme(panel.border = element_blank())+
  ggtitle("Observed HFP Hotspots and Coldspots")
ggsave("figures/PD_hot_cold.png", width=7, height=4, units = "in", dpi = 600, bg = "white")


