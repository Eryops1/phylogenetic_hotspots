
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
cmat <- cor_mat(dat_no.na[,-grep("LEVEL|hfp|deforest|bio|change|PC|FC|pd_rand_mean",names(dat_no.na))], method = "s")
cpmat <- cor_pmat(dat_no.na[,-grep("LEVEL|hfp|deforest|bio|change|PC|FC|pd_rand_mean",names(dat_no.na))], method = "s")
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

library(spatialRF)

# subset spatial object to dat_no.na variables
shp <- shp[,which(names(shp) %in% names(dat_no.na))]
shp2 <- na.omit(shp)

#coordinates of the cases
shp2$centroids <- st_centroid(shp2) %>% 
  st_coordinates()
shp2$y <- shp2$centroids[,1] # x=lng
shp2$x <- shp2$centroids[,2]
xy <- st_drop_geometry(shp2[, c("x", "y")])

#distance matrix
dm <- dist(xy, method = "euclidean", diag = TRUE, upper = TRUE)
dm <- as.matrix(dm)

#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 100000, 500000, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000)

#random seed for reproducibility
random.seed <- 1

# spatialRF::plot_training_df(
#   data = shp2,
#   dependent.variable.name = SES.PD,
#   predictor.variable.names = predictor.variable.names,
#   ncol = 5,
#   point.color = viridis::viridis(100, option = "F"),
#   line.color = "gray30"
# )
# 
# ## assess autocorrelation
# 
# spatialRF::plot_training_df_moran(
#   data = shp2,
#   dependent.variable.name = dependent.variable.name,
#   predictor.variable.names = predictor.variable.names,
#   distance.matrix = dm,
#   distance.thresholds = distance.thresholds,
#   fill.color = viridis::viridis(
#     100,
#     option = "F",
#     direction = -1
#   ),
#   point.color = "gray40"
# )

library(spdep)
# create neighbor list
nb <- spdep::poly2nb(shp2, row.names = shp2$LEVEL3_COD)
names(nb) <- shp2$LEVEL3_COD

# distance object
col.W <- nb2listw(nb, style="W", zero.policy = TRUE)


# *** Moran's I -----------------------------------------------------------

#moran.test(tmp$richness, col.W, zero.policy = TRUE)
tmpm <- st_drop_geometry(shp2)
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

# use this to check where PD hotspots and threats overlap?

nb <- spdep::poly2nb(tmp, row.names = tmp$LEVEL3_COD)
col.W <- nb2listw(nb, style="W", zero.policy = TRUE)


lee.test(tmpm$SES.PD, tmpm$SES.PD, listw=col.W, zero.policy = T)
lee.mc(tmpm$SES.PD, tmpm$SES.PD, nsim=99, col.W, zero.policy=TRUE)

les <- apply(tmpm[,-which(names(tmpm)=="pd_rand_mean")], 2, lee.test, y=tmpm$SES.PD, 
             listw=col.W, zero.policy = TRUE, alternative="two.sided")
les.df <- sapply(les, "[[", "estimate")
les.df <- rbind(les.df, sapply(les, "[[", "p.value"))
row.names(les.df) <- c("Lee", "expect", "var", "pvalue")
les.df <- as.data.frame(t(les.df))
les.df$twoSD <- 2*sqrt(les.df$var)

# sort the cor strength
les.df <- les.df[order(les.df$Lee, decreasing = T),]
co <- row.names(les.df)
row.names(les.df) <- factor(row.names(les.df), levels=co)
les.df$vars <- row.names(les.df)
les.df$vars <- factor(les.df$vars, levels=co)
ggplot(les.df, aes(y=expect, x=vars))+
  geom_linerange(aes(ymin=expect-twoSD, ymax=expect+twoSD), col="grey70")+
  geom_point(aes(y=Lee, x=row.names(les.df), col=factor(pvalue<0.05)), show.legend = F)+
  ylab("Lee's L with SES.PD")+
  scale_color_manual("p<0.05",values=c("grey20", "red"))+
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
