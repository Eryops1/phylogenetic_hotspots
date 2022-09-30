
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
gc()
library(sf)
library(ggplot2)
theme_set(theme_bw()+theme(text=element_text(size=8), panel.grid=element_blank()))
library(cowplot)
library(beepr)
library(biscale)
library(spdep)
library(scico)
library(ggpattern)
library(data.table)
library(terra)
if(!dir.exists("figures"))dir.create("figures")
source("99_functions.R")


gallpeters_projection <- "+proj=cea +lon_0=0 +x_0=0 +y_0=0 +lat_ts=45 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
# Gall-Peters. Horizontally compressed version of the Lambert equal-area.
# Standard parallels at 45°N/S. Aspect ratio of ~1.6. Similar is Balthasar
# projection with standard parallels at 50°N/S.
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

my_projection <- gallpeters_projection
shp <- readRDS("fin_shp.rds")

# remove not needed data
shp <- shp[,-grep("obs_p|obs_rank|reps|LEVEL2|LEVEL1|LEVEL_3_CO|ID|\\.3|_rw|CONTI|REGION", names(shp))]
names(shp)<- gsub("\\.1", "_mean", names(shp))
names(shp)<- gsub("\\.2", "_sd", names(shp))
#st_write(shp, "fin_shape_for_gis_checks.gpkg")




# Standardization effects ------------------------------------------------
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
bonbon_PD <-
  ggplot(df2[grep("PD", df2$name),], aes(name, value, color = richness, group=LEVEL_NAME))+
  geom_bump(size = 1, position = "identity", smooth=sm) +
  scale_color_scico("Species richness", trans="sqrt", alpha = alevel, palette="batlow", guide = "none")+
  scale_x_discrete("", expand = c(0.05,0.05))+
  scale_y_discrete("rank")+
  theme_minimal()

bonbon_PE <- ggplot(df2[grep("PE", df2$name),], aes(name, value, color = richness, group=LEVEL_NAME))+
  geom_bump(size = 1, position = "identity", smooth=sm) +
  scale_color_scico("Species richness", trans="sqrt", alpha = alevel, palette="batlow", guide = "none")+
  scale_x_discrete("", expand = c(0.05,0.05))+
  scale_y_discrete("rank")+
  geom_point(size = 0, aes(fill=richness)) +
  scale_fill_scico("Species richness", trans="sqrt", alpha = 1, palette="batlow")+
  theme_minimal()

plot_grid(bonbon_PD, bonbon_PE, ncol=2, rel_widths = c(.41,.59))  
ggsave("figures/standardization_effects.png", width=8, height=4, units = "in", dpi = 600, bg = "white")

# quantitative:
pd_rank_changes <- sum(abs(tapply(df2$value[grep("PD", df2$name)], df2$LEVEL_NAME[grep("PD", df2$name)], diff)))
pe_rank_changes <- sum(abs(tapply(df2$value[grep("PE", df2$name)], df2$LEVEL_NAME[grep("PE", df2$name)], diff)))
pd_rank_changes/length(unique(df2$LEVEL_NAME)); pe_rank_changes/length(unique(df2$LEVEL_NAME))
## --> rank changes more pronounced in PD than in PE   

df$PD_changes <- df$SES.PD_rank-df$PD_rank
df$PE_changes <- df$SES.PE_rank-df$PE_rank
df$hs <- NA
df$hs[df$LEVEL_NAME %in% c("Borneo", "China South-Central", "China Southeast", "Cape Provinces", 
                           "Queensland", "Thailand", "Western Australia", "East Himalaya", "Guatemala", "India", "Malaya", 
                           "Mexico Gulf", "Myanmar", "Philippines", "Sumatera", "Vietnam")] <-"hs" 
ggplot(df, aes(x=richness, label=LEVEL_NAME))+
  geom_point(aes(y=PD_changes), col="#52548D", alpha=.3)+
  #  geom_label(data=df[df$hs=="hs",], aes(x=richness, y=PD_changes), size=3, label.size=0)+
  geom_point(aes(y=PE_changes), col="#EFB984", alpha=.3)+
  #  geom_label(data=df[df$hs=="hs",], aes(x=richness, y=PE_changes), size=3, label.size=0)
  geom_smooth(aes(y=PD_changes), col="#52548D")+
  geom_smooth(aes(y=PE_changes), col="#EFB984")+
  scale_x_continuous(trans="sqrt")+coord_cartesian(expand=F, xlim=c(30,23000))+
  scale_y_continuous("Rank changes")+
  geom_hline(yintercept=0, size=0.1)
ggsave("figures/rank_changes.png", width=3, height=3, units = "in", dpi = 300, bg = "white")

# SES.PD vs PD_obs colored for SR
ggplot(shp2, aes(x=PD_obs, y=SES.PD, size=3))+
  geom_point(aes(color=richness))+
  scale_size(guide="none")+
  scale_color_scico(trans="sqrt")
ggsave("figures/ses_PD_vs_PD_obs.png", width=5, height=4, units = "in", dpi = 600, bg = "white")


### compare to brody sandels rarefaction approach ####
rar  <- readRDS("brody_rarefaction.rds")
rar.df <- data.frame(pd_obs_rar = tapply(rar$PD_obs, rar$grids, median),
                     zscore_rar = tapply(rar$zscore, rar$grids, median),
                     zscore_sd_rar = tapply(rar$zscore, rar$grids, sd),
                     zscore_rand_rar = tapply(rar$pd_rand_mean, rar$grids, median))
rar.df$LEVEL3_COD <- row.names(rar.df)
shp2 <- merge(shp2, rar.df, all.x=T)

plot(shp2$PD_obs, shp2$pd_obs_rar)
plot(shp2$SES.PD, shp2$zscore_rar)
cor.test(shp2$SES.PD, shp2$zscore_rar, method="s")

# sesPD_rarefaction
thicc_lines <- shp2[which(shp2$area<min.area),]
(pd_ses_rarefaction_map <- ggplot(shp2) + 
    geom_sf(aes(fill=zscore_rar),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=zscore_rar), show.legend=F)+
    scale_colour_scico("zscore_rar", palette="batlow")+
    scale_fill_scico("zscore_rar", palette="batlow")+  
    theme_void()+
    theme(legend.position = c(0.22, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")+
    ggtitle("I dont think its reasonable to subsample 23k to 41, no matter how many times?")
)


#botta-dukat: I recommend checking symmetry of
#null-distribution before calculating the SES value. If the distribution is
#skewed, I recommend either log-transformation of the test statistic, or using
#probit-transformed p-value as effect size measure.

hist(shp2$pd_rand_mean)
library(psych)
describe(shp2$pd_rand_mean)
# null distribution is skewed. transform:
describe(log(shp2$pd_rand_mean))
describe(sqrt(shp2$pd_rand_mean)) # works

# load results from sqrt transformed stuff
sesPD_norm <- readRDS("sesPD_norm.rds")
norm.df <- data.frame(pd_obs_norm = sesPD_norm$PD_obs,
                      zscore_norm = sesPD_norm$zscore, 
                      LEVEL3_COD=sesPD_norm$grids, 
                      zscore_rand_norm=sesPD_norm$pd_rand_mean)
shp2 <- merge(shp2, norm.df, all.x=T)

plot(shp2$PD_obs, shp2$pd_obs_norm)
plot(shp2$SES.PD, shp2$zscore_norm)
abline(a=0,b=1)
cor.test(shp2$SES.PD, shp2$zscore_norm, method="s")

plot(shp2$richness, shp2$SES.PD, ylim=c(-90,5))
points(shp2$richness, shp2$zscore_norm, col="red")
plot(shp2$richness, shp2$zscore_rar, col="blue")

thicc_lines <- shp2[which(shp2$area<min.area),]
(pd_ses_norm_map <- ggplot(shp2) + 
    geom_sf(aes(fill=zscore_norm),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=zscore_norm), show.legend=F)+
    scale_colour_scico("zscore_norm", palette="batlow")+
    scale_fill_scico("zscore_norm", palette="batlow")+  
    theme_void()+
    theme(legend.position = c(0.22, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")
)

hist(shp2$pd_rand_mean)
hist(shp2$pd_obs_rar)

# does the ratio of sd to PD_rand_mean change with SR??
names(shp2)
ggplot(shp2)+
  geom_point(aes(x=richness, y=pd_rand_sd))
ggplot(shp2)+
  geom_point(aes(x=richness, y=pd_rand_mean))

# ok, but what if the distance between PD_obs and PD_rand_mean also scales with richness?
ggplot(shp2)+
  geom_point(aes(x=richness, y=PD_obs - pd_rand_mean))
ggplot(shp2)+
  geom_point(aes(x=richness, y=SES.PD))

ggplot(shp2)+
  geom_point(aes(x=richness, y=pd_rand_mean))+
  scale_x_continuous(trans="sqrt")+
  geom_point(aes(x=richness, y=PD_obs))


# shp2$SES.PD.norm <- normalized(shp2$SES.PD)
# shp2$PD_obs.norm <- normalized(shp2$PD_obs)
# 
# thicc_lines <- shp2[which(shp2$area<min.area),]
# lcol <- min(thicc_lines$SES.PD.norm)/max(shp2$SES.PE.norm)
# ucol <- max(thicc_lines$SES.PD.norm)/max(shp2$SES.PE.norm)
# (pd_standardization_effects_map <- ggplot(shp2) + 
#     geom_sf(aes(fill=PD_obs.norm-SES.PD.norm)) + 
#     geom_sf(data=thicc_lines, lwd=1.5, aes(col=PD_obs.norm-SES.PD.norm), show.legend=F)+
#     scale_color_gradient2("stand. differences",low = "red", mid = "white", high = "blue")+
#     scale_fill_gradient2("stand. differences",low = "red", mid = "white", high = "blue")+
#     theme_void()+
#     ggtitle("PD[0,1] - SES.PD[0,1]. Blue=SES is smaller, red=SES is bigger than raw values")
# )
