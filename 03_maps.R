
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
gc()
library(sf)
library(ggplot2)
theme_set(theme_bw()+theme(text=element_text(size=8), panel.grid=element_blank()))
library(cowplot)
#library(rgdal)
library(beepr)
library(biscale)
library(spdep)
library(scico)
library(ggpattern)
library(data.table)
if(!dir.exists("figures"))dir.create("figures")
source("99_functions.R")


# Load data ---------------------------------------------------------------

gallpeters_projection <- "+proj=cea +lon_0=0 +x_0=0 +y_0=0 +lat_ts=45 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
# Gall-Peters. Horizontally compressed version of the Lambert equal-area.
# Standard parallels at 45°N/S. Aspect ratio of ~1.6. Similar is Balthasar
# projection with standard parallels at 50°N/S.
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
# chose:
my_projection <- gallpeters_projection
shp <- readRDS("fin_shp.rds")

# remove not needed data
shp <- shp[,-grep("obs_p|obs_rank|reps|LEVEL2|LEVEL1|LEVEL_3_CO|ID|\\.3|_rw|CONTI|REGION", names(shp))]
names(shp)<- gsub("\\.1", "_mean", names(shp))
names(shp)<- gsub("\\.2", "_sd", names(shp))

#st_write(shp, "fin_shape_for_gis_checks.gpkg")



# Maps SR, SES.PD,  WE, SES.PE ------------------------------------------------

# transform projection
shp <- st_transform(shp, my_projection)

min.area <- 6e+9
thicc_lines <- shp[which(shp$area<min.area),]


lcol <- min(thicc_lines$PD_obs)/max(shp$PD_obs)
ucol <- max(thicc_lines$PD_obs)/max(shp$PD_obs)
(pd_map <- ggplot(shp) + 
    geom_sf(data=shp, aes(fill=PD_obs),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=PD_obs), show.legend=F)+
    scale_colour_scico("PD", palette = "batlow", trans = "sqrt", 
                            begin = lcol, end = sqrt(ucol))+
    scale_fill_scico("PD", palette = "batlow", trans="sqrt")+ #, 
    theme(legend.position = c(0.22, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10))+
    xlab(" ")
)
lcol <- min(thicc_lines$PE_obs)/max(shp$PE_obs)
ucol <- max(thicc_lines$PE_obs)/max(shp$PE_obs)
(pe_map <- ggplot(shp) + 
    geom_sf(aes(fill=PE_obs),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=PE_obs), show.legend=F)+
    scale_colour_scico("PE", palette="batlow", trans="sqrt",
                           begin = lcol, end = ucol)+
    scale_fill_scico("PE", palette="batlow",trans="sqrt")+ #, 
    theme(legend.position = c(0.18, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")
)
# Simple SR (excluding Antarctica)
lcol <- min(thicc_lines$richness)/max(shp$richness)
ucol <- max(thicc_lines$richness)/max(shp$richness)
(sr_map <- ggplot(shp[!shp$LEVEL3_COD=="ANT",]) + 
    geom_sf(aes(fill=richness),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=richness), show.legend=F)+
    scale_colour_scico("SR", palette="batlow", trans = "sqrt", 
                           begin = lcol, end = sqrt(ucol))+
    scale_fill_scico("SR", palette="batlow", trans = "sqrt")+ #, 
    theme_void()+
    theme(legend.position = c(0.22, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )
)


shp2 <- shp[!is.na(shp$SES.PD),]
thicc_lines <- shp2[which(shp2$area<min.area),]

# simple WE
lcol <- min(thicc_lines$WE)/max(shp$WE)
ucol <- max(thicc_lines$WE)/max(shp$WE)
(we_map <- ggplot(shp2) + 
    geom_sf(aes(fill=WE),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=WE), show.legend=F)+
    scale_colour_scico("WE", palette="batlow", trans = "sqrt", 
                           begin = lcol, end = sqrt(ucol))+
    scale_fill_scico("WE", palette="batlow", trans = "sqrt")+ #, 
    theme_void()+
    theme(legend.position = c(0.22, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )+
    xlab(" ")
)


lcol <- 1-min(thicc_lines$SES.PD)/(min(shp2$SES.PD)-max(shp2$SES.PD))
ucol <- max(thicc_lines$SES.PD)/max(shp2$SES.PD)
(pd_ses_map <- ggplot(shp2) + 
    geom_sf(aes(fill=SES.PD),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=SES.PD), show.legend=F)+
    scale_colour_scico("SES.PD", palette="batlow", begin=lcol, end=ucol)+
    scale_fill_scico("SES.PD", palette="batlow")+  
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

lcol <- min(thicc_lines$SES.PE+abs(min(shp2$SES.PE)))/diff(range(shp2$SES.PE)) 
ucol <- max(thicc_lines$SES.PE)/max(shp2$SES.PE)
(pe_ses_map <- ggplot(shp2) + 
    geom_sf(aes(fill=SES.PE),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=SES.PE), show.legend=F)+
    scale_colour_scico("SES.PE", palette="batlow",
                           begin = lcol, end = ucol)+
    scale_fill_scico("SES.PE", palette="batlow")+ 
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



plot_grid(sr_map, we_map, pd_ses_map, pe_ses_map, ncol = 2, 
          labels=c("A","B","C","D"), label_fontface=1)
ggsave("figures/maps.png", width=12.5, height=7, units = "in", dpi = 600, bg = "white")















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
















# Myer hotspot proportions ------------------------------------

h <- st_read("hotspots_fixed.gpkg")
h <- st_wrap_dateline(h, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
h <- st_transform(h, my_projection)

s <- as(st_geometry(shp2), "Spatial")
m3 <- as(st_geometry(h), "Spatial")

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
# saveRDS(shp2, "including_hotspot_coverage.rds")
shp2 <- readRDS("including_hotspot_coverage.rds")

# how many myers hotspots covered?
tmp <- shp2[shp2$PD_hotspot=="3-3" |shp2$PE_hotspot %in% c("4-3", "3-4", "4-4"), ]
s <- as(st_geometry(tmp), "Spatial")

res <- raster::intersect(m3, s) 
# first = myer hotspot, second=evol hotspot
myer_number <- gsub(" [0-9]{1,2}", "", names(res))
evol_number <- gsub("[0-9]{1,2} ", "", names(res))
unique(h$NAME[as.numeric(myer_number)])
unique(tmp$LEVEL_NAME[as.numeric(evol_number)])


library(scico)
ggplot(shp2) + 
  geom_sf(aes(col=hotspot_coverage), show.legend=F)+
  geom_sf(aes(fill=hotspot_coverage), col=NA)+
  scale_fill_scico("hotspot \ncoverage", palette="batlow", end=0.7)+
  scale_color_scico(palette="batlow", end=0.7)+
  theme(legend.position = c(0.18, 0.3),
        legend.key.height = unit(6,"mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title=element_text(size=9))+
  coord_sf(expand=F)
ggsave("figures/hotspot_coverage.png", units="in", dpi=300, width=7, height=4.1)

## correlations with PD,PE,SR,WE ----
cor.test(shp2$richness, shp2$hotspot_coverage, method="p") 
cor.test(shp2$SES.PE, shp2$hotspot_coverage, method="p")
cor.test(shp2$SES.PD, shp2$hotspot_coverage, method="p")
cor.test(shp2$area, shp2$hotspot_coverage, method="s")
cor.test(shp2$WE, shp2$hotspot_coverage, method="p")

plot_grid(
ggplot(shp2, aes(x=hotspot_coverage, y=richness))+
  geom_density_2d()+
  geom_smooth()+
  annotate("text", x=0.9,y=7500, label="rho=ns"),
ggplot(shp2, aes(x=hotspot_coverage, y=WE))+
  geom_density_2d_filled(show.legend=F)+
  scale_fill_scico_d(palette="lajolla")+
  geom_smooth() + coord_cartesian(ylim=c(0,3500))+
  annotate("text", x=0.9,y=3000, label="rho=0.17"),
ggplot(shp2, aes(x=hotspot_coverage, y=SES.PD))+
  geom_rug()+
  geom_smooth()+
  annotate("text", x=0.9,y=-10, label="rho=0.26"),
ggplot(shp2, aes(x=hotspot_coverage, y=SES.PE))+
  geom_point(aes(col=SES.PE), show.legend=F)+  scale_color_scico(palette="lajolla")+
  geom_smooth()+
  annotate("text", x=0.9,y=25, label="rho=0.34"), 
nrow=2)

## Lee's L ----
tmp <- shp2[,grep("LEVEL|hotspot_coverage|SES.PE$|SES.PD$|WE|richness", names(shp2))]
nb <- spdep::poly2nb(tmp, row.names = tmp$LEVEL3_COD)
col.W <- spdep::nb2listw(nb, style="W", zero.policy = TRUE)
tmp <- st_drop_geometry(tmp)
leeHS <- apply(tmp[,grep("SES.PD|SES.PE|WE|richness", names(tmp))], 2, spdep::lee.test, y=tmp$hotspot_coverage, 
               listw=col.W, zero.policy = TRUE, alternative="two.sided", na.action=na.omit)
leeHS.df <- sapply(leeHS, "[[", "estimate")
leeHS.df <- rbind(leeHS.df, sapply(leeHS, "[[", "p.value"))
row.names(leeHS.df) <- c("Lee", "expect", "var", "pvalue")
leeHS.df <- as.data.frame(t(leeHS.df))
leeHS.df$twoSD <- 2*sqrt(leeHS.df$var)

ggplot(leeHS.df, aes(y=expect, x=row.names(leeHS.df)))+
  geom_linerange(aes(ymin=expect-twoSD, ymax=expect+twoSD), col="grey70")+
  geom_point(aes(y=Lee, x=row.names(leeHS.df), col=factor(pvalue<0.05)), show.legend = F)+
  scale_color_manual("p<0.05",values=c("#e4974a"))+
  ylab("Lee's L with hotspot coverage")+
  xlab("")+  coord_flip()
ggsave("figures/LeesL_HS_coverage_and_PD_vars.png", width=3, height=2, units = "in", dpi = 300)

save.image("workspace_point1.RData")














# Choropleth threats + PE / PD -------------------------------
load("workspace_point1.RData")

#### plot setup ----
my_pal1 <- c("#d3d3d3", "#a8aec3", "#7f88b2", "#5665a3", 
             "#d3bc9a", "#a89b8e", "#7f7a82", "#565a77", 
             "#d3a45e", "#a88756", "#7f6a4f", "#564e48", 
             "#d38819", "#a87017", "#7f5816", "#564114")
my_pal1 <- c(
  "1-1" = as.character(my_pal1[1]), # low x, low y, etc.....
  "2-1" = as.character(my_pal1[2]), 
  "3-1" = as.character(my_pal1[3]), 
  "4-1" = as.character(my_pal1[4]),
  "1-2" = as.character(my_pal1[5]),
  "2-2" = as.character(my_pal1[6]),
  "3-2" = as.character(my_pal1[7]),
  "4-2" = as.character(my_pal1[8]), 
  "1-3" = as.character(my_pal1[9]), 
  "2-3" = as.character(my_pal1[10]), 
  "3-3" = as.character(my_pal1[11]), 
  "4-3" = as.character(my_pal1[12]), 
  "1-4" = as.character(my_pal1[13]), 
  "2-4" = as.character(my_pal1[14]), 
  "3-4" = as.character(my_pal1[15]), 
  "4-4" = as.character(my_pal1[16]))
my_pal <- my_pal1



#### bi_classes ----
dim <- 4
# PE_deforest_bi <- bi_wrap(dat=shp2, x=deforestation2, y=SES.PE, style="jenks", dim=4)

shp2$bi_class_PE_deforest <- bi_class(shp2, x = deforestation2, y = SES.PE, style = "jenks", dim = 4)$bi_class
shp2$bi_class_PD_deforest <- bi_class(shp2, x = deforestation2, y = SES.PD, style = "jenks", dim = 4)$bi_class
shp2$bi_class_PE_hfp <- bi_class(shp2, x = hfp_mean, y = SES.PE, style = "jenks", dim = 4)$bi_class
shp2$bi_class_PD_hfp <- bi_class(shp2, x = hfp_mean, y = SES.PD, style = "jenks", dim = 4)$bi_class
# pre_change fix: that one neg value fucks up clusters, move closer to the rest
shp2$pre_change[which.min(shp2$pre_change)] <- -2000
shp2$bi_class_PE_pre_change <- bi_class(shp2, x=pre_change, y=SES.PE, style="jenks", dim=4)$bi_class
shp2$bi_class_PD_pre_change <- bi_class(shp2, x=pre_change, y=SES.PD, style="jenks", dim=4)$bi_class
shp2$bi_class_PE_mat_change <- bi_class(shp2, x=mat_change, y=SES.PE, style="jenks", dim=4)$bi_class
shp2$bi_class_PD_mat_change <- bi_class(shp2, x=mat_change, y=SES.PD, style="jenks", dim=4)$bi_class
thicc_lines <- shp2[which(shp2$area<min.area),]


#### facet wrap ----
df <- shp2[, grepl("bi_class|SES\\.PD$|SES\\.PE$|area$", names(shp2))]
df <- df[,!grepl("bi_class$|bi_class_ses$", names(df))]  # remove unwanted columns
df.mlt <- tidyr::pivot_longer(df, cols=grep("bi_class", names(df)), names_to="group")
thicc.mlt <- df.mlt[which(df.mlt$area<min.area),]

fw <- ggplot() +
  geom_sf(df.mlt, mapping = aes(fill=value), color = NA, size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal=my_pal, dim=dim, na.value="white") +
  geom_sf(data=thicc.mlt, lwd=1, aes(col=value), show.legend=F)+
  bi_scale_color(pal=my_pal, dim=dim, na.value="white")+
  coord_sf(expand=F)+
  facet_wrap(~group, ncol=2)+
  theme(strip.background=element_blank())

#### legends ----
legs <- unique(df.mlt$group)
legs.names <- unique(df.mlt$group)
na <- gsub("bi_class_", "", legs)
na.x <- gsub("^.*?_", "", na)
na.y <- paste0("SES.", gsub("_.*", "", na))
for(i in 1:8){
  tmp <- bi_legend(pal=my_pal, dim=dim, xlab=na.x[i], ylab=na.y[i], size=5)+
    theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text = element_blank(), legend.background=element_blank())+
    annotate("text", x=rep(1:4,each=4), y=rep(1:4,4), 
             label=class_col(df.mlt$value[df.mlt$group==legs[i]]), col="white", size=2, alpha=.7)
  assign(legs.names[i], tmp)
}

## maps ----
ggdraw() + draw_plot(fw, 0, 0, 1, 1) + draw_plot(bi_class_PD_deforest, 0.095, .75, .09, .09)+
  draw_plot(bi_class_PD_mat_change, 0.095, .51, .09, .09) + 
  draw_plot(bi_class_PE_deforest, 0.095, .265, .09, .09) + 
  draw_plot(bi_class_PE_mat_change, 0.095, .02, .09, .09) + 
  draw_plot(bi_class_PD_hfp, 0.565, .75, .09, .09) + 
  draw_plot(bi_class_PD_pre_change, 0.565, .51, .09, .09) + 
  draw_plot(bi_class_PE_hfp, 0.565, .265, .09, .09) +
  draw_plot(bi_class_PE_pre_change, 0.565, .02, .09, .09) 

ggsave("figures/facet_maps.png", height=8, width=6.3, unit="in", dpi=300)

save.image("workspace_point2.RData")












# Our hotspots  ---------------------------
load("workspace_point2.RData")

dim <- 4
shp2$PE_hotspot <- bi_class(shp2, x = WE, y = SES.PE, style = "jenks", dim = dim)$bi_class
shp2$PD_hotspot <- bi_class(shp2, x = richness, y = SES.PD, style = "jenks", dim = dim)$bi_class
thicc_lines <- shp2[which(shp2$area<min.area),]

df <- shp2[, grepl("PE_hotspot|PD_hotspot|SES\\.PD$|SES\\.PE$|area$", names(shp2))]
df.mlt <- tidyr::pivot_longer(df, cols=grep("hotspot", names(df)), names_to="group")
thicc.mlt <- df.mlt[which(df.mlt$area<min.area),]

hs <- ggplot() +
  geom_sf(df.mlt, mapping = aes(fill=value), color = NA, size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal=my_pal, dim=dim, na.value="white") +
  geom_sf(data=thicc.mlt, lwd=1, aes(col=value), show.legend=F)+
  bi_scale_color(pal=my_pal, dim=dim, na.value="white")+
  coord_sf(expand=F)+
  facet_wrap(~group, ncol=2)+
  theme(strip.background=element_blank())+
  

#### create legends ----
legs <- unique(df.mlt$group)
legs.names <- unique(df.mlt$group)
na.x <- c("WE", "SR")
na.y <- c("SES.PE", "SES.PD")
for(i in 1:length(legs)){
  tmp <- bi_legend(pal=my_pal, dim=dim, xlab=na.x[i], ylab=na.y[i], size=5)+
    theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text = element_blank(), legend.background=element_blank())+
    annotate("text", x=rep(1:4,each=4), y=rep(1:4,4), 
             label=class_col(df.mlt$value[df.mlt$group==legs[i]]), col="white", size=2, alpha=.7)
  assign(legs.names[i], tmp)
}

## draw maps ----
ggdraw() + draw_plot(hs, 0, 0, 1, 1) + 
  draw_plot(PD_hotspot, 0.095, .385, .09, .09)+
  draw_plot(PE_hotspot, 0.565, .385, .09, .09)+
  draw_figure_label("A", "top.left")+
  draw_figure_label("B", "top")
  
ggsave(paste0("figures/choropleth_hotspots.png"), height=4, width=6.3, unit="in", dpi=300)



## TABLE 1 -------------------------------------------------
table(shp2$PD_hotspot)
table(shp2$PE_hotspot)

table(shp2$hotspot_coverage>0)

pdspots <- "3-3"
pespots <- c("3-4", "4-3", "4-4")
table(shp2$PD_hotspot==pdspots, shp2$hotspot_coverage>0)
table(shp2$PE_hotspot%in%pespots, shp2$hotspot_coverage>0)

tabs <- st_drop_geometry(shp2)

tabs[tabs$PD_hotspot==pdspots,c("LEVEL_NAME", "hotspot_coverage")]
tabs[tabs$PE_hotspot%in%pespots,c("LEVEL_NAME", "hotspot_coverage")]
tmp <- apply(tabs[,c("deforestation2", "hfp_mean", "mat_change", "pre_change")], 2, function(x){rank(x)/length(x)})
colnames(tmp) <- paste0(colnames(tmp), "_rank")
tabs <- cbind(tabs, tmp)

(tab1 <- tabs[tabs$PD_hotspot==pdspots|tabs$PE_hotspot%in%pespots,c("LEVEL_NAME", "richness", "deforestation2",
                                "hfp_mean", "mat_change", "pre_change",
                                "hotspot_coverage")])
tab1[,-1] <- round(tab1[,-1], 2)
tab1 <- tab1[order(tab1$LEVEL_NAME),]
knitr::kable(tab1, digits=2, format="simple", row.names=F)
knitr::kable(tab1, digits=2, format="html", row.names=F) %>% kableExtra::kable_styling("hover")
library(DT)
datatable(tab1)

## get piecharts for biomes
biomes <- readRDS("../DATA/PDiv/biomes_olson_ALL.rds")
biome_names <- c("(Sub)Tropical moist broadleaf forests", 
                 "(Sub)Tropical dry broadleaf forests", 
                 "(Sub)Tropical coniferous forests",
                 "Temperate broadleaf and mixed forests",
                 "Temperate coniferous forests", "Boreal_forests/taiga",
                 "(Sub)Tropical grasslands, savannas, shrublands",
                 "Temperate grasslands, savannas, shrublands",
                 "Flooded_grasslands and savannas",
                 "Montane grasslands and shrublands", "Tundra",
                 "Mediterranean forests, woodlands, scrub",
                 "Deserts and xeric shrublands", "Mangroves", "Lakes", "Rock and ice")

names(biomes)[grepl("^X", names(biomes))] <- biome_names

# remove everything less than 3%
for(i in 1:nrow(biomes)){
  tmp <- biomes[i,-1]
  tmp[which(tmp<0.03)] <- 0
  biomes[i,-1] <- tmp
}

rowSums(biomes[,-1], na.rm=T)
# scale to 100% to account for minor boundary inaccuracies
biomes[,-1] <- t(apply(biomes[,-1], 1, function(x){x/sum(x, na.rm=T)}) )
rowSums(biomes[,-1], na.rm=T)

hab.df <- tabs[tabs$PD_hotspot%in%pdspots | shp2$PE_hotspot%in%pespots, c("LEVEL_NAME", "LEVEL3_COD")]
hab.df <- merge(hab.df, biomes, all.x=TRUE, by.x="LEVEL3_COD", by.y="country")
hab.mlt <- melt(setDT(hab.df))
hab.mlt <- hab.mlt[order(LEVEL_NAME, value),]
hab.mlt$sort <- rep(16:1, length(unique(hab.mlt$LEVEL_NAME)))
# remove empty biomes
hab.mlt <- hab.mlt[hab.mlt$value!=0,]

#hab.mlt <- hab.mlt[hab.mlt$sort<=3,]
hab.mlt <- droplevels(hab.mlt)
hab.mlt <- as.data.frame(hab.mlt)

# ggplot(hab.mlt, aes(x="", y=value, fill=variable)) +
#   geom_bar(stat="identity", width=.1) +
#   coord_polar("y", start=0)+
#   scale_fill_scico_d("biome", palette="batlow", begin=.1, end=.6, alpha=.7)+
#   facet_wrap(~LEVEL_NAME)+
#   theme_void() + theme(strip.background=element_blank())+
# 
# scico_palette_show("batlow")

# treemap
library(treemapify)
ggplot(hab.mlt, aes(area=value, fill=variable)) +
  geom_treemap()+
  scale_fill_scico_d("biome", palette="batlow", begin=0, end=1, alpha=.7)+
  facet_wrap(~LEVEL_NAME, ncol=1, strip.position="left")+
  theme_void() + theme(strip.background=element_blank(), strip.text.y.left=element_text(angle=0, hjust=1))
ggsave()  





save.image("workspace_point3.RData")














# Threats + hotspots -----------------------------------------------------
load("workspace_point3.RData")

thr <- shp2[,grep("LEVEL|SES|hfp_mean|deforestation2|^pre_change|^mat_change", names(shp2))]
nb <- spdep::poly2nb(thr, row.names = thr$LEVEL3_COD)
col.W <- spdep::nb2listw(nb, style="W", zero.policy = TRUE)
thr <- st_drop_geometry(thr)

tesPD <- apply(thr[,grep("hfp|deforest|_change", names(thr))], 2, spdep::lee.test, y=thr$SES.PD, 
             listw=col.W, zero.policy = TRUE, alternative="two.sided", na.action=na.omit)
tesPD.df <- sapply(tesPD, "[[", "estimate")
tesPD.df <- rbind(tesPD.df, sapply(tesPD, "[[", "p.value"))
row.names(tesPD.df) <- c("Lee", "expect", "var", "pvalue")
tesPD.df <- as.data.frame(t(tesPD.df))
tesPD.df$twoSD <- 2*sqrt(tesPD.df$var)


tesPE <- apply(thr[,grep("hfp|deforest|_change", names(thr))], 2, spdep::lee.test, y=thr$SES.PE, 
               listw=col.W, zero.policy = TRUE, alternative="two.sided", na.action=na.omit)
tesPE.df <- sapply(tesPE, "[[", "estimate")
tesPE.df <- rbind(tesPE.df, sapply(tesPE, "[[", "p.value"))
row.names(tesPE.df) <- c("Lee", "expect", "var", "pvalue")
tesPE.df <- as.data.frame(t(tesPE.df))
tesPE.df$twoSD <- 2*sqrt(tesPE.df$var)

tesPD.df$response <- "SES.PD"; tesPE.df$response <- "SES.PE";
tes <- rbind(tesPD.df, tesPE.df)
tes$x <- gsub("1", "", row.names(tes))

# plots
ggplot(tes, aes(y=expect, x=x))+
  geom_linerange(aes(ymin=expect-twoSD, ymax=expect+twoSD), col="grey70")+
  geom_point(aes(y=Lee, x=x, col=factor(pvalue<0.05)), show.legend = F)+
  scale_color_scico_d("p<0.05", end=.8)+
  ylab("Lee's L")+
  xlab("")+  coord_flip() + facet_wrap(~response) + theme(strip.background=element_blank())

ggsave("figures/LeesL_threats.png", width=3, height=2, units = "in", dpi = 300)
#A positive Lee’s L indicates that clusters match for the two variables. A
#negative value indicates that the clusters have an opposite spatial
#distribution . A value around zero indicates that the spatialstructures of the
#two variables do not match. 

## one more trivariate map test (Brian)
# ggplot(shp2) + 
#   theme_void()+
#   geom_sf(aes(fill="FFE400", alpha=normalized(SES.PE)),lwd=0, col=NA)+  
#  geom_sf(aes(fill="FF0090", alpha=normalized(SES.PD)),lwd=0, col=NA)
# #   geom_sf(aes(fill="00ACE9", alpha=normalized(sub_trop_mbf)),lwd=0, col=NA)+
# # geom_sf(aes(fill="FF0090", alpha=normalized(SES.PE)),lwd=0, col=NA) 

blend.colors <- function(col) {
  m <- col2rgb(col, alpha=TRUE)
  vals <- apply(m, 1, mean)
  rgb(vals[1], vals[2], vals[3], vals[4], maxColorValue=255)
}
col2rgb("green", alpha=0.4)
blend.colors(c("green", "blue", "red"))

ggplot(shp2) + 
  theme_void()+
  geom_sf(aes(alpha=normalized(SES.PD)), fill="blue", lwd=0, col=NA) + 
  geom_sf(aes(alpha=normalized(SES.PE)), fill="red", lwd=0, col=NA)
# this is not ideal, the order you plot it makes a difference.... the colors are not merging
#  geom_sf(aes(alpha=normalized(sub_trop_mbf)), fill="green", lwd=0, col=NA)

##


bat <- scico(palette="batlow", n=4, begin=0.1, end=0.7, alpha=0.7)
plot(rep(1,4), col=bat, pch=21, lwd=30)
# vir <- viridis::turbo(4, begin=0.05, end=.9, alpha=.6)
# plot(rep(1,4), col=vir, pch=21, lwd=30)
deforestation_map <- ggplot(shp2) + 
  theme_void()+
    geom_sf(aes(fill=deforestation2),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=deforestation2), show.legend=F)+
  scale_color_gradient(low="grey95", high=bat[2])+
  scale_fill_gradient("deforestation", low="grey95", high=bat[2])+
  theme(legend.position = c(0.2, 0.3),
        legend.key.height = unit(6,"mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 8))+
    xlab(" ")+ coord_sf(expand = F)
deforestation_pd <- deforestation_map + 
  geom_sf_pattern(data=shp2[shp2$PD_hotspot==pdspots,], pattern = 'stripe', fill="grey30",
                  colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.005, pattern_spacing=0.02)+
  coord_sf(expand = F)
deforestation_pe <- deforestation_map + 
  geom_sf_pattern(data=shp2[shp2$PE_hotspot%in%pespots,], pattern = 'stripe', fill="grey30",
                  colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.005, pattern_spacing=0.02)+
  coord_sf(expand = F)


hfp_map <- ggplot(shp2) + 
  theme_void()+
  geom_sf(aes(fill=hfp_mean),lwd=0, col=NA) + 
  geom_sf(data=thicc_lines, lwd=1, aes(col=hfp_mean), show.legend=F)+
  scale_color_gradient(low="grey95", high=bat[3])+
  scale_fill_gradient("hfp", low="grey95", high=bat[3])+
  theme(legend.position = c(0.2, 0.3),
        legend.key.height = unit(6,"mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 8))+
  xlab(" ")+ coord_sf(expand = F)
hfp_pd  <- hfp_map + 
  geom_sf_pattern(data=shp2[shp2$PD_hotspot==pdspots,], pattern = 'stripe', fill="grey30",
                  colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.005, pattern_spacing=0.02)+
  coord_sf(expand = F)
hfp_pe  <- hfp_map + 
  geom_sf_pattern(data=shp2[shp2$PE_hotspot%in%pespots,], pattern = 'stripe', fill="grey30",
                  colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.005, pattern_spacing=0.02)+
  coord_sf(expand = F)


mat_change_map <- ggplot(shp2) + 
  theme_void()+
  geom_sf(aes(fill=mat_change),lwd=0, col=NA) + 
  geom_sf(data=thicc_lines, lwd=1, aes(col=mat_change), show.legend=F)+
  scale_color_gradient(low="grey95", high=bat[4])+
  scale_fill_gradient("mat_change", low="grey95", high=bat[4])+
  theme(legend.position = c(0.2, 0.3),
        legend.key.height = unit(6,"mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 8))+
  xlab(" ")+ coord_sf(expand = F)
mat_change_pd <- mat_change_map + 
  geom_sf_pattern(data=shp2[shp2$PD_hotspot==pdspots,], pattern = 'stripe', fill="grey30",
                  colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.005, pattern_spacing=0.02)+
  coord_sf(expand = F)
mat_change_pe <- mat_change_map + 
  geom_sf_pattern(data=shp2[shp2$PE_hotspot%in%pespots,], pattern = 'stripe', fill="grey30",
                  colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.005, pattern_spacing=0.02)+
  coord_sf(expand = F)


pre_change_map <- ggplot(shp2) + 
  theme_void()+
  geom_sf(aes(fill=pre_change),lwd=0, col=NA) + 
  geom_sf(data=thicc_lines, lwd=1, aes(col=pre_change), show.legend=F)+
  scale_color_gradient2(low="brown", high=bat[1])+
  scale_fill_gradient2("pre_change", low="brown", high=bat[1])+
  theme(legend.position = c(0.2, 0.3),
        legend.key.height = unit(6,"mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 8))+
  xlab(" ")+ coord_sf(expand = F)
pre_change_pd <- pre_change_map + 
  geom_sf_pattern(data=shp2[shp2$PD_hotspot==pdspots,], pattern = 'stripe', fill="grey30",
                  colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.005, pattern_spacing=0.02)+
  coord_sf(expand = F)
pre_change_pe <- pre_change_map + 
  geom_sf_pattern(data=shp2[shp2$PE_hotspot%in%pespots,], pattern = 'stripe', fill="grey30",
                  colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.005, pattern_spacing=0.02)+
  coord_sf(expand = F)


plot_grid(deforestation_pd, hfp_pd, mat_change_pd, pre_change_pd, 
          deforestation_pe, hfp_pe, mat_change_pe, pre_change_pe, ncol=2)
ggsave("figures/threats_with_PD_hotspots.png", width=14, height=18, units = "in", dpi = 300, bg="white")





# Our hotspots on the spectrum -----------------------
### reshape dataframe for plotting
plot.df <- st_drop_geometry(shp2[,c("LEVEL_NAME", "deforestation2", "hfp_mean", 
                                    "pre_change", "mat_change", "PD_hotspot", "PE_hotspot")])
plot.df$PD_hotspot[plot.df$PD_hotspot!=pdspots] <- "n"
plot.df$PD_hotspot[plot.df$PD_hotspot==pdspots] <- "y"
plot.df$PE_hotspot[!plot.df$PE_hotspot%in%pespots] <- "n"
plot.df$PE_hotspot[plot.df$PE_hotspot%in%pespots] <- "y"
plot.df <- reshape::melt(plot.df)
#plot.df$value <- unlist(tapply(plot.df$value, plot.df$variable, normalized))

ypos <- rep(rev(seq(0.2,0.80,0.04)), 4)

library(ggtext)

ggplot(plot.df, aes(x=value, y=..scaled..))+
  geom_density(fill="grey85", col=NA)+
  geom_vline(data=plot.df[plot.df$PD_hotspot=="y"|plot.df$PE_hotspot=="y",], aes(xintercept=value, col=LEVEL_NAME), 
             size=.5, show.legend=F)+
  geom_richtext(data=plot.df[plot.df$PD_hotspot=="y"|plot.df$PE_hotspot=="y",], 
                aes(x=value, y=ypos[rank(value)], label=LEVEL_NAME, col=LEVEL_NAME),
                size=2, angle=40, hjust=0, label.padding=unit(0,"mm"), 
                show.legend=F, label.color=NA)+
  scale_color_scico_d(end=0.9)+
  scale_y_continuous("Scaled density")+
  facet_wrap("variable", scales="free", labeller=as_labeller(c("deforestation2"="deforestation",
                                                               "hfp_mean"="HFP",
                                                               "pre_change"="future PRE change",
                                                               "mat_change"="future MAT change")))+
  coord_cartesian(expand=T, ylim=c(0.05,1.1))+
  # mean value indicators
  theme(strip.background=element_blank(), axis.text.y=element_blank(), 
        axis.ticks=element_blank(), panel.grid=element_blank())
ggsave("figures/all_hotspots_on_the_spectrum.pdf", width=7, height=5, units = "in", dpi = 300, bg="white")
#ggsave("figures/all_hotspots_on_the_spectrum.png", width=8, height=7, units = "in", dpi = 300, bg="white")

# alt visualisation
plot.df <- plot.df[order(plot.df$variable, plot.df$value),]
plot.df$x <- rep(1:353, 4)
plot.df$lab <- plot.df$LEVEL_NAME
plot.df$lab[which(plot.df$PD_hotspot!="y"|plot.df$PE_hotspot!="y")] <- NA
plot.df.sub <- plot.df[which(plot.df$PD_hotspot=="y"|plot.df$PE_hotspot=="y"),]

ggplot(plot.df, aes(y=value, x=x))+
  geom_area(fill="grey85")+
  geom_bar(data=plot.df.sub, aes(x=x, y=value, fill=LEVEL_NAME, col=LEVEL_NAME), 
           stat="identity", width=.9, show.legend=F)+
  facet_wrap(~variable, scales="free")+
  geom_richtext(data=plot.df.sub, aes(x=x, y=value, label=LEVEL_NAME, col=LEVEL_NAME),
               size=2, angle=90, hjust=0, label.padding=unit(.1,"mm"),
               show.legend=F, label.color=NA)+
  scale_color_scico_d(end=0.9)+scale_fill_scico_d(end=0.9)+
  scale_y_continuous("Value")+
  facet_wrap("variable", scales="free", 
             labeller=as_labeller(c("deforestation2"="deforestation",
                              "hfp_mean"="HFP", "pre_change"="future PRE change",
                              "mat_change"="future MAT change")))+
  coord_cartesian(expand=F, clip="off", xlim=c(0,354))+
  theme(strip.background=element_blank(), axis.text.x=element_blank(), 
        axis.ticks=element_blank(), panel.grid=element_blank())+
  ggtitle("Inverse cumulative distribution function")

#ggsave("figures/all_hotspots_on_the_spectrum_alt_layout.pdf", width=6, height=4, units = "in", dpi = 300, bg="white")

# (spectrum_PD <- ggplot(plot.df, aes(x=value, y=..scaled..))+
#   geom_density(fill="grey", col=NA)+
#   # geom_vline(data=data.frame(variable=names(tapply(plot.df$value, plot.df$variable, mean, na.rm=T)), 
#   #                            median=as.numeric(tapply(plot.df$value, plot.df$variable, mean, na.rm=T))),
#   #            aes(xintercept=median), col="white")+
#   geom_vline(data=plot.df[plot.df$PD_hotspot=="y",], aes(xintercept=value, col=LEVEL_NAME), 
#              size=1, show.legend=F)+
#   geom_richtext(data=plot.df[plot.df$PD_hotspot=="y",], 
#                 aes(x=value, y=ypos[rank(value)], label=LEVEL_NAME, col=LEVEL_NAME),
#                 size=2.5, angle=40, hjust=0, label.padding=unit(0,"mm"), 
#                 show.legend=F, label.color=NA)+
#   scale_color_scico_d(end=0.9)+
#   facet_wrap("variable", scales="free", labeller=as_labeller(c("deforestation2"="deforestation",
#                                                                   "hfp_mean"="HFP",
#                                                                   "pre_change"="future PRE change",
#                                                                   "mat_change"="future MAT change")))+
#   coord_cartesian(expand=F, ylim=c(0,1.1), xlim=c(0,1.1))+
#     # mean value indicators
#   theme(strip.background=element_blank(), axis.text=element_blank(), 
#         axis.ticks=element_blank(), panel.grid=element_blank()))
# 
# ## more facets! ----
# # library(tidytext)
# # ggplot(plot.df[order(plot.df$value),], aes(x=seq(nrow(plot.df)), y=value))+
# #   geom_line(fill="grey")+
# #   facet_wrap(~variable, scales="free")
# # 
# #   geom_vline(data=data.frame(variable=names(tapply(plot.df$value, plot.df$variable, mean, na.rm=T)), 
# #                              median=as.numeric(tapply(plot.df$value, plot.df$variable, mean, na.rm=T))),
# #              aes(xintercept=median), col="white")+
# #   # geom_vline(data=plot.df[plot.df$PD_hotspot=="y",], aes(xintercept=value, col=LEVEL_NAME), 
# #   #            size=1, show.legend=F)+
# #   geom_richtext(data=plot.df[plot.df$PD_hotspot=="y",],
# #                 aes(x=value, y=1, label=LEVEL_NAME, col=LEVEL_NAME),
# #                 size=2.5, angle=90, hjust=0, label.padding=unit(0,"mm"),
# #                 show.legend=F, label.color=NA)+
# #   facet_wrap("variable", scales="free", labeller=as_labeller(c("deforestation2"="deforestation",
# #                                                                "hfp_mean"="HFP",
# #                                                                "pre_change"="future PRE change",
# #                                                                "mat_change"="future MAT change")))+
# #   geom_rug(data=plot.df[plot.df$PD_hotspot=="y",], aes(col=LEVEL_NAME)8
# #            length=unit(5,"mm"), lwd=1, show.legend=F)+
# #     coord_cartesian(expand=F, xlim=c(0,1.05))+
# #   #coord_flip(expand=F, xlim=c(0,1.1))+
# #   # mean value indicators
# #   theme(strip.background=element_blank(), axis.text=element_blank(), 
# #         axis.ticks=element_blank(), panel.grid=element_blank())
# 
# 
# 
# 
# spectrum_PE <- ggplot(plot.df, aes(x=value, y=..scaled..))+
#   geom_density(fill="grey", col=NA)+
#   # geom_vline(data=data.frame(variable=names(tapply(plot.df$value, plot.df$variable, mean, na.rm=T)), 
#   #                            median=as.numeric(tapply(plot.df$value, plot.df$variable, mean, na.rm=T))),
#   #            aes(xintercept=median), col="white")+
#   geom_vline(data=plot.df[plot.df$PE_hotspot=="y",], aes(xintercept=value, col=LEVEL_NAME), 
#              size=1, show.legend=F)+
#   geom_richtext(data=plot.df[plot.df$PE_hotspot=="y",], 
#                 aes(x=value, y=ypos[rank(value)], label=LEVEL_NAME, col=LEVEL_NAME),
#                 size=2.5, angle=40, hjust=0, label.padding=unit(0,"mm"), 
#                 show.legend=F, label.color=NA)+
#   facet_wrap("variable", scales="free", labeller=as_labeller(c("deforestation2"="deforestation",
#                                                                "hfp_mean"="HFP",
#                                                                "pre_change"="future PRE change",
#                                                                "mat_change"="future MAT change")))+
#   scale_color_scico_d(end=0.9)+
#   coord_cartesian(expand=F, ylim=c(0, 1.1))+
#   #   # PD hotspot indicators
#   # geom_rug(data=plot.df[plot.df$PD_hotspot=="y",], aes(col=LEVEL_NAME), 
#   #          length=unit(10,"mm"), lwd=2, lty=2)+
#   theme(strip.background=element_blank(), axis.text=element_blank(), 
#         axis.ticks=element_blank(), panel.grid=element_blank())
# 
# plot_grid(spectrum_PD, spectrum_PE, nrow=2, labels=c("PD", "PE"), label_size=10, label_fontface=1)
ggsave("figures/hotspots_on_the_spectrum.pdf", width=8, height=7, units = "in", dpi = 300, bg="white")
ggsave("figures/hotspots_on_the_spectrum.png", width=8, height=7, units = "in", dpi = 300, bg="white")







# Our hotspots with Myer ------------

## adding Myer hotspots -------------------------------------
ha <- st_read("hotspot_area.gpkg")
pm <- st_read("Poly-Micronesia.gpkg")
gc()

# ha_united <- st_union(ha)
# st_write(ha_united, "hotspot_area_united.gpkg")
ha_united <- st_read("hotspot_area_united.gpkg")

# change projection + solve dateline issue
ha_united <- st_wrap_dateline(ha_united, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
pm <- st_wrap_dateline(pm, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
ha_united <- st_transform(ha_united, my_projection)
pm <- st_transform(pm, my_projection)


## Ordination of threats ---------
dat <- st_drop_geometry(shp2[,c("deforestation2", "hfp_mean", "pre_change", "mat_change")])
dat <- na.omit(dat)
tmp <- apply(dat, 2, normalized)
par(mfrow=c(2,2))
apply(tmp, 2, hist)
dat$deforestation2 <- sqrt(dat$deforestation2)
dat$hfp_mean <- sqrt(dat$hfp_mean)
tmp <- apply(dat, 2, normalized)

pca1 <- prcomp(tmp, center = TRUE, scale. = TRUE)
summary(pca1)
biplot(pca1, cex=0.9)

## Kernel PCA (non-linear pca)
library(kernlab)
kpca1 <- kpca(~., data = as.data.frame(tmp), kernel = 'rbfdot')
pcv(kpca1)
plot(rotated(kpca1),
     xlab="1st Principal Component",ylab="2nd Principal Component")


st_bbox(shp2)
# (PD_hotspot_bi <- ggplot() +
#   geom_sf(shp2, mapping = aes(fill=PD_hotspot==pdspots), color="black", size=0.1, show.legend=FALSE) +
#   geom_sf(data=thicc_lines, lwd=1, col="grey95", show.legend=F)+
#   scale_fill_manual(values=c(NA, "red"), na.value="white")+
#   geom_sf(data=ha_united, col=NA, fill="grey", alpha=.7)+
#   theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm")) + 
#   coord_sf(expand = F, ylim=c(-6031217, 5558740), xlim=c(-10092430, 12592430)) + theme_void())
# PD_hotspot_map <- PD_hotspot_bi
# ggsave("figures/choropleth_PD_richness_bi.png", width=10, height=5.74, units = "in", dpi = 600, bg = "white")
# 
# (PE_hotspot_bi <- ggplot() +
#   geom_sf(shp2, mapping = aes(fill = PE_hotspot%in%pespots), color = "black", size = 0.1, show.legend = FALSE) +
#   geom_sf(data=thicc_lines, lwd=1, col="grey95", show.legend=F)+
#   scale_fill_manual(values=c(NA, "red"), na.value="white")+
#   geom_sf(data=ha_united, col=NA, fill="grey", alpha=.7)+
#   # geom_sf_pattern(data=ha_united, pattern = 'stripe', fill="grey30",
#   #                 colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.005, pattern_spacing=0.02)+
#   # geom_sf_pattern(data=pm, pattern = 'stripe', fill="grey30",
#   #                 colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.01, pattern_spacing=0.02)+
#   theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm")) + 
#   coord_sf(expand = F, ylim=c(-6031217, 5558740), xlim=c(-10092430, 12592430)) + theme_void())
# PE_hotspot_map <- PE_hotspot_bi
# ggsave("figures/choropleth_PE_richness_bi.png", width=10, height=5.74, units = "in", dpi = 600, bg = "white")
# 
# plot_grid(PD_hotspot_map, PE_hotspot_map, ncol=2, labels=c("A", "B"), label_fontface=1)
# ggsave("figures/our_hotspots.png", width=12.5, height=3.5, units = "in", dpi = 600, bg = "white")
# ggsave("figures/our_hotspots.pdf", width=12.5, height=3.5, units = "in", dpi = 600, bg = "white")

# one map to rule them all
shp2$hotspot_type <- NA
shp2$hotspot_type[shp2$PD_hotspot==pdspots] <- "PD"
shp2$hotspot_type[shp2$PE_hotspot%in%pespots] <- "PE"
shp2$hotspot_type[shp2$PE_hotspot%in%pespots&shp2$PD_hotspot==pdspots] <- "PD & PE"

bc <- c("#52548D", "#C57391", "#EFB984")
(hotspots_bi <- ggplot() +
    geom_sf(shp2, mapping = aes(fill=hotspot_type), color=NA, size=0.1) +
    geom_sf(data=thicc_lines, lwd=1, col="grey95", show.legend=F)+
    scale_fill_manual("Hotspot\ntype", values=bc, na.value="grey95")+
    geom_sf(data=ha_united, col=NA, fill="grey", alpha=.6)+
    theme_void()+ coord_sf(expand = F)+
    theme(legend.position = c(0.2, 0.3),
      legend.key.height = unit(5,"mm"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.background = element_blank()))
ggsave(plot=hotspots_bi, "figures/one_hotspots_map.pdf", width=6.75, height=4, dpi=300, bg="white")  




# Get values for tiny bars
# remember to not make the new subset the max values, draws a wrong picture (eg
# Cape povinces having hfp of 100%)

plot.df <- st_drop_geometry(shp2[,c("LEVEL_NAME", "deforestation2", "hfp_mean", 
                                    "pre_change", "mat_change", "PD_hotspot", "PE_hotspot")])
plot.df$PD_hotspot[plot.df$PD_hotspot!=pdspots] <- "n"
plot.df$PD_hotspot[plot.df$PD_hotspot==pdspots] <- "y"
plot.df$PE_hotspot[!plot.df$PE_hotspot%in%pespots] <- "n"
plot.df$PE_hotspot[plot.df$PE_hotspot%in%pespots] <- "y"
plot.hp <- plot.df[plot.df$PD_hotspot=="y" | plot.df$PE_hotspot=="y" ,]

# plot.hp$deforestation2 <- normalized(plot.hp$deforestation2)
# plot.hp$hfp_mean <- normalized(plot.hp$hfp_mean)
# plot.hp$mat_change <- normalized(plot.hp$mat_change)
# plot.hp$pre_change <- normalized(plot.hp$pre_change)
plot.hp$deforestation2 <- plot.hp$deforestation2/max(shp2$deforestation2, na.rm=T)
plot.hp$hfp_mean <- plot.hp$hfp_mean/max(shp2$hfp_mean, na.rm=T)
plot.hp$mat_change <- plot.hp$mat_change/max(shp2$mat_change, na.rm=T)
plot.hp$pre_change <- plot.hp$pre_change/max(shp2$pre_change, na.rm=T)

# plot a facet for each country showing each variable
test <- melt(plot.hp)
#x <- rep(c(1,1,2,2), 16)
#y <- rep(c(1,2,1,2), 16)
#test$ycord <- y
#test$xcord <- x
# ggplot(test, aes(x=variable, y=1, col=variable, alpha=value))+
#   geom_point(shape=15, size=12)+
#   facet_wrap(~LEVEL_NAME+PD_hotspot, scales="free")+
#   theme_void()

bat <- scico(palette="hawaii", n=4, begin=0, alpha=1)
# arrange colors (bat[2] = defore)
bat2 <- c(bat[3], bat[2], bat[4], bat[1])
ggplot(test, aes(fill=variable, x=1, y=value))+
  geom_bar(position=position_dodge(), stat="identity")+
  scale_fill_manual(values=bat2)+
  facet_wrap(~LEVEL_NAME+PD_hotspot, scales="free")+
  theme_void()+coord_cartesian(expand=F, ylim=c(-0.3, 1))+theme(panel.border=element_rect(fill=NA))
ggsave("figures/little_bars_alpha09.pdf", width=3, height=3)





save.image("workspace_maps.RData")




# MORE SI PLOTS #####
load("workspace_maps.RData")

ggplot(shp2, aes(x=SES.PD, y=SES.PE, label=LEVEL_NAME))+
  geom_point()+
  geom_label(nudge_x=1, hjust=0, label.size=0, label.padding=unit(0.5,"mm"), 
             alpha=0, color = alpha('black', .5), size=3)+
  scale_x_continuous(limits=c(-80,20))
ggsave("figures/SES_PE_vs_SES_PD.png", width=5, height=5)



+# C A C H E ------------------------------------------------------------------
# 
# 

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
#   #scale_colour_scico("SR", palette="batlow", trans = "sqrt", 
#   #                        begin = lcol, end = sqrt(ucol))+
#   scale_fill_scico("threat", palette="batlow")+ #, 
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


# # all hotspots in one map
# thicc_lines <- shp2[which(shp2$area<min.area),]
# al <- 0.3
# ggplot()+
#   # world layer
#   geom_sf(data=shp2, col=NA, fill="grey95", lwd=0, show.legend = F) +
#   # PE layer
#   geom_sf(data=shp2[shp2$SES.PE_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
#   geom_sf(data=thicc_lines[thicc_lines$SES.PE_spot==1, ], col="#6575B1", show.legend=F, lwd=2) +
#   # PD layer
#   geom_sf(data=shp2[shp2$SES.PD_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
#   geom_sf(data=thicc_lines[thicc_lines$SES.PD_spot==1, ], col=alpha("#6575B1",0.3), show.legend=F, lwd=2) +
#   # SR layer
#   geom_sf(data=shp2[shp2$SR_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
#   geom_sf(data=thicc_lines[thicc_lines$SR_spot==1, ], col="#6575B1", show.legend=F, lwd=2) +
#   # style settings
#   ggtitle("All hotspots, layered")#+
# #  theme_void()+
# #  # add Myer layer
# #  geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
# #  geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)
# 
# # get 
# x <- st_drop_geometry(shp2[,c("SR_spot","SES.PD_spot", "SES.PE_spot")])
# x$SR_spot <- ifelse(x$SR_spot<0, NA,x$SR_spot)
# x$SES.PE_spot <- ifelse(x$SES.PE_spot<0, NA,x$SES.PE_spot)
# x$SES.PD_spot <- ifelse(x$SES.PD_spot<0, NA,x$SES.PD_spot)
# shp2$hotspot_sum <- rowSums(x, na.rm = T)
# 
# ggplot()+
#   # world layer
#   geom_sf(data=shp2, aes(fill=factor(hotspot_sum)),col=NA, lwd=0) +
#   geom_sf(data=thicc_lines[which(thicc_lines$LEVEL3_COD %in% shp2$LEVEL3_COD[shp2$hotspot_sum!=0]), ],
#           col="#6575B1", show.legend=F, lwd=2) +
#   # hotspot sum laye...
#   # geom_sf(aes(fill=hotspot_sum), col=NA, lwd=0, show.legend = F, alpha=al) + 
#   # geom_sf(data=thicc_lines[thicc_lines$SES.PE_spot==1, ], col="#6575B1", show.legend=F, lwd=2) +
#   # # PD layer
#   # geom_sf(data=shp2[shp2$SES.PD_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
#   # geom_sf(data=thicc_lines[thicc_lines$SES.PD_spot==1, ], col=alpha("#6575B1",0.3), show.legend=F, lwd=2) +
#   # # SR layer
#   # geom_sf(data=shp2[shp2$SR_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
#   # geom_sf(data=thicc_lines[thicc_lines$SR_spot==1, ], col="#6575B1", show.legend=F, lwd=2) +
#   # # style settings
#   scale_fill_manual(values = c("grey95", "#6575B1"))+
#   # scale_color_manual(values = c("grey90", "#6575B1"))+
#   ggtitle("Sum all hotspots")+
#   # add Myer layer
#   geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
#   geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)
# 
# 
# 

# 
# # Hotspot definition ------------------------------------------------------
# # Myer2000 uses endemism (at least x number of endemics) and habitat loss(at least 70% habitat loss)
# 
# shp2 <- readRDS("including_hotspot_coverage.rds")
# 
# # get a reasonable percent:
# kmean_calc <- function(df, ...){
#   kmeans(df, scaled = ..., nstart = 30)
# }
# km2 <- kmean_calc(shp2$SES.PE, 2)
# km3 <- kmean_calc(shp2$SES.PE, 3)
# km4 <- kmeans(shp2$SES.PE, 4)
# km5 <- kmeans(shp2$SES.PE, 5)
# km6 <- kmeans(shp2$SES.PE, 6)
# km7 <- kmeans(shp2$SES.PE, 7)
# km8 <- kmeans(shp2$SES.PE, 8)
# km9 <- kmeans(shp2$SES.PE, 9)
# km_perc <- function(x)x$betweenss/x$totss
# plot(unlist(lapply(list(km2,km3,km4,km5,km6,km7,km8,km9),km_perc)))
# plot(shp2$SES.PE~shp2$richness, col=km6$cluster) # top 8 look good, matches 2.5% well
# 
# # highest 2.5% of everything:
# ## Endemism hotspots ---------------------
# C <- coldspots(shp2$SES.PE) # coldspots
# H <- hotspots(shp2$SES.PE) # hotspots
# ## Merge endemism values to shapefile of grid cells.
# DF <- data.frame(LEVEL3_COD=shp2$LEVEL3_COD, cold=C, hot=H)
# DF$SES.PE_spot <- 0
# DF$SES.PE_spot[DF$hot==1] <- 1
# DF$SES.PE_spot[DF$cold==1] <- -1
# shp2 <- merge(shp2, DF[,c("LEVEL3_COD", "SES.PE_spot")], by = "LEVEL3_COD", all = TRUE)
# shp2$SES.PE_spot
# 
# ggplot(shp2, aes(x=factor(SES.PE_spot), y=hotspot_coverage, label=LEVEL3_COD))+
#   #geom_point() +
#   geom_text(position="jitter")
# 
# 
# ## PD hotspots --------------------------
# C <- coldspots(shp2$SES.PD) # coldspots
# H <- hotspots(shp2$SES.PD) # hotspots
# DF <- data.frame(LEVEL3_COD=shp2$LEVEL3_COD, cold=C, hot=H)
# DF$SES.PD_spot <- 0
# DF$SES.PD_spot[DF$hot==1] <- 1
# DF$SES.PD_spot[DF$cold==1] <- -1
# shp2 <- merge(shp2, DF[,c("LEVEL3_COD", "SES.PD_spot")], by = "LEVEL3_COD", all = TRUE)
# 
# 
# 
# ## Richness hotspots -----------------------
# C <- coldspots(shp2$richness) # coldspots
# H <- hotspots(shp2$richness) # hotspots
# DF <- data.frame(LEVEL3_COD=shp2$LEVEL3_COD, cold=C, hot=H)
# DF$SR_spot <- 0
# DF$SR_spot[DF$hot==1] <- 1
# DF$SR_spot[DF$cold==1] <- -1
# shp2 <- merge(shp2, DF[,c("LEVEL3_COD", "SR_spot")], by = "LEVEL3_COD", all = TRUE)
# 
# 
# 
# ## Deforestation hotspots -----------------------
# H <- hotspots(shp2$deforestation2) # hotspots
# shp2$deforest_spot <- H
# 
# # Myers & PD hotspot matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage>0 & shp2$SES.PD_spot==1]
# # Myers & PD hotspot no matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage==0 & shp2$SES.PD_spot==1]
# 
# # Myers & PE hotspot matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage>0 & shp2$SES.PE_spot==1]
# # Myers & PE hotspot no matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage==0 & shp2$SES.PE_spot==1]
# 
# # Myers & SR hotspot matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage>0 & shp2$SR_spot==1]
# # Myers & SR hotspot no matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage==0 & shp2$SR_spot==1]
# 
# # Myers & Deforestation hotspot matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage>0 & shp2$deforest_spot==1]
# # Myers & Deforestation hotspot no matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage==0 & shp2$deforest_spot==1]
# 
# 
# sort factor to match facet wrap (currently adjusted to 3 classes)
# shp2$bi_class_ses_l <- factor(shp2$bi_class_ses, levels=c("1-3", "2-3", "3-3", "1-2", "2-2", "3-2", "1-1", "2-1" ,"3-1"))
# ggplot(shp2)+
#   geom_boxplot(aes(x=bi_class_ses_l, y=richness, col=bi_class_ses_l, fill=bi_class_ses_l), varwidth = F)+
#   scale_fill_manual(values=coropleth_palette[c(3,6,9,2,5,8,1,4,7)])+
#   scale_color_manual(values=coropleth_palette[c(3,6,9,2,5,8,1,4,7)])+
#   facet_wrap("bi_class_ses_l", ncol=3, scales="free_x")+
#   scale_x_discrete("", expand = c(0.2,0.2))+
#   theme_minimal_hgrid()+
#   scale_y_continuous(trans="sqrt", expand=c(0,0))+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#     panel.spacing.y = unit(0, "lines"),
#     axis.text = element_text(size=8),
#     legend.position = "none"
#     )
# 
# ## Maps hotspots ---------------------------
# min.area <- 1.5e+9
# thicc_lines <- shp2[which(shp2$area<min.area),]
# 
# # extract overlapping parts for another color
# hf <- st_read("hotspots_fixed.gpkg")
# wrld_wrap <- st_wrap_dateline(hf, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
# m <- st_transform(wrld_wrap, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
# # fix sf trouble...
# #sf::sf_use_s2(FALSE)
# #PE.int <- st_intersection(shp2[shp2$SES.PE_spot==1,], m)
# # subtract countries from hotspot areas to not cover them
# #hs.int <- st_(m, shp2[shp2$SES.PE_spot==1,])
# 
# 
# # put myers map on top
# thicc_lines$SES.PE_spot[thicc_lines$SES.PE_spot!=0] # no hot/cold zones in small areas
# (PE_hotspot_map <- ggplot(shp2)+
#     geom_sf(aes(fill=factor(SES.PE_spot==1)),lwd=.1, col=NA, show.legend = F) + 
#     geom_sf(data=thicc_lines, aes(col=factor(SES.PE_spot==1)), show.legend=F, lwd=2)+
#     scale_fill_manual(values = c("grey90", "#6575B1"))+
#     scale_color_manual(values = c("grey90", "#6575B1"))+
#     theme_void()+
#     # add Myer layer
#     geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
#     geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)+
#     # add hotspot intersection area
#     #geom_sf(data=PE.int, fill="#B16575", color=NA, show.legend=F)+  
#     # add text
#     geom_hline(yintercept=c(-2343636,0,2343636), size=0.1, lty=c("dashed", "solid","dashed"))+
#     annotate("text", x= -17067530, y= 2043636, label= "Tropic of Cancer", hjust=0)+
#     annotate("text", x= -17067530, y= -330000, label= "Equator", hjust=0)+
#     annotate("text", x= -17067530, y= -2643636, label= "Tropic of Capricorn", hjust=0)+
#     ggtitle("Top SES.PE 2.5% with Myers biodiv hotspots"))
# 
# 
# PD_hotspot_map <- ggplot(shp2)+
#   geom_sf(aes(fill=factor(SES.PD_spot==1)),lwd=0, col=NA, show.legend = F) + 
#   geom_sf(data=thicc_lines, aes(col=factor(SES.PD_spot==1)), show.legend=F, lwd=2)+
#   scale_fill_manual(values = c("grey90", "red"))+
#   scale_color_manual(values = c("grey90", "red"))+
#   theme(panel.border = element_blank())+
#   ggtitle("SES.PD Hotspots")+
#   theme_void()+
#   # add Myer layer
#   geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
#   geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)
# 
# table(thicc_lines$SR_spot) # only cold zones in small areas, adjust color
# SR_hotspot_map <- ggplot(shp2)+
#   geom_sf(aes(fill=factor(SR_spot==1)),lwd=0, col=NA, show.legend = F) + 
#   geom_sf(data=thicc_lines, aes(col=factor(SR_spot==1)), show.legend=F, lwd=2)+
#   scale_fill_manual(values = c("grey90", "red"))+
#   scale_color_manual(values = c("grey90"))+
#   theme(panel.border = element_blank())+
#   ggtitle("SR Hotspots")+
#   theme_void()+
#   # add Myer layer
#   geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
#   geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)
# 
# plot_grid(PE_hotspot_map, PD_hotspot_map, SR_hotspot_map, ncol=2)
# ggsave("figures/single_hotspots.png", width=10, height=6, units = "in", dpi = 600, bg = "white")
# 
# 
# table(thicc_lines$deforest_spot)
# (deforest_hotspot_map <- ggplot(shp2)+
#     geom_sf(aes(fill=factor(deforest_spot)),lwd=0, col=NA) + 
#     geom_sf(data=thicc_lines, aes(col=factor(deforest_spot)), show.legend=F, lwd=2)+
#     scale_fill_manual(values = c("grey90", "#6575B1"))+
#     scale_color_manual(values = c("grey90", "#6575B1"))+
#     theme(panel.border = element_blank())+
#     ggtitle("Top 2.5% deforestation countries")+
#     theme_void()+
#     # add Myer layer
#     geom_sf(data=m[m$Type=="hotspot area",], fill="#75B165", color=NA, show.legend=F, alpha=.5)+
#     geom_sf(data=m[m$NAME=="Polynesia-Micronesia" & m$Type=="outer limit",], fill="#75B165", alpha=.5, color=NA, show.legend=F)
# )

# bonbon_area <- ggplot(df2[grep("PD", df2$name),], aes(name, value, color = area, group=LEVEL_NAME))+
#   geom_bump(size = 1, position = "identity", smooth=sm) +
#   scale_color_scico(trans="sqrt", alpha = alevel, palette="batlow", guide = "none")+
#   scale_x_discrete("", expand = c(0.05,0.05))+
#   scale_y_discrete("rank")+
#   geom_point(size = 0, aes(fill=area)) +
#   scale_fill_scico("Area", trans="sqrt", alpha = 1, palette="batlow")+
#   theme_minimal()


# bonbon_area_PE <- ggplot(df2[grep("PE", df2$name),], aes(name, value, color = area, group=LEVEL_NAME))+
#   geom_bump(size = 2, position = "identity", smooth=sm) +
#   scale_color_scico(trans="sqrt", alpha = alevel, palette="batlow", guide = "none")+
#   scale_x_discrete("", expand = c(0.05,0.05))+
#   scale_y_discrete("rank")+
#   geom_point(size = 0, aes(fill=area)) +
#   scale_fill_scico("Area", trans="sqrt", alpha = 1, palette="batlow")+
#   theme_minimal()s


# lsize = 1
# font.size=2.5
# library(ggtext)  
# a <- shp2$deforestation2[shp2$PD_hotspot=="3-3"]
# b <- rev(seq(45,90,5))
# (deforst_hist_PD <- ggplot(shp2, aes(x=deforestation2))+
#     geom_histogram()+
#     geom_vline( 
#       aes(xintercept=deforestation2, col=LEVEL_3_CO), lty=1, size=lsize,show.legend=F)+
#     geom_richtext(data=data.frame(a=a, region=shp2$LEVEL_NAME[shp2$PD_hotspot=="3-3"]), 
#                   aes(x=a+0.5, y=rev(seq(20,65,5))[rank(a)], label=region, col=region), size=font.size, angle=40, hjust=0, 
#                   label.padding=unit(0,"mm"), show.legend=F, label.color=NA)+
#     coord_cartesian(expand=F, ylim=c(0,100))
# )
# a <- shp2$hfp_mean[shp2$PD_hotspot=="3-3"]
# b <- c(30,30,rev(seq(22,36,2)))
# (hfp_hist_PD <- ggplot(shp2, aes(x=hfp_mean))+
#     geom_histogram()+
#     theme(axis.title.y.left=element_blank())+
#     geom_vline(data=data.frame(a=a, region=shp2$LEVEL_3_CO[shp2$PD_hotspot=="3-3"]), 
#                aes(xintercept=a, col=region), lty=1, size=lsize)+ 
#     coord_cartesian(expand=F, ylim=c(0,40))+
#     theme(legend.key.height=unit(0.4,"cm")))
# a <- shp2$mat_change[shp2$PD_hotspot=="3-3"]
# b <- rev(seq(40,67,3))
# (mat_hist_PD <- ggplot(shp2, aes(x=mat_change))+
#     geom_histogram()+
#     theme(axis.title.y.left=element_blank())+
#     geom_vline(data=data.frame(a=a, region=shp2$LEVEL_3_CO[shp2$PD_hotspot=="3-3"]), 
#                aes(xintercept=a, col=region), lty=1, show.legend=F)+
#     coord_cartesian(expand=F, ylim=c(0,70)))
# a <- shp2$pre_change[shp2$PD_hotspot=="3-3"]
# b <- rev(seq(55,100,5))
# (pre_hist_PD <- ggplot(shp2, aes(x=pre_change))+
#     geom_histogram()+
#     theme(axis.title.y.left=element_blank())+
#     coord_cartesian(expand=F, ylim=c(0,120))+
#     geom_vline(data=data.frame(a=a, region=shp2$LEVEL_3_CO[shp2$PD_hotspot=="3-3"]), 
#                aes(xintercept=a, col=region), lty=1, show.legend=F))
# 
# a <- shp2$deforestation2[shp2$PE_hotspot%in%pespots]
# b <- rev(seq(45,90,5))
# (deforst_hist_PE <- ggplot(shp2, aes(x=deforestation2))+
#     geom_histogram()+
#     coord_cartesian(expand=F, ylim=c(0,110))+
#     geom_vline(data=data.frame(a=a, region=shp2$LEVEL_3_CO[shp2$PE_hotspot%in%pespots]), 
#                aes(xintercept=a, col=region), lty=1, show.legend=F))
# a <- shp2$hfp_mean[shp2$PE_hotspot%in%pespots]
# b <- rev(seq(18,32,1))
# (hfp_hist_PE <- ggplot(shp2, aes(x=hfp_mean))+
#     geom_histogram()+
#     theme(axis.title.y.left=element_blank())+
#     coord_cartesian(expand=F, ylim=c(0,45))+
#     geom_vline(data=data.frame(a=a, region=shp2$LEVEL_3_CO[shp2$PE_hotspot%in%pespots]), 
#                aes(xintercept=a, col=region), lty=1, show.legend=T)+
#     theme(legend.key.height=unit(0.4,"cm")))
# a <- shp2$mat_change[shp2$PE_hotspot%in%pespots]
# b <- rev(seq(46,55,1))
# (mat_hist_PE <- ggplot(shp2, aes(x=mat_change))+
#     geom_histogram()+
#     theme(axis.title.y.left=element_blank())+
#     coord_cartesian(expand=F, ylim=c(0,70))+
#     geom_vline(data=data.frame(a=a, region=shp2$LEVEL_3_CO[shp2$PE_hotspot%in%pespots]), 
#                aes(xintercept=a, col=region), lty=1, show.legend=F))
# a <- shp2$pre_change[shp2$PE_hotspot%in%pespots]
# b <- rev(seq(55,100,5))
# (pre_hist_PE <- ggplot(shp2, aes(x=pre_change))+
#     geom_histogram()+
#     theme(axis.title.y.left=element_blank())+
#     coord_cartesian(expand=F, ylim=c(0,120), xlim=c(-1000,1500))+
#     geom_vline(data=data.frame(a=a, region=shp2$LEVEL_3_CO[shp2$PE_hotspot%in%pespots]), 
#                aes(xintercept=a, col=region), lty=1, show.legend=F))
# 
# 
# 
# plot_grid(deforst_hist_PD, hfp_hist_PD, mat_hist_PD, pre_hist_PD, 
#           deforst_hist_PE, hfp_hist_PE, mat_hist_PE, pre_hist_PE, ncol=2, 
#           labels=c("PD", "", "", "", "PE", "", "", ""), label_size=9, label_fontface="plain")
#load("workspace_point1.RData")

# dim <- 4
# #breaks=bi_class_breaks(shp2, x=SES.PD, y=SES.PE, style="jenks", dim=dim)
# shp2$bi_class <- bi_class(shp2, x = richness, y = WE, style = "jenks", dim = dim)$bi_class
# shp2$bi_class_ses <- bi_class(shp2, x = SES.PD, y = SES.PE, style = "jenks", dim = dim)$bi_class
# thicc_lines <- shp2[which(shp2$area<min.area),]

# (choropl1 <- ggplot() +
#   geom_sf(na.omit(shp2), mapping = aes(fill = bi_class_ses), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = my_pal, dim=dim, na.value="white") +
#   geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class_ses), show.legend=F)+
#   bi_scale_color(pal = my_pal, dim=dim, na.value="white") +
#   bi_theme()+ theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm"))+
#   coord_sf(expand = F)) 
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="sesPD", ylab="sesPE ", size=8,
#                     breaks=bi_class_breaks(shp2, x=SES.PD, y=SES.PE, style="jenks", dim=dim))+
#   theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_text(angle=40, hjust=1),
#         axis.text=element_text(size=5))
# (pd_pe_map <- ggdraw() + draw_plot(choropl1, 0, 0, 1, 1) + draw_plot(legend, 0.05, 0, 0.3, 0.3))
# 
# 
# choropl2 <- ggplot() +
#     geom_sf(shp2, mapping = aes(fill = bi_class), color = NA, size = 0.1, show.legend = FALSE) +
#     bi_scale_fill(pal = my_pal, dim = 4, na.value="white") +
#     geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class), show.legend=F)+
#     bi_scale_color(pal = my_pal, dim = 4, na.value="white") +
#     bi_theme()+ theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm"))+
#     coord_sf(expand = F)# t,r,b,l bi_theme()+ 
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="SR", ylab="WE ", size=8,
#                       breaks=bi_class_breaks(shp2, x=richness, y=WE, style = "jenks", dim = dim))+
#     theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_text(angle=40, hjust=1),
#           axis.text=element_text(size=5))
# (sr_we_map <- ggdraw() + draw_plot(choropl2, 0, 0, 1, 1) + draw_plot(legend, 0.04, 0, 0.33, 0.33))
# 
# plot_grid(pd_pe_map, sr_we_map, ncol=2)
# ggsave("figures/choropleth_map.png", width=12.5, height=3.5, units = "in", dpi = 600, bg = "white")
# 
# 
# ## Threat boxplots for groups ------------------------------------------
# # threats: HFP, future clim change, deforestation
# 
# # boxplot for deforest vs biclass groups -a more quantitative plot
# #coropleth_palette <- read.csv("figures/coropleth_palette.txt", header = F)$V1 # actual order
# 
# (biclass_vs_deforest <- ggplot(shp2)+
#   geom_boxplot(aes(x=bi_class_ses, y=deforestation2), varwidth = T,  show.legend=F)+ #, fill=bi_class_ses
#   #scale_fill_manual(values=coropleth_palette)+
#   xlab("SES.PD-SES.PE group"))
# # it does not look like there is a strong pattern in deforestation. might be single hotspots in it though. gotta inspect visually. check for hfp and future climate change as well:
# (biclass_vs_hfp <- ggplot(shp2)+
#     geom_boxplot(aes(x=bi_class_ses, y=hfp_mean), varwidth = T,  show.legend=F)+ 
#     xlab("SES.PD-SES.PE group"))
# (biclass_vs_fut_prec <- ggplot(shp2)+
#     geom_boxplot(aes(x=bi_class_ses, y=pre_change), varwidth = T,  show.legend=F)+
#     xlab("SES.PD-SES.PE group"))
# (biclass_vs_fut_mat <- ggplot(shp2)+
#     geom_boxplot(aes(x=bi_class_ses, y=mat_change), varwidth = T,  show.legend=F)+
#     xlab("SES.PD-SES.PE group"))
# plot_grid(biclass_vs_deforest, biclass_vs_hfp, biclass_vs_fut_prec, biclass_vs_fut_mat, ncol=2)
# ggsave("figures/boxplot_PD_PE_threats.png", width=7, height=7, units = "in", dpi = 600, bg = "white")
# 
# ## continuous plot
# y <- c("deforestation2", "hfp_mean", "pre_change", "mat_change")
# for(i in 1:4){
#   response <- y[i] 
#   g <- ggplot(shp2, aes_string(x="SES.PD", y=response))+
#     geom_point()+
#     geom_smooth()
#   assign(response, g) #generate an object for each plot
# }
# for(i in 1:4){
#   response <- y[i] 
#   g <- ggplot(shp2, aes_string(x="SES.PE", y=response))+
#     geom_point()+
#     geom_smooth()
#   assign(paste0(response, "_PE"),g) #generate an object for each plot
# }
# plot_grid(deforestation2, hfp_mean, pre_change, mat_change, 
#           deforestation2_PE, hfp_mean_PE, pre_change_PE, mat_change_PE, ncol=4)
# 
# 
# 
# ## adding Myer hotspots -------------------------------------
# ha <- st_read("hotspot_area.gpkg")
# pm <- st_read("Poly-Micronesia.gpkg")
# gc()
# 
# # ha_united <- st_union(ha)
# # st_write(ha_united, "hotspot_area_united.gpkg")
# ha_united <- st_read("hotspot_area_united.gpkg")
# 
# # change projection + solve dateline issue
# ha_united <- st_wrap_dateline(ha_united, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
# pm <- st_wrap_dateline(pm, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
# ha_united <- st_transform(ha_united, my_projection)
# pm <- st_transform(pm, my_projection)
# 
# 
# # plotting with triangle pattern takes a minute, be patient!
# choropl3 <- choropl1+
#     geom_sf_pattern(data=ha_united, pattern = "circle", fill="grey30", #pattern_type="triangle", 
#                     colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.01, pattern_spacing=0.02)+
#     geom_sf_pattern(data=pm, pattern = "circle", fill="grey30", #pattern_type="triangle",  
#                     colour="grey30", lwd=0.1, pattern_size=0.2, alpha=0, pattern_density=0.01, pattern_spacing=0.02)+
#   theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm")) + coord_sf(expand = F)
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="sesPD", ylab="sesPE ", size=8,
#                     breaks=bi_class_breaks(shp2, x=SES.PD, y=SES.PE, style="jenks", dim=dim))+
#   theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_text(angle=45, hjust=1))
# pd_pe_myer <- ggdraw() + draw_plot(choropl3, 0, 0, 1, 1) + draw_plot(legend, 0.05, 0, 0.3, 0.3)
# ggsave("figures/PD_PE_Myer.png", width=10, height=5.74, units = "in", dpi = 600, bg = "white")

# OLD SINGLE CHOROLPLETH PLOTS THREATS ----
# # legend position and size
# l.x=0.07; l.y=0.01; l.width=0.25; l.height=0.25
# 
# ## PLOTS ----
# PE_deforest <- ggplot() +
#   geom_sf(shp2, mapping = aes(fill = bi_class_PE_deforest), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = my_pal, dim=dim, na.value="white") +
#   geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class_PE_deforest), show.legend=F)+
#   bi_scale_color(pal = my_pal, dim=dim, na.value="white") +
#   bi_theme()+
#   coord_sf(expand = F)
# (legend <- bi_legend(pal=my_pal, dim=dim, xlab="deforestation", ylab="sesPE ", size=8,
#                      breaks=bi_class_breaks(shp2, x=deforestation2, y=SES.PE, style="jenks", dim=dim))+
#     theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text = element_blank())+
#     annotate("text", x=rep(1:4,each=4), y=rep(1:4,4), 
#              label=class_col(shp2$bi_class_PE_deforest), col="white", size=4, alpha=.7))
# ses_pe_deforest_map <- ggdraw() + draw_plot(PE_deforest, 0, 0, 1, 1) + draw_plot(legend, l.x, l.y, l.width, l.height)
# ggsave(plot=ses_pe_deforest_map, "figures/choropleth_PE_deforest.png", 
#        width=8, height=4.62, units = "in", dpi = 300, bg = "white")
# 
# 
# PD_deforest <- ggplot() +
#   geom_sf(shp2, mapping = aes(fill = bi_class_PD_deforest), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = my_pal, dim = dim, na.value="white") +
#   geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class_PD_deforest), show.legend=F)+
#   bi_scale_color(pal = my_pal, dim = dim, na.value="white") +
#   bi_theme()+ theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm")) + coord_sf(expand = F) 
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="deforestation", ylab="SES.PD ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=deforestation2, y=SES.PD, style="jenks", dim=dim))+
#   theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_blank())+
#   annotate("text", x=rep(1:4,each=4), y=rep(1:4,4), 
#            label=class_col(shp2$bi_class_PD_deforest), col="white", size=4, alpha=.7)
# ses_pd_deforest_map <- ggdraw() + draw_plot(PD_deforest, 0, 0, 1, 1) + draw_plot(legend, l.x, l.y, l.width, l.height)
# ggsave(plot=ses_pd_deforest_map,"figures/choropleth_PD_deforest.png", width=10, height=5.74, units = "in", dpi = 600, bg = "white")
# 
# 
# 
# ### Human footprint index
# hfp_pe <- ggplot() +
#   geom_sf(shp2, mapping = aes(fill = bi_class_PE_hfp), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = my_pal, dim = dim, na.value="white") +
#   geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class_PE_hfp), show.legend=F)+
#   bi_scale_color(pal = my_pal, dim = dim, na.value="white") +
#   bi_theme()+ theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm")) + coord_sf(expand = F)
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="HFP", ylab="SES.PE ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=hfp_mean, y=SES.PE, style="jenks", dim=dim))+
#   theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_blank())+
#   annotate("text", x=rep(1:4,each=4), y=rep(1:4,4), 
#            label=class_col(shp2$bi_class_PE_hfp), col="white", size=4, alpha=.7)
# ggdraw() + draw_plot(hfp_pe, 0, 0, 1, 1) + draw_plot(legend, l.x, l.y, l.width, l.height)
# #ggsave("figures/choropleth_PE_hfp.png", width=10, height=5.74, units="in", dpi=300, bg="white")
# 
# hfp_pd <- ggplot() +
#   geom_sf(shp2, mapping = aes(fill = bi_class_PD_hfp), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = my_pal, dim = dim, na.value="white") +
#   geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class_PD_hfp), show.legend=F)+
#   bi_scale_color(pal = my_pal, dim = dim, na.value="white") +
#   bi_theme()+ theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm")) + coord_sf(expand = F)
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="HFP", ylab="SES.PD ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=hfp_mean, y=SES.PD, style="jenks", dim=dim))+
#   theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_text(angle=45, hjust=1))
# ggdraw() + draw_plot(hfp_pd, 0, 0, 1, 1) + draw_plot(legend, l.x, l.y, l.width, l.height)
# ggsave("figures/choropleth_PD_hfp.png", width=10, height=5.74, units="in", dpi=300, bg="white")
# 
# ### Future climate
# prechange_pe <- ggplot() +
#   geom_sf(shp2, mapping = aes(fill = bi_class_PE_pre_change), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = my_pal, dim = dim, na.value="white") +
#   geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class_PE_pre_change), show.legend=F)+
#   bi_scale_color(pal = my_pal, dim = dim, na.value="white") +
#   bi_theme()+ theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm")) + coord_sf(expand = F)
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="pre_change", ylab="SES.PE ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=pre_change, y=SES.PE, style="jenks", dim=dim))+
#   theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_text(angle=45, hjust=1))
# ggdraw() + draw_plot(prechange_pe, 0, 0, 1, 1) + draw_plot(legend, l.x, l.y, l.width, l.height)
# ggsave("figures/choropleth_PE_prechange.png", width=10, height=5.74, units="in", dpi=300, bg="white")
# 
# prechange_pd <- ggplot() +
#   geom_sf(shp2, mapping = aes(fill = bi_class_PD_pre_change), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = my_pal, dim = dim, na.value="white") +
#   geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class_PD_pre_change), show.legend=F)+
#   bi_scale_color(pal = my_pal, dim = dim, na.value="white") +
#   bi_theme()+ theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm")) + coord_sf(expand = F)
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="pre_change", ylab="SES.PD ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=pre_change, y=SES.PD, style="jenks", dim=dim))+
#   theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_text(angle=45, hjust=1))
# ggdraw() + draw_plot(prechange_pd, 0, 0, 1, 1) + draw_plot(legend, l.x, l.y, l.width, l.height)
# ggsave("figures/choropleth_PD_prechange.png", width=10, height=5.74, units="in", dpi=300, bg="white")
# 
# matchange_pe <- ggplot() +
#   geom_sf(shp2, mapping = aes(fill = bi_class_PE_mat_change), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = my_pal, dim = dim, na.value="white") +
#   geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class_PE_mat_change), show.legend=F)+
#   bi_scale_color(pal = my_pal, dim = dim, na.value="white") +
#   bi_theme()+ theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm")) + coord_sf(expand = F)
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="mat_change", ylab="SES.PE ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=mat_change, y=SES.PE, style="jenks", dim=dim))+
#   theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_text(angle=45, hjust=1))
# ggdraw() + draw_plot(matchange_pe, 0, 0, 1, 1) + draw_plot(legend, l.x, l.y, l.width, l.height)
# ggsave("figures/choropleth_PE_matchange.png", width=10, height=5.74, units="in", dpi=300, bg="white")
# 
# matchange_pd <- ggplot() +
#   geom_sf(shp2, mapping = aes(fill = bi_class_PD_mat_change), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = my_pal, dim = dim, na.value="white") +
#   geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class_PD_mat_change), show.legend=F)+
#   bi_scale_color(pal = my_pal, dim = dim, na.value="white") +
#   bi_theme()+ theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm")) + coord_sf(expand = F)
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="mat_change", ylab="SES.PD ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=mat_change, y=SES.PD, style="jenks", dim=dim))+
#   theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_text(angle=45, hjust=1))
# ggdraw() + draw_plot(matchange_pd, 0, 0, 1, 1) + draw_plot(legend, l.x, l.y, l.width, l.height)
# ggsave("figures/choropleth_PD_matchange.png", width=10, height=5.74, units="in", dpi=300, bg="white")
