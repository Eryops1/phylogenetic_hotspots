# plots and stats for global phylodiversity patterns



wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
gc()
library(sf)
library(ggplot2)
library(cowplot)
library(biscale)
library(spdep)
library(scico)
library(data.table)
library(terra)
if(!dir.exists("figures"))dir.create("figures")
source("99_functions.R")
library(extrafont)
library(broom)
library(scales)
library(lwgeom) # for st_transform_proj
#font_import() # takes 5 minutes if loading all, select wisely
loadfonts()
theme_set(theme_bw()+theme(text=element_text(size=7, family="Helvetica"), 
                           panel.grid=element_blank(), 
                           legend.text.align=1))


# Load data ---------------------------------------------------------------


gallpeters_projection <- "+proj=cea +lon_0=0 +x_0=0 +y_0=0 +lat_ts=45 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
winkel_tripel <- "+proj=wintri +datum=WGS84 +no_defs +over"
# Gall-Peters. Horizontally compressed version of the Lambert equal-area.
# Standard parallels at 45°N/S. Aspect ratio of ~1.6. 

my_projection <- winkel_tripel
shp <- readRDS("data/fin_shp.rds")
shp <- st_as_sf(shp)

names(shp)[grep("SES\\.PD", names(shp))] <- "sesPD"
names(shp)[grep("SES\\.PE", names(shp))] <- "sesPE"


# transform projection
shp <- st_transform(shp, crs=my_projection)

# remove not needed data
shp <- shp[!shp$LEVEL3_COD=="ANT",]
shp <- shp[,-grep("obs_p|obs_rank|reps|LEVEL2|LEVEL1|LEVEL_3_CO|ID|\\.3|_rw|CONTI|REGION|AvTD|TTD|mpd", names(shp))]
names(shp)<- gsub("\\.1", "_mean", names(shp))
names(shp)<- gsub("\\.2", "_sd", names(shp))

# Hotspot shapefile
h <- st_read("data/hotspots_fixed.gpkg")
h <- st_wrap_dateline(h, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
h <- st_transform(h, my_projection)



# Plotting settings
min.area <- 6e+9
thicc_lines <- shp[which(shp$area<min.area),]
shp2 <- shp[!is.na(shp$sesPD),]





# Global patterns  ------------------------------------------------
# # graticule:
# grat_wintri <- 
#   st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9)) %>%
#   st_transform_proj(crs = winkel_tripel)

lcol <- min(thicc_lines$PD_obs)/max(shp$PD_obs)
ucol <- max(thicc_lines$PD_obs)/max(shp$PD_obs)
# new upper limit for second scale:
upl = 1/ucol # upper limit

(pd_map <- ggplot(shp) + 
#    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
    geom_sf(data=shp, aes(fill=PD_obs),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=PD_obs), show.legend=F)+
    # geom_point(data=thicc_lines, aes(color = PD_obs, geometry = geometry),
    #            stat = "sf_coordinates")+
    scale_color_distiller("PD", palette="BuGn", direction=1, values=c(0,upl), trans = "sqrt")+
    scale_fill_distiller("PD", palette="BuGn", direction=1, trans = "sqrt")+
    # scale_colour_scico("PD", palette = "batlow", trans = "sqrt", begin = lcol, end = sqrt(ucol))+
    # scale_fill_scico("PD", palette = "batlow", trans="sqrt")+ #,
    theme_void()+
    theme(legend.position = c(0.2, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align=1,
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10))+
    xlab(" ")+
    coord_sf(datum = NULL) # needed to keep sf from generating graticule (this would fail)
)
lcol <- min(thicc_lines$PE_obs)/max(shp$PE_obs)
ucol <- max(thicc_lines$PE_obs)/max(shp$PE_obs)
# new upper limit for second scale:
upl = 1/ucol # upper limit

(pe_map <- ggplot(shp) + 
    geom_sf(aes(fill=PE_obs),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=PE_obs), show.legend=F)+
    theme_void()+
    scale_color_distiller("PE", palette="BuGn", direction=1, values=c(0,upl), trans = "sqrt")+
    scale_fill_distiller("PE", palette="BuGn", direction=1, trans = "sqrt")+
    # scale_colour_scico("PE", palette = "batlow", trans = "sqrt", begin = lcol, end = sqrt(ucol))+
    # scale_fill_scico("PE", palette = "batlow", trans="sqrt")+ #,
    theme(legend.position = c(0.2, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align=1,
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")+
    coord_sf(datum = NULL)
)
# Simple SR 
lcol <- min(thicc_lines$richness)/max(shp$richness)
ucol <- max(thicc_lines$richness)/max(shp$richness)
# new upper limit for second scale:
upl = 1/ucol # upper limit
(sr_map <- ggplot(shp) + 
    geom_sf(aes(fill=richness),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=richness), show.legend=F)+
    scale_color_distiller("SR", palette="BuGn", direction=1, values=c(0,upl), trans = "sqrt")+
    scale_fill_distiller("SR", palette="BuGn", direction=1, trans = "sqrt")+
    # scale_colour_scico("SR", palette = "batlow", trans = "sqrt", begin = lcol, end = sqrt(ucol))+
    # scale_fill_scico("SR", palette = "batlow", trans="sqrt")+ #,
    theme_void()+
    coord_sf(expand=F, datum=NULL)+
    theme(legend.position = c(0.2, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align=1,
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )
)


# simple WE
thicc_lines <- shp2[which(shp2$area<min.area),]
lcol <- min(thicc_lines$WE)/max(shp$WE)
ucol <- max(thicc_lines$WE)/max(shp$WE)
# new upper limit for second scale:
upl = 1/ucol # upper limit
(we_map <- ggplot(shp2) + 
    geom_sf(aes(fill=WE),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=WE), show.legend=F)+
    scale_color_distiller("WE", palette="BuGn", direction=1, values=c(0,upl), trans = "sqrt")+
    scale_fill_distiller("WE", palette="BuGn", direction=1, trans = "sqrt")+
    # scale_colour_scico("WE", palette = "batlow", trans = "sqrt", begin = lcol, end = sqrt(ucol))+
    # scale_fill_scico("WE", palette = "batlow", trans="sqrt")+ #,
    theme_void()+coord_sf(expand=F, datum=NULL)+
    theme(legend.position = c(0.2, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align=1,
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )+
    xlab(" ")
)


lcol <- 1-min(thicc_lines$sesPD)/(min(shp2$sesPD)-max(shp2$sesPD))
ucol <- max(thicc_lines$sesPD)/max(shp2$sesPD)
# new lower limit for second scale:
lol = -1/(1-lcol) # lower limit
(pd_ses_map <- ggplot(shp2) + 
    geom_sf(aes(fill=sesPD),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=sesPD), show.legend=F)+
    scale_color_distiller(bquote("PD"[std]), palette="BuGn", direction=1, values=c(lol,1))+
    scale_fill_distiller(bquote("PD"[std]), palette="BuGn", direction=1)+
    # scale_colour_scico("PDstd", palette = "batlow", trans = "sqrt", begin = lcol, end = sqrt(ucol))+
    # scale_fill_scico("PDstd", palette = "batlow", trans="sqrt")+ #,
    theme_void()+coord_sf(expand=F, datum=NULL)+
    theme(legend.position = c(0.2, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align=1,
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")
)

#lcol <- min(thicc_lines$sesPE+abs(min(shp2$sesPE)), na.rm=T)/diff(range(shp2$sesPE)) 
lcol <- min(thicc_lines$sesPE, na.rm=T)/diff(range(shp2$sesPE))
ucol <- max(thicc_lines$sesPE, na.rm=T)/diff(range(shp2$sesPE))
# new limits for second scale:
lol = min(shp2$sesPE, na.rm=T) / min(thicc_lines$sesPE, na.rm=T) # lower limit
ul = 1 / ucol # upper limit
# adjusting lower and upper limits at the same time is tricky ^____^

(pe_ses_map <- ggplot(shp2) + 
    geom_sf(aes(fill=sesPE),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=sesPE), show.legend=F)+
    scale_color_distiller(bquote("PE"[std]), palette="BuGn", direction=1, values=c(-lol-1,ul))+
    scale_fill_distiller(bquote("PE"[std]), palette="BuGn", direction=1)+
    # scale_colour_scico("PEstd", palette = "batlow", trans = "sqrt", begin = lcol, end = sqrt(ucol))+
    # scale_fill_scico("PEstd", palette = "batlow", trans="sqrt")+
    theme_void()+coord_sf(expand=F, datum=NULL)+
    theme(legend.position = c(0.2, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.background = element_blank(),
          legend.text.align = 1,
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")
)

# PD, PE
png("figures/maps_unstandardised.png", width=12.5, height=3.5, units = "in", res = 600, bg = "white")
plot_grid(pd_map, pe_map, ncol = 2,
          labels=c("A","B"), label_fontface=1)
#ggsave("figures/maps_unstandardised.png", width=12.5, height=3.5, units = "in", dpi = 600, bg = "white")
dev.off()

# SR, WE, PDstd, PEstd
png("figures/maps.png", width=12.5, height=7.5, units = "in", res = 600, bg = "white")
plot_grid(sr_map+ggtitle("Species richness\n")+theme(plot.title = element_text(hjust = 0.5, size=8)), 
          we_map+ggtitle("Weighted endemism\n")+theme(plot.title = element_text(hjust = 0.5, size=8)), 
          pd_ses_map+ggtitle("Phylogenetic diversity, standardized effect size\n")+theme(plot.title = element_text(hjust = 0.5, size=8)), 
          pe_ses_map+ggtitle("Phylogenetic endemism, standardized effect size\n")+theme(plot.title = element_text(hjust = 0.5, size=8)),
          ncol = 2, labels=c("A","B","C","D"), label_fontface=1, label_fontfamily="Helvetica", 
          scale=1)
#ggsave("figures/maps.png", width=12.5, height=7.5, units = "in", dpi = 300, bg = "white")
dev.off()

# ALL
png("figures/maps_all.png", width=12.5, height=10.3, units = "in", res = 600, bg = "white")
plot_grid(sr_map+ggtitle("Species richness\n")+theme(plot.title = element_text(hjust = 0.5, size=8)), 
          we_map+ggtitle("Weighted endemism\n")+theme(plot.title = element_text(hjust = 0.5, size=8)),
          pd_map+ggtitle("Phylogenetic diversity\n")+theme(plot.title = element_text(hjust = 0.5, size=8)),
          pe_map+ggtitle("Phylogenetic endemism\n")+theme(plot.title = element_text(hjust = 0.5, size=8)),
          pd_ses_map+ggtitle("Phylogenetic diversity, standardized effect size\n")+theme(plot.title = element_text(hjust = 0.5, size=8)), 
          pe_ses_map+ggtitle("Phylogenetic endemism, standardized effect size\n")+theme(plot.title = element_text(hjust = 0.5, size=8)),
          ncol = 2, labels=c("A","B","C","D","E","F"), label_fontface=1, label_fontfamily="Helvetica", 
          scale=1)
dev.off()



# Conservation hotspot coverage 
ggplot(shp2) + 
  geom_sf(aes(col=hotspot_coverage), show.legend=T)+
  geom_sf(aes(fill=hotspot_coverage), col="grey80", lwd=.1)+
  # scale_fill_gradient("hotspot \ncoverage", low="white", high="red")+
  # scale_color_gradient(low="white", high="red")+
  scale_color_distiller("hotspot \ncoverage", palette="BuGn", direction=1)+
  scale_fill_distiller("hotspot \ncoverage", palette="BuGn", direction=1)+
  
  theme_void()+
  theme(legend.position = c(0.2, 0.3),
        legend.key.height = unit(6,"mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title=element_text(size=9))+
  coord_sf(expand=F)
ggsave("figures/hotspot_coverage.png", units="in", dpi=300, width=7, height=4.1)







# Spatial correlations --------------------------------------------------

# sesPD + sesPE VS taxonomic measures
# reduce to relevant variables:
dat <- shp[,c("LEVEL3_COD", "richness", "WE", "sesPD", "sesPE")]
names(dat)
dat <- na.omit(dat)

nb <- spdep::poly2nb(dat, row.names = dat$LEVEL3_COD)
col.W <- nb2listw(nb, style="W", zero.policy = TRUE)
lee.dat <- st_drop_geometry(dat[,-which(names(dat)=="LEVEL3_COD")])

les <- apply(lee.dat[,!names(lee.dat)%in%"sesPD"], 2, lee.test, y=lee.dat$sesPD, 
             listw=col.W, zero.policy = TRUE, alternative="two.sided")


cor.test(lee.dat$sesPD, lee.dat$sesPE, method="p")
lee.test(lee.dat$sesPD, lee.dat$sesPE, listw=col.W, zero.policy = TRUE, alternative="two.sided")
lee.mc(lee.dat$sesPD, lee.dat$sesPE, listw=col.W, zero.policy = TRUE, alternative="greater", nsim=999)

cor.test(lee.dat$sesPD, lee.dat$richness, method="p")
lee.test(lee.dat$sesPD, lee.dat$richness, listw=col.W, zero.policy = TRUE, alternative="two.sided")
lee.mc(lee.dat$sesPD, lee.dat$richness, listw=col.W, zero.policy = TRUE, alternative="less", nsim=999)

cor.test(lee.dat$sesPE, lee.dat$WE, method="p")
lee.test(lee.dat$sesPE, lee.dat$richness, listw=col.W, zero.policy = TRUE, alternative="two.sided")
lee.mc(lee.dat$sesPE, lee.dat$richness, listw=col.W, zero.policy = TRUE, alternative="greater", nsim=999)






# Scaling effects area + SR -----------------------------------------------

## PD ----
stand_fun <- function(x){(x-mean(x))/sd(x)}

tmp <- st_drop_geometry(shp2[,c("area", "PD_obs", "richness", "sesPD", "LEVEL_NAME")])
tmp$PD_richness <- tmp$PD/tmp$richness
names(tmp) <- c("area", "PD", "SR", "sesPD", "LEVEL_NAME", "PD_richness")

tmp[,c(2,3)] <- apply(tmp[,c(2,3)], 2, normalized)
setDT(tmp)
# reshape
tmp <- data.table::melt(tmp, id.var=c("area", "LEVEL_NAME"))

# plot
bc <- c("#52548D", "#C57391", "#EFB984") # PD, PD+PE, PE
(p1 <- ggplot(tmp[grep("sesPD", tmp$variable),], aes(x=area, y=value, col=variable))+
  geom_point(alpha=0.1)+
  geom_smooth(se=F, method="lm", lwd=0.7)+
  scale_x_continuous(trans="log10", #breaks = trans_breaks("log10", function(x) 10^x),
                     labels = math_format(
                      format = function(x){number(log10(x), accuracy = 1)}))+
  labs(y=expression('PD'[std]))+
  scale_color_manual(values=c(bc[1], "grey20"))+
  theme(legend.position="none"))

(p2 <- ggplot(shp2, aes(x=richness, y=sesPD))+
    geom_point(alpha=0.1)+
    geom_smooth(se=F, method="lm", lwd=0.7, col="grey20")+
    scale_x_continuous(trans="log10", 
                       labels = math_format(
                         format = function(x){number(log10(x), accuracy = 1)}))+
    scale_color_manual(values=c("grey20"))+
    labs(y=expression('PD'[std]),
         x="species richness")+
    theme(legend.position="none"))

(pempty <- ggplot(tmp[grep("PD_richness", tmp$variable),], aes(x=area, y=value, col=variable))+
    geom_point(alpha=0)+ theme_void()+ theme(legend.position="none"))

(p3 <- ggplot(tmp[grep("^PD$|^SR", tmp$variable),], aes(x=area, y=value, col=variable))+
  geom_point(alpha=0.1)+
  scale_x_continuous(trans="log10", #breaks = trans_breaks("log10", function(x) 10^x),
                     labels = math_format(
                       format = function(x){number(log10(x), accuracy = 1)}))+
  ylab("standardised value")+
  scale_color_manual(values=c(bc[2], "grey20"))+
  geom_smooth(se=F, method="lm", alpha=0.5, lwd=.7)+
    scale_y_log10()+
  theme(legend.position=c(0.12, 0.8),
        legend.background=element_blank(), 
        legend.spacing.y=unit(1,"mm"),
        legend.title=element_blank(),
        legend.key.width=unit(3, "mm"),
        strip.background=element_blank()))

plot_grid(plot_grid(p3, pempty, labels=c("A","B"),label_y=1,  label_size=11, label_fontface="plain"),
          plot_grid(p1, p2, labels=c("C","D"),label_y=1.02,  label_size=11, label_fontface="plain"), 
          ncol=1, label_size=11, label_fontface="plain")

ggsave("figures/scaling_log2.pdf", height=3.4, width=3.4, unit="in", dpi=300)

## PD regression coefficients ----
shp2$richness_norm_log <- log(normalized(shp2$richness))
shp2$PD_obs_norm_log <- log(normalized(shp2$PD_obs))

tmp <- st_drop_geometry(shp2[,c("richness_norm_log", "area")])
tmp <- tmp[!is.infinite(rowSums(tmp)),]
m_SR <- lm(data=tmp, richness_norm_log~log(area))

tmp <- st_drop_geometry(shp2[,c("PD_obs_norm_log", "area")])
tmp <- tmp[!is.infinite(rowSums(tmp)),]
m_PD <- lm(data=tmp, PD_obs_norm_log~log(area))

mlist <- list(m_SR, m_PD)
names(mlist) <- c("SR ~ area", "PD ~ area")
(tab <- huxtable::huxreg(mlist,
                         error_format = "({std.error})",
                         error_pos = "below",
                         number_format = "%.2f",
                         align = ".",
                         stars = c(`***` = 0.001),
                         statistics = c(R2 = "r.squared")))
huxtable::font_size(tab) <- 13
huxtable::quick_html(file="PD_regression_table", tab, open=F) 


## PE ----
tmp <- st_drop_geometry(shp2[,c("area", "PE_obs", "richness", "sesPE", "LEVEL_NAME")])
names(tmp) <- c("area", "PE", "SR", "PEstd", "LEVEL_NAME")

tmp[,c(2,3)] <- apply(tmp[,c(2,3)], 2, normalized)
setDT(tmp)
# reshape
tmp <- data.table::melt(tmp, id.var=c("area", "LEVEL_NAME"))
bc <- c("#52548D", "#C57391", "#EFB984") # PD, PD+PE, PE
(p12 <- ggplot(tmp[grep("PEstd", tmp$variable),], aes(x=area, y=value, col=variable))+
    geom_point(alpha=0.1)+
    geom_smooth(se=F, method="lm", lwd=0.7)+
    scale_x_continuous(trans="log", 
                       labels = math_format(
                         format = function(x){number(log10(x), accuracy = 1)}))+
    ylab(expression('PE'[std]))+
    scale_color_manual(values=c(bc[1], "grey20"))+
    theme(legend.position="none"))

(p22 <- ggplot(shp2, aes(x=richness, y=sesPE))+
    geom_point(alpha=0.1)+
    geom_smooth(se=F, method="lm", lwd=0.7, col="grey20")+
    scale_x_continuous("species richness", trans="log10", labels = math_format(
      format = function(x){number(log10(x), accuracy = 1)}))+
    scale_color_manual(values=c("grey20"))+
    labs(y=expression('PE'[std]),
         x="species richness")+
    theme(legend.position="none"))

(pempty <- ggplot(tmp[grep("PD_richness", tmp$variable),], aes(x=area, y=value, col=variable))+
    geom_point(alpha=0)+ theme_void()+ theme(legend.position="none"))

(p32 <- ggplot(tmp[grep("^PE$|^SR", tmp$variable),], aes(x=area, y=value, col=variable))+
    geom_point(alpha=0.1)+
    scale_x_continuous(trans="log", 
                       labels = math_format(
                         format = function(x){number(log10(x), accuracy = 1)}))+
    ylab("standardised value")+
    scale_color_manual(values=c(bc[2], "grey20"))+
    geom_smooth(se=F, method="lm", alpha=0.5, lwd=.7)+
    scale_y_log10()+
    theme(legend.position=c(0.12, 0.8),
          legend.background=element_blank(), 
          legend.spacing.y=unit(1,"mm"),
          legend.title=element_blank(),
          legend.key.width=unit(3, "mm"),
          strip.background=element_blank()))

plot_grid(plot_grid(p32, pempty, labels=c("A","B"),label_y=1,  label_size=11, label_fontface="plain"),
          plot_grid(p12, p22, labels=c("C","D"),label_y=1.02,  label_size=11, label_fontface="plain"), 
          ncol=1, label_size=11, label_fontface="plain")

ggsave("figures/scaling_PE.png", height=3.4, width=3.4, unit="in", dpi=300)



##  PE regression coefficients ------
shp2$richness_norm_log <- log(normalized(shp2$richness))
shp2$PE_obs_norm_log <- log(normalized(shp2$PE_obs))

tmp <- st_drop_geometry(shp2[,c("richness_norm_log", "area")])
tmp <- tmp[!is.infinite(rowSums(tmp)),]
m_SR <- lm(data=tmp, richness_norm_log~log(area))

tmp <- st_drop_geometry(shp2[,c("PE_obs_norm_log", "area")])
tmp <- tmp[!is.infinite(rowSums(tmp)),]
m_PE <- lm(data=tmp, PE_obs_norm_log~log(area))

mlist <- list(m_SR, m_PE)
names(mlist) <- c("SR ~ area", "PE ~ area")
(tab <- huxtable::huxreg(mlist,
                         error_format = "({std.error})",
                         error_pos = "below",
                         number_format = "%.2f",
                         align = ".",
                         stars = c(`***` = 0.001),
                         statistics = c(R2 = "r.squared")))
huxtable::font_size(tab) <- 13
huxtable::quick_html(file = "PE_model.html", tab, open = F)


save(list=c("shp", "shp2"), file="data/workspace_point0.RData")





## plot sesPE vs sesPD ----
ggplot(shp2, aes(x=sesPD, y=sesPE, label=LEVEL_NAME))+
  geom_point()+
  geom_label(nudge_x=1, hjust=0, label.size=0, label.padding=unit(0.5,"mm"), 
             alpha=0, color = alpha('black', .5), size=3)+
  scale_x_continuous(limits=c(-80,20))
#ggsave("figures/SES_PE_vs_SES_PD.png", width=5, height=5)



