
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
if(!dir.exists("figures"))dir.create("figures")
source("99_functions.R")



load(file="data/workspace_point1.RData")
min.area <- 6e+9
gallpeters_projection <- "+proj=cea +lon_0=0 +x_0=0 +y_0=0 +lat_ts=45 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
my_projection <- gallpeters_projection
shp <- readRDS("data/fin_shp.rds")

# Conservation hotspots overlap ----

tmp <- shp2[shp2$PD_hotspot=="3-3" | shp2$PE_hotspot %in% c("4-3", "3-4", "4-4"), ]
s <- as(st_geometry(tmp), "Spatial")

h <- st_read("data/hotspots_fixed.gpkg")
h <- st_wrap_dateline(h, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
h <- st_transform(h, my_projection)
m3 <- as(st_geometry(h), "Spatial")

s2 <- vect(s)
m32 <- vect(m3)

res <- list() # list of length hotspots
for(i in 1:nrow(s2)){
  s3 <- s2[i,]
  res[[i]] <- relate(m32, s3, "intersects")
}
res <- as.data.frame(res)
names(res) <- tmp$LEVEL_NAME
res$sum <- rowSums(res)
res$conservation_hotspots <- h$NAME

res$conservation_hotspots[res$sum>0] # covered by at least one phylohotspot
colSums(res[,-ncol(res)]) # all phylohotspots overlap with conservation hotspots





# Hotspot map  ------------------------------------------

ha_united <- st_read("data/hotspots_merged_cleaned.gpkg")
# loads the conservation hotspots that have been processed in QGIS to merge
# inner and outher boundaries.

# change projection + solve dateline issue
ha_united <- st_wrap_dateline(ha_united, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
ha_united <- st_transform(ha_united, my_projection)


## Hotspot map ----------------------------------------

diss <- st_read("data/level3_dissolved.gpkg")
diss <- st_simplify(diss)

thicc_lines <- shp2[which(shp2$area<6e+08),]

bc <- c("#52548D", "#C57391", "#EFB984")
shp3 <- st_simplify(shp2)
thicc_lines3 <- st_simplify(thicc_lines)
ha_united3 <- st_simplify(ha_united)

shp2$hotspot_type[shp2$hotspot_type=="no hotspot"] = NA

ggplot() +
  geom_sf(shp3, mapping = aes(fill=hotspot_type), color=NA, size=0.1) +
  geom_sf(data=ha_united3, col=alpha("white",0.9), fill="grey", alpha=.3, size=.5, lty="dotted")+
  geom_sf(data=diss, lwd=0.1, col="grey40", fill=NA, show.legend=F)+
  scale_fill_manual("Hotspot type", values=bc, na.value=NA, na.translate=F, 
                    breaks=c("PD", "PD & PE", "PE"),
                    labels=c(expression('PD'[std]), bquote(PD[std]+PE[std]), expression('PE'[std])))+
 theme_void()+ coord_sf(ylim=c(-7001217, 8928740))+
  theme(legend.position = c(0.15, 0.15),
        legend.key.height = unit(5,"mm"),
        legend.background = element_rect(fill="white", colour="white", size=3),
        #        legend.key = element_blank(),
        panel.background = element_blank(), 
        legend.text=element_text(size=7),
        legend.title=element_text(size=8))
ggsave("figures/one_hotspots_map__no_text.pdf", width=160, height=90, units="mm", dpi=300, bg="white")  



### alt panel layout ----
labsize <- 2

p1 <- ggplot() +
  geom_sf(shp3, mapping = aes(fill=hotspot_type), color=NA, size=0.1) +
  geom_sf(data=thicc_lines3, lwd=1, col=NA, show.legend=F)+
  geom_sf(data=diss, lwd=0.1, col="black", fill=NA, show.legend=F)+
  scale_fill_manual("Hotspot\ntype", values=bc, na.value="white")+
  geom_sf(data=ha_united3, col=NA, fill="grey80", alpha=0.5, lwd=0)+
  geom_sf_text(data=ha_united3, aes(label=NAME), size=labsize)+
  theme_void()+ coord_sf(xlim=c(-8292430, -5002430), ylim=c(2031217, 3728740))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        plot.margin=margin(1,1,1,1, "pt"),
        plot.background = element_rect(colour = "black", fill=NA, size=.5))
p2 <- ggplot() +
  geom_sf(shp3, mapping = aes(fill=hotspot_type), color=NA, size=0.1) +
  geom_sf(data=thicc_lines3, lwd=1, col=NA, show.legend=F)+
  geom_sf(data=diss, lwd=0.1, col="black", fill=NA, show.legend=F)+
  scale_fill_manual("Hotspot\ntype", values=bc, na.value="white")+
  geom_sf(data=ha_united3, col=NA, fill="grey80", alpha=0.5, lwd=0)+
  geom_sf_text(data=ha_united, aes(label=NAME), size=labsize)+
  theme_void()+ coord_sf(xlim=c(802430, 3602430), ylim=c(-5231217, -3728740))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        plot.margin=margin(1,1,1,1, "pt"),
        plot.background = element_rect(colour = "black", fill=NA, size=.5))
p3 <- ggplot() +
  geom_sf(shp3, mapping = aes(fill=hotspot_type), color=NA, size=0.1) +
  geom_sf(data=thicc_lines3, lwd=1, col=NA, show.legend=F)+
  geom_sf(data=diss, lwd=0.1, col="black", fill=NA, show.legend=F)+
  scale_fill_manual("Hotspot\ntype", values=bc, na.value="white")+
  geom_sf(data=ha_united3, col=NA, fill="grey80", alpha=0.5, lwd=0)+
  geom_sf_text(data=ha_united3, aes(label=NAME), size=labsize)+
  theme_void()+ coord_sf(xlim=c(8002430, 13502430), ylim=c(-5531217, -1528740))+
  theme(legend.position = "none", 
        plot.margin=margin(1,1,1,1, "pt"),
        plot.background = element_rect(colour = "black", fill=NA, size=.5))
p4 <- ggplot() +
  geom_sf(shp3, mapping = aes(fill=hotspot_type), color=NA, size=0.1) +
  geom_sf(data=thicc_lines3, lwd=1, col=NA, show.legend=F)+
  geom_sf(data=diss, lwd=0.1, col="black", fill=NA, show.legend=F)+
  scale_fill_manual("Hotspot\ntype", values=bc, na.value="white")+
  geom_sf(data=ha_united3, col=NA, fill="grey80", alpha=0.5, lwd=0)+
  geom_sf_text(data=ha_united3, aes(label=NAME), size=labsize)+
  theme_void()+ coord_sf(xlim=c(4802430, 10502430), ylim=c(-1328740, 5228740))+
  theme(legend.position = c(0.1, 0.05),
        legend.key.height = unit(5,"mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(), 
        plot.margin=margin(1,1,1,1, "pt"),
        plot.background = element_rect(colour = "black", fill=NA, size=.5))

plot_grid(plot_grid(p1,p2,p3, ncol=1, scale = 0.96, labels=c("A","B", "C"), label_size=10, label_fontface="plain"),
          p4, ncol=2, scale = 0.96, labels=c("", "D"), label_size=11, label_fontface="plain", 
          rel_widths=c(1, 2))

ggsave("figures/fig3_panel_layout.pdf", width=150, height=118, units="mm", dpi=300, bg="white")



