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


# Load data ---------------------------------------------------------------

shp <- readRDS("data/fin_shp.rds")
shp <- st_as_sf(shp)


# transform projection
my_projection <- "+proj=wintri +datum=WGS84 +no_defs +over"
shp <- st_transform(shp, crs=my_projection)



# remove not needed data + rename
keep = c("LEVEL3_COD", "LEVEL_3_CO", "LEVEL_NAME", "area", "richness", "PD_obs",
         "PE_obs", "WE", "SE", "PDE", "mrd", "sub_trop_mbf", "sub_trop_dbf", "sub_trop_cf",
         "temp_bmf", "temp_cf", "boreal_f_taiga", "sub_trop_gss", "temp_gss", 
         "flooded_gs", "mont_gs", "tundra", "medit_fws", "deserts_x_shrub", 
         "hfp.1", "hfp.2", "hfp.3", "deforestation2", "hotspot_coverage", 
         "mat_change", "pre_change", "sr_complementarity", "sr_added", 
         "pd_complementarity", "pd_added", "geometry", "pde_complementarity", 
         "sr_added_perc", "pd_added_perc", "pde_added_perc", "sr_comp_top10", 
         "sr_comp_50", "pd_comp_50", "pd_comp_top10", "pde_comp_50", "topsr", "toppd")
shp <- shp[,grep(paste(keep, collapse="|"), names(shp))]
names(shp)[grep("\\.1", names(shp))] = "hfp_mean"


# Conservation hotspot shapefile
h <- st_read("data/hotspots_fixed.gpkg")
h <- st_wrap_dateline(h, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
h <- st_transform(h, my_projection)



# Min area highlights and shp2 (no NA)
min.area <- 6e+9
min.area <- units::as_units(min.area, "m2")
thicc_lines <- shp[which(shp$area<min.area),]






# Global patterns  ------------------------------------------------
library(dplyr)

# graticule
grat_wintri <-
  st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9)) %>%
  st_transform_proj(crs = my_projection) #%>%
#  st_convex_hull() %>%
#  summarise(geometry = st_union(geometry))


grat_wintri <-
  st_graticule(lat = c(-98.9, 89.9),
               lon = c(-179.9, 179.9)) %>%
  st_transform_proj(crs = my_projection) #%>%
#  st_convex_hull() %>%
#  summarise(geometry = st_union(geometry))


# colors
bc <- c("#35abc4", "#4b9e31", "#eeea40")
bc <- c("#0FC4EA", "#4b9e31", "#eeea40")

# set map theme
theme_set(theme_void()+
            theme(legend.position = c(0.1, 0.3),
                  legend.key.height = unit(6,"mm"),
                  legend.key.width = unit(4,"mm"),
                  legend.background = element_blank(),
                  legend.key = element_blank(),
                  legend.text.align=1,
                  panel.background = element_blank(),
                  panel.border = element_blank(),
                  text = element_text(size = 10)))




# FIG 1 ------------------------------------------------------------------

# PD

# upper color for islands
n2 = max(thicc_lines$PD_obs)/max(shp$PD_obs)
x <- seq(0, 1, length.out=diff(range(shp$PD_obs)))
n3 = which.min(abs(x-n2))
cols = seq_gradient_pal("white", bc[1])(x)
maxcol = cols[n3]

(pd_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
    geom_sf(data=shp, aes(fill=PD_obs),lwd=0.25/.pt, col="grey50") + 
    scale_fill_gradient("PD", low="gray95", high=bc[1])+
    geom_sf(data=thicc_lines, lwd=1, aes(col=as.numeric(PD_obs)), show.legend=F)+
    scale_color_gradient(low="gray95", high=maxcol)+                   # islands
    coord_sf(expand=F, datum=NULL) # needed to keep sf from generating graticule (this would fail)
)

# PE

n2 = max(thicc_lines$PE_obs)/max(shp$PE_obs)
x <- seq(0, 1, length.out=diff(range(shp$PE_obs)))
n3 = which.min(abs(x-n2))
cols = seq_gradient_pal("white", "grey50")(x)
maxcol = cols[n3]

(pe_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
    geom_sf(data=shp, aes(fill=PE_obs),lwd=0.25/.pt, col="gray50") + 
    scale_fill_gradient("PE", low="gray95", high="grey20", trans="sqrt")+
    geom_sf(data=thicc_lines, lwd=1, aes(col=as.numeric(PE_obs)), show.legend=F)+
    scale_color_gradient(low="gray95", high=maxcol)+
    coord_sf(expand=F, datum=NULL) # needed to keep sf from generating graticule (this would fail)
)

# SR 

n2 = max(thicc_lines$richness)/max(shp$richness)
n3 = which.min(abs(x-n2))
x <- seq(0, 1, length.out=diff(range(shp$richness)))
cols = seq_gradient_pal("white", "grey50")(x)
maxcol = cols[n3]

(sr_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
    geom_sf(aes(fill=richness),lwd=0.25/.pt, col="gray50") + 
    scale_fill_gradient("SR", low="gray95", high="grey20")+
    geom_sf(data=thicc_lines, lwd=1, aes(col=as.numeric(richness)), show.legend=F)+
    scale_color_gradient(low="gray95", high=maxcol)+
    coord_sf(expand=F, datum=NULL)
)

# WE

n2 = max(thicc_lines$WE)/max(shp$WE)
n3 = which.min(abs(x-n2))
x <- seq(0, 1, length.out=diff(range(shp$WE)))
cols = seq_gradient_pal("white", "grey50")(x)
maxcol = cols[n3]

(we_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
    geom_sf(aes(fill=WE),lwd=0.25/.pt, col="gray50") + 
    scale_fill_gradient("WE", low="gray95", high="grey20", trans="sqrt")+
    geom_sf(data=thicc_lines, lwd=1, aes(col=as.numeric(WE)), show.legend=F)+
    scale_color_gradient(low="gray95", high=maxcol)+
    coord_sf(expand=F, datum=NULL)
)

# SE

n2 = max(thicc_lines$SE)/max(shp$SE)
n3 = which.min(abs(x-n2))
x <- seq(0, 1, length.out=diff(range(shp$SE)))
cols = seq_gradient_pal("white", "grey50")(x)
maxcol = cols[n3]

(se_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
    geom_sf(aes(fill=SE),lwd=0.25/.pt, col="gray50") + 
    scale_fill_gradient("SE", low="gray95", high="#8E6496", trans="sqrt")+
    geom_sf(data=thicc_lines, lwd=1, aes(col=as.numeric(WE)), show.legend=F)+
    scale_color_gradient(low="gray95", high=maxcol)+
    coord_sf(expand=F, datum=NULL)
)

# PD endemism

n2 = max(thicc_lines$PDE)/max(shp$PDE)
n3 = which.min(abs(x-n2))
x <- seq(0, 1, length.out=diff(range(shp$PDE)))
cols = seq_gradient_pal("white", bc[2])(x)
maxcol = cols[n3]

(pde_map <- ggplot(shp) + 
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
    geom_sf(data=shp, aes(fill=PDE),lwd=0.25/.pt, col="grey50") + 
    scale_fill_gradient("PD\nendemism", low="gray95", high=bc[2], trans="sqrt")+
    geom_sf(data=thicc_lines, lwd=1, aes(col=as.numeric(PDE)), show.legend=F)+
    scale_color_gradient(low="gray95", high=maxcol)+
    coord_sf(expand=F, datum=NULL)
)


# 
# plot_grid(sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5)), 
#           #we_map+ggtitle("Weighted endemism\n")+theme(plot.title = element_text(hjust = 0.5)), 
#           pd_map+ggtitle("Phylogenetic diversity")+theme(plot.title = element_text(hjust = 0.5)), 
#           pde_map+ggtitle("PD endemism")+theme(plot.title = element_text(hjust = 0.5)),
#           ncol = 1, labels=c("A","B","C"), label_fontface=1, label_fontfamily="Helvetica", 
#           scale=1)
# ggsave("figures/fig1.png", width=5, height=10, units = "in", dpi = 300, bg = "white")

plot_grid(sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5)), 
          pd_map+ggtitle("Phylogenetic diversity")+theme(plot.title = element_text(hjust = 0.5)),
          se_map+ggtitle("Species endemism")+theme(plot.title = element_text(hjust = 0.5)), 
          pde_map+ggtitle("PD endemism")+theme(plot.title = element_text(hjust = 0.5)),
          ncol = 2, labels=c("(a)","(b)","(c)","(d)"), label_fontface=2, label_fontfamily="Helvetica", 
          scale=1)
ggsave("figures/fig1.png", width=10, height=6.5, units = "in", dpi = 300, bg = "white")


# FIG 1, area-corrected  ------------------------------------------------------------------

# standardize for area
shp_area = shp[,c("LEVEL3_COD" ,"richness", "PD_obs", "PDE", "SE", "area")]
shp_area$richness = as.numeric(shp_area$richness/shp_area$area)
shp_area$PD_obs = as.numeric(shp_area$PD_obs/shp_area$area)
shp_area$PDE = as.numeric(shp_area$PDE/shp_area$area)
shp_area$SE = as.numeric(shp_area$SE/shp_area$area)

# PD
(pd_map <- ggplot(shp_area) + 
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
    geom_sf(data=shp_area, aes(fill=PD_obs, col=PD_obs),lwd=3/.pt) + 
    scale_fill_gradient("PD", low="gray95", high=bc[1], trans="log")+
    scale_color_gradient("PD", low="gray95", high=bc[1], trans="log")+
    geom_sf(data=shp_area, aes(fill=PD_obs), lwd=0.25/.pt, col="grey50") + 
    coord_sf(expand=F, datum=NULL)
)

# SR 
(sr_map <- ggplot(shp_area) + 
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
    geom_sf(aes(fill=richness, col=richness),lwd=3/.pt) + 
    scale_fill_gradient("SR", low="gray95", high="grey20", trans="log")+
    scale_color_gradient("SR", low="gray95", high="grey20", trans="log")+
    geom_sf(aes(fill=richness),lwd=0.25/.pt, col="gray50") + 
    coord_sf(expand=F, datum=NULL)
)

# SE
(se_map <- ggplot(shp_area) + 
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
    geom_sf(aes(fill=SE, col=SE),lwd=3/.pt) + 
    scale_fill_gradient("Species\nendemism", low="gray95", high="#8E6496", trans="log")+
    scale_color_gradient("Species\nendemism", low="gray95", high="#8E6496", trans="log")+
    geom_sf(aes(fill=SE),lwd=0.25/.pt, col="gray50") + 
    coord_sf(expand=F, datum=NULL)
)

# PD endemism
(pde_map <- ggplot(shp_area) + 
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
    geom_sf(aes(fill=PDE, col=PDE),lwd=3/.pt) + 
    scale_fill_gradient("PD\nendemism", low="gray95", high=bc[2], trans="log")+
    scale_color_gradient("PD\nendemism", low="gray95", high=bc[2], trans="log")+
    geom_sf(aes(fill=PDE),lwd=0.25/.pt, col="gray50") + 
    coord_sf(expand=F, datum=NULL)
)

png("figures/fig1_area.png", width=10, height=6.5, units = "in", res = 300, bg = "white")
plot_grid(sr_map+ggtitle("Species richness/area")+theme(plot.title = element_text(hjust = 0.5)), 
          pd_map+ggtitle("Phylogenetic diversity/area")+theme(plot.title = element_text(hjust = 0.5)),
          se_map+ggtitle("Species endemism/area")+theme(plot.title = element_text(hjust = 0.5)), 
          pde_map+ggtitle("PD endemism/area")+theme(plot.title = element_text(hjust = 0.5)),
          ncol = 2, labels=c("(a)","(b)","(c)","(d)"), label_fontface=2, label_fontfamily="Helvetica", 
          scale=1)
dev.off()
#ggsave("figures/fig1_extended_area.png", width=10, height=6.5, units = "in", dpi = 300, bg = "white")






# ALT color schemes

# library(scico)
# library(MetBrewer)
# MetBrewer::colorblind_palettes
# mets = c("Archambault", "Cassatt1","Cassatt2" ,"Demuth","Derain","Egypt","Greek",
#          "Hiroshige","Hokusai2" ,"Hokusai3","Ingres","Isfahan1","Isfahan2","Java",
#          "Johnson","Kandinsky","Morgenstern","OKeeffe1","OKeeffe2","Pillement",
#          "Tam","Troy", "VanGogh3", "Veronese")

plot_grid(sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
            scale_fill_gradient("Species richness", low="gray95", high="grey15"),
          pd_map+ggtitle("Phylogenetic diversity")+theme(plot.title = element_text(hjust = 0.5))+
            scale_fill_gradient("PD", low="gray95", high="#184C57"),
          pde_map+ggtitle("PD endemism")+theme(plot.title = element_text(hjust = 0.5))+
            scale_fill_gradient("PD endemism", trans="sqrt", low="gray95", high="#3C7E27"),
          ncol = 1, labels=c("A","B","C"), label_fontface=1, label_fontfamily="Helvetica", 
          scale=1)
ggsave("figures/fig1_outlines.png", width=5, height=10, units = "in", dpi = 300, bg = "white")

# plot_grid(sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             #scale_fill_viridis_c("Species richness", option="plasma")
#             scale_fill_gradientn(colors = rev(met.brewer("Cassatt1"))), 
#           pd_map+ggtitle("Phylogenetic diversity")+theme(plot.title = element_text(hjust = 0.5))+
#             #scale_fill_scico("PD", palette="batlow")
#             scale_fill_gradientn(colors = rev(met.brewer(mets[15]))), 
#           pde_map+ggtitle("PD endemism")+theme(plot.title = element_text(hjust = 0.5))+
#             #scale_fill_viridis_c("PD endemism", trans="sqrt"),
#             scale_fill_gradientn(colors = met.brewer("Cassatt1")),
#           ncol = 1, labels=c("A","B","C"), label_fontface=1, label_fontfamily="Helvetica", 
#           scale=1)

# plot_grid(sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[1])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[2])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[3])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[4])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[5])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[6])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[7])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[8])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[9])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[10])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[11])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[12])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[13])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[14])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[15])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[16])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[17])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[18])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[19])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[20])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[21])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[22])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[23])),
#           sr_map+ggtitle("Species richness")+theme(plot.title = element_text(hjust = 0.5))+
#             scale_fill_gradientn(colors = met.brewer(mets[24])),
#           
#           ncol=6)

# Spatial correlations --------------------------------------------------


dat <- shp[,c("LEVEL_NAME", "richness", "PD_obs", "PDE", 
                               "WE", "PE_obs")]
nb <- spdep::poly2nb(dat, row.names = dat$LEVEL3_COD)
col.W <- nb2listw(nb, style="W", zero.policy = TRUE)
lee.dat <- st_drop_geometry(dat[,-which(names(dat)=="LEVEL_NAME")])

cor.test(lee.dat$richness, lee.dat$PD_obs, method="s")
lee.test(lee.dat$richness, lee.dat$PD_obs, listw=col.W, zero.policy = TRUE, alternative="two.sided")
lee.mc(lee.dat$richness, lee.dat$PD_obs, listw=col.W, zero.policy = TRUE, alternative="greater", nsim=999)

cor.test(lee.dat$richness, lee.dat$PDE, method="s")
lee.test(lee.dat$richness, lee.dat$PDE, listw=col.W, zero.policy = TRUE, alternative="two.sided")
lee.mc(lee.dat$richness, lee.dat$PDE, listw=col.W, zero.policy = TRUE, alternative="greater", nsim=999)

cor.test(lee.dat$PDE, lee.dat$PE_obs, method="s") # 0.87
plot(sqrt(lee.dat$PDE), lee.dat$PE_obs)
lee.test(sqrt(lee.dat$PDE), lee.dat$PE_obs, listw=col.W, zero.policy = TRUE, alternative="two.sided")
lee.test(lee.dat$PDE, lee.dat$PE_obs, listw=col.W, zero.policy = TRUE, alternative="two.sided")
lee.mc(lee.dat$PDE, lee.dat$PE_obs, listw=col.W, zero.policy = TRUE, alternative="less", nsim=999)

cor.test(lee.dat$PDE, lee.dat$WE, method="s")
lee.test(lee.dat$PDE, lee.dat$WE, listw=col.W, zero.policy = TRUE, alternative="two.sided")
lee.mc(lee.dat$PDE, lee.dat$WE, listw=col.W, zero.policy = TRUE, alternative="greater", nsim=999)





View(st_drop_geometry(dat))





# FIG 5 Threat spectrum -------------------------------------------------------------

# hotspots are PD compt top 10
plot.df <- st_drop_geometry(shp[,c("LEVEL_NAME", "deforestation2", "hfp_mean", 
                                    "pre_change", "mat_change", "pd_comp_top10")])
plot.df <- reshape::melt(plot.df, id.vars=c("LEVEL_NAME", "pd_comp_top10"))
ypos <- rep(rev(seq(0.2,0.80,0.04)), 4)
#plot.df$value[which(plot.df$value=="NaN")] = NA

plot.df <- plot.df[order(plot.df$variable, plot.df$value),]
plot.df$x <- rep(1:(nrow(plot.df)/lunique(plot.df$variable)), lunique(plot.df$variable))
plot.df$lab <- plot.df$LEVEL_NAME
#plot.df$pd_comp_top10[plot.df$hotspot_type=="no hotspot"] <- NA

plot.df$value[plot.df$variable=="pre_change" & plot.df$value <= (-2000) & !is.na(plot.df$value)] <- (-2000)
plot.df.sub <- plot.df[which(!is.na(plot.df$pd_comp_top10)),]


# cosmetics!!!

# remove non hotspot level names
library(ggtext)
plot.df.sub$LEVEL_NAME[!plot.df.sub$pd_comp_top10] <- NA

p.spectrum <- ggplot(plot.df, aes(y=value, x=x))+
  geom_area(fill="grey85")+
  geom_bar(data=plot.df.sub, aes(x=x, y=value, fill=pd_comp_top10, col=pd_comp_top10), 
           stat="identity", width=.9)+
  xlab("rank")+
  scale_color_manual("PD top 10", values=c("grey70", bc[1]))+
  scale_y_continuous("Value")+
  scale_fill_manual("PD top 10", values=c("grey70", bc[1]))+
  facet_wrap("variable", scales="free", ncol=2, strip.position="top",
             labeller=as_labeller(c("deforestation2"="deforestation",
                                    "hfp_mean"="human footprint", "pre_change"="precipitation change until 2070",
                                    "mat_change"="temperature change until 2070")))+
  coord_cartesian(xlim=c(0,367))+
  geom_richtext(data=plot.df.sub, aes(x=x, y=value, label=LEVEL_NAME, col=pd_comp_top10),
                size=4, angle=90, hjust=-0.05, label.padding=unit(.05,"mm"),
                show.legend=F, label.color=NA, fill=NA)+
  theme_classic()+
  theme(strip.background=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks=element_blank(), 
        panel.grid=element_blank(),
        legend.position="none", 
        legend.background = element_blank(), 
        plot.background = element_blank(),
        strip.text=element_text(size=10))
p.spectrum
ggsave("figures/fig5.pdf", width=8, height=6, units = "in", dpi = 300, bg = "white")
ggsave("figures/fig5.png", width=8, height=6, units = "in", dpi = 300, bg = "white")














# Scaling effects area + SR -----------------------------------------------

## PD ----
stand_fun <- function(x){(x-mean(x))/sd(x)}

tmp <- st_drop_geometry(shp[,c("area", "PD_obs", "richness", "LEVEL_NAME")])
tmp$PD_richness <- tmp$PD/tmp$richness
names(tmp) <- c("area", "PD", "SR", "LEVEL_NAME", "PD_richness")

tmp[,c(2,3)] <- apply(tmp[,c(2,3)], 2, normalized)
setDT(tmp)
# reshape
tmp <- data.table::melt(tmp, id.var=c("area", "LEVEL_NAME"))

# plot

# (p1 <- ggplot(tmp[grep("sesPD", tmp$variable),], aes(x=area, y=value, col=variable))+
#   geom_point(alpha=0.1)+
#   geom_smooth(se=F, method="lm", lwd=0.7)+
#   scale_x_continuous(trans="log10", #breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = math_format(
#                       format = function(x){number(log10(x), accuracy = 1)}))+
#   labs(y=expression('PD'[std]))+
#   scale_color_manual(values=c(bc[1], "grey20"))+
#   theme(legend.position="none"))
# 
# (p2 <- ggplot(shp2, aes(x=richness, y=sesPD))+
#     geom_point(alpha=0.1)+
#     geom_smooth(se=F, method="lm", lwd=0.7, col="grey20")+
#     scale_x_continuous(trans="log10", 
#                        labels = math_format(
#                          format = function(x){number(log10(x), accuracy = 1)}))+
#     scale_color_manual(values=c("grey20"))+
#     labs(y=expression('PD'[std]),
#          x="species richness")+
#     theme(legend.position="none"))

(pempty <- ggplot(tmp[grep("PD_richness", tmp$variable),], aes(x=as.numeric(area), y=value, col=variable))+
    geom_point(alpha=0)+ theme_void()+ theme(legend.position="none"))

(p3 <- ggplot(tmp[grep("^PD$|^SR", tmp$variable),], aes(x=as.numeric(area), y=value, col=variable))+
  geom_point(alpha=0.1)+
  scale_x_continuous(trans="log10", #breaks = trans_breaks("log10", function(x) 10^x),
                     labels = math_format(
                       format = function(x){number(log10(x), accuracy = 1)}))+
  ylab("standardised value")+
  scale_color_manual(values=c(bc[1], "grey50"))+
  geom_smooth(se=F, method="lm", alpha=0.5, lwd=.7)+
  scale_y_log10()+
  xlab("Area")+
  theme_classic()+
  theme(legend.position=c(0.16, 0.9),
        legend.background=element_blank(), 
        legend.spacing.y=unit(1,"mm"),
        legend.title=element_blank(),
        legend.key.width=unit(3, "mm"),
        strip.background=element_blank()))

plot_grid(p3, pempty, labels=c("A","B"), label_fontface="plain")
ggsave("figures/scaling_log2.pdf", height=3.4, width=5, unit="in", dpi=300)
ggsave("figures/scaling_log2.png", height=2.5, width=5, unit="in", dpi=300, bg="white")




## PD regression coefficients ----
shp$richness_norm_log <- log(normalized(shp$richness))
shp$PD_obs_norm_log <- log(normalized(shp$PD_obs))

tmp <- st_drop_geometry(shp[,c("richness_norm_log", "area")])
tmp <- tmp[!is.infinite(rowSums(tmp)),]
m_SR <- lm(data=tmp, richness_norm_log~log(area))

tmp <- st_drop_geometry(shp[,c("PD_obs_norm_log", "area")])
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
huxtable::quick_html(file="figures/PD_regression_table", tab, open=F) 



# genus size in plants ----------
tmp = fread("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt")
sub = tmp[taxon_status=="Accepted" & taxon_rank=="Species", ]
res = data.table(species_number = tapply(sub$plant_name_id, sub$genus, lunique),
                 genus = names(tapply(sub$plant_name_id, sub$genus, lunique))) 
ggplot(res)+
  geom_histogram(aes(x=species_number))+
  scale_x_continuous(trans="log1p", breaks=c(1,10,100,1000,5000))






# Maximum genus size 	173 
max(res$species_number) # 3054
# Species in monotypic genera 	9.8%
sum(res$species_number[res$species_number==1])/sum(res$species_number) # 0.01
# Proportion of genera monotypic 	43.3% 
length(res$genus[res$species_number==1])/nrow(res) # 28%


# SI data table -----------------------------
rem = grep("_3_CO|_complementarity|_added|area|_added_perc|norm_log|\\.[0-9]|hotspot_coverage", names(shp))
tabsi = st_drop_geometry(shp[, -rem])
tabsi <- tabsi[order(tabsi$LEVEL_NAME),]
View(tabsi)
fwrite(tabsi, "data/data_summary_table_global_PD.csv")













# ## PE 
# tmp <- st_drop_geometry(shp2[,c("area", "PE_obs", "richness", "sesPE", "LEVEL_NAME")])
# names(tmp) <- c("area", "PE", "SR", "PEstd", "LEVEL_NAME")
# 
# tmp[,c(2,3)] <- apply(tmp[,c(2,3)], 2, normalized)
# setDT(tmp)
# # reshape
# tmp <- data.table::melt(tmp, id.var=c("area", "LEVEL_NAME"))
# bc <- c("#52548D", "#C57391", "#EFB984") # PD, PD+PE, PE
# (p12 <- ggplot(tmp[grep("PEstd", tmp$variable),], aes(x=area, y=value, col=variable))+
#     geom_point(alpha=0.1)+
#     geom_smooth(se=F, method="lm", lwd=0.7)+
#     scale_x_continuous(trans="log", 
#                        labels = math_format(
#                          format = function(x){number(log10(x), accuracy = 1)}))+
#     ylab(expression('PE'[std]))+
#     scale_color_manual(values=c(bc[1], "grey20"))+
#     theme(legend.position="none"))
# 
# (p22 <- ggplot(shp2, aes(x=richness, y=sesPE))+
#     geom_point(alpha=0.1)+
#     geom_smooth(se=F, method="lm", lwd=0.7, col="grey20")+
#     scale_x_continuous("species richness", trans="log10", labels = math_format(
#       format = function(x){number(log10(x), accuracy = 1)}))+
#     scale_color_manual(values=c("grey20"))+
#     labs(y=expression('PE'[std]),
#          x="species richness")+
#     theme(legend.position="none"))
# 
# (pempty <- ggplot(tmp[grep("PD_richness", tmp$variable),], aes(x=area, y=value, col=variable))+
#     geom_point(alpha=0)+ theme_void()+ theme(legend.position="none"))
# 
# (p32 <- ggplot(tmp[grep("^PE$|^SR", tmp$variable),], aes(x=area, y=value, col=variable))+
#     geom_point(alpha=0.1)+
#     scale_x_continuous(trans="log", 
#                        labels = math_format(
#                          format = function(x){number(log10(x), accuracy = 1)}))+
#     ylab("standardised value")+
#     scale_color_manual(values=c(bc[2], "grey20"))+
#     geom_smooth(se=F, method="lm", alpha=0.5, lwd=.7)+
#     scale_y_log10()+
#     theme(legend.position=c(0.12, 0.8),
#           legend.background=element_blank(), 
#           legend.spacing.y=unit(1,"mm"),
#           legend.title=element_blank(),
#           legend.key.width=unit(3, "mm"),
#           strip.background=element_blank()))
# 
# plot_grid(plot_grid(p32, pempty, labels=c("A","B"),label_y=1,  label_size=11, label_fontface="plain"),
#           plot_grid(p12, p22, labels=c("C","D"),label_y=1.02,  label_size=11, label_fontface="plain"), 
#           ncol=1, label_size=11, label_fontface="plain")
# 
# ggsave("figures/scaling_PE.png", height=3.4, width=3.4, unit="in", dpi=300)
# 
# 
# 
# ##  PE regression coefficients
# shp2$richness_norm_log <- log(normalized(shp2$richness))
# shp2$PE_obs_norm_log <- log(normalized(shp2$PE_obs))
# 
# tmp <- st_drop_geometry(shp2[,c("richness_norm_log", "area")])
# tmp <- tmp[!is.infinite(rowSums(tmp)),]
# m_SR <- lm(data=tmp, richness_norm_log~log(area))
# 
# tmp <- st_drop_geometry(shp2[,c("PE_obs_norm_log", "area")])
# tmp <- tmp[!is.infinite(rowSums(tmp)),]
# m_PE <- lm(data=tmp, PE_obs_norm_log~log(area))
# 
# mlist <- list(m_SR, m_PE)
# names(mlist) <- c("SR ~ area", "PE ~ area")
# (tab <- huxtable::huxreg(mlist,
#                          error_format = "({std.error})",
#                          error_pos = "below",
#                          number_format = "%.2f",
#                          align = ".",
#                          stars = c(`***` = 0.001),
#                          statistics = c(R2 = "r.squared")))
# huxtable::font_size(tab) <- 13
# huxtable::quick_html(file = "PE_model.html", tab, open = F)
# 
# 
# save(list=c("shp", "shp2"), file="data/workspace_point0.RData")











# old

# # Spatial correlations 
# 
# # sesPD + sesPE VS taxonomic measures
# # reduce to relevant variables:
# dat <- shp[,c("LEVEL3_COD", "richness", "WE", "sesPD", "sesPE")]
# names(dat)
# dat <- na.omit(dat)
# 
# nb <- spdep::poly2nb(dat, row.names = dat$LEVEL3_COD)
# col.W <- nb2listw(nb, style="W", zero.policy = TRUE)
# lee.dat <- st_drop_geometry(dat[,-which(names(dat)=="LEVEL3_COD")])
# 
# les <- apply(lee.dat[,!names(lee.dat)%in%"sesPD"], 2, lee.test, y=lee.dat$sesPD, 
#              listw=col.W, zero.policy = TRUE, alternative="two.sided")
# 
# 
# cor.test(lee.dat$sesPD, lee.dat$sesPE, method="p")
# lee.test(lee.dat$sesPD, lee.dat$sesPE, listw=col.W, zero.policy = TRUE, alternative="two.sided")
# lee.mc(lee.dat$sesPD, lee.dat$sesPE, listw=col.W, zero.policy = TRUE, alternative="greater", nsim=999)
# 
# cor.test(lee.dat$sesPD, lee.dat$richness, method="p")
# lee.test(lee.dat$sesPD, lee.dat$richness, listw=col.W, zero.policy = TRUE, alternative="two.sided")
# lee.mc(lee.dat$sesPD, lee.dat$richness, listw=col.W, zero.policy = TRUE, alternative="less", nsim=999)
# 
# cor.test(lee.dat$sesPE, lee.dat$WE, method="p")
# lee.test(lee.dat$sesPE, lee.dat$richness, listw=col.W, zero.policy = TRUE, alternative="two.sided")
# lee.mc(lee.dat$sesPE, lee.dat$richness, listw=col.W, zero.policy = TRUE, alternative="greater", nsim=999)

# (pd_map <- ggplot(shp) + 
#     geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
#     geom_sf(data=shp, aes(fill=PD_obs),lwd=0.25/.pt, col="gray95") + 
#     #    geom_sf(data=thicc_lines, lwd=1, aes(col=PD_obs), show.legend=F)+
#     # geom_point(data=thicc_lines, aes(color = PD_obs, geometry = geometry),
#     #            stat = "sf_coordinates")+
#     #    scale_color_manual("PD", palette="BuGn", direction=1, values=c(0,upl), trans = "sqrt")+
#     scale_fill_gradient("PD", low="white", high=bc[1])+
#     #    scale_color_gradient("PD", low="white", high=bc[1])+
#     # scale_colour_scico("PD", palette = "batlow", trans = "sqrt", begin = lcol, end = sqrt(ucol))+
#     # scale_fill_scico("PD", palette = "batlow", trans="sqrt")+ #,
#     theme_void()+
#     theme(legend.position = c(0.2, 0.3),
#           legend.key.height = unit(6,"mm"),
#           legend.key.width = unit(4,"mm"),
#           legend.background = element_blank(),
#           legend.key = element_blank(),
#           legend.text.align=1,
#           panel.background = element_blank(),
#           panel.border = element_blank(),
#           text = element_text(size = 10))+
#     xlab(" ")+
#     coord_sf(expand=F, datum=NULL) # needed to keep sf from generating graticule (this would fail)
# )

# # lcol <- 1-min(thicc_lines$sesPD)/(min(shp2$sesPD)-max(shp2$sesPD))
# # ucol <- max(thicc_lines$sesPD)/max(shp2$sesPD)
# # # new lower limit for second scale:
# # lol = -1/(1-lcol) # lower limit
# (pd_ses_map <- ggplot(shp2) + 
#     geom_sf(aes(fill=sesPD),lwd=0, col=NA) + 
#     geom_sf(data=thicc_lines, lwd=1, aes(col=sesPD), show.legend=F)+
#     scale_color_distiller(bquote("PD"[std]), palette="BuGn", direction=1, values=c(lol,1))+
#     scale_fill_distiller(bquote("PD"[std]), palette="BuGn", direction=1)+
#     # scale_colour_scico("PDstd", palette = "batlow", trans = "sqrt", begin = lcol, end = sqrt(ucol))+
#     # scale_fill_scico("PDstd", palette = "batlow", trans="sqrt")+ #,
#     theme_void()+coord_sf(expand=F, datum=NULL)+
#     theme(legend.position = c(0.2, 0.3),
#           legend.key.height = unit(6,"mm"),
#           legend.key.width = unit(4,"mm"),
#           legend.background = element_blank(),
#           legend.key = element_blank(),
#           legend.text.align=1,
#           panel.background = element_blank(),
#           panel.border = element_blank(),
#           text = element_text(size = 10),
#     )+
#     xlab(" ")
# )
# 
# #lcol <- min(thicc_lines$sesPE+abs(min(shp2$sesPE)), na.rm=T)/diff(range(shp2$sesPE)) 
# lcol <- min(thicc_lines$sesPE, na.rm=T)/diff(range(shp2$sesPE))
# ucol <- max(thicc_lines$sesPE, na.rm=T)/diff(range(shp2$sesPE))
# # new limits for second scale:
# lol = min(shp2$sesPE, na.rm=T) / min(thicc_lines$sesPE, na.rm=T) # lower limit
# ul = 1 / ucol # upper limit
# # adjusting lower and upper limits at the same time is tricky ^____^
# 
# (pe_ses_map <- ggplot(shp2) + 
#     geom_sf(aes(fill=sesPE),lwd=0, col=NA) + 
#     geom_sf(data=thicc_lines, lwd=1, aes(col=sesPE), show.legend=F)+
#     scale_color_distiller(bquote("PE"[std]), palette="BuGn", direction=1, values=c(-lol-1,ul))+
#     scale_fill_distiller(bquote("PE"[std]), palette="BuGn", direction=1)+
#     # scale_colour_scico("PEstd", palette = "batlow", trans = "sqrt", begin = lcol, end = sqrt(ucol))+
#     # scale_fill_scico("PEstd", palette = "batlow", trans="sqrt")+
#     theme_void()+coord_sf(expand=F, datum=NULL)+
#     theme(legend.position = c(0.2, 0.3),
#           legend.key.height = unit(6,"mm"),
#           legend.key.width = unit(4,"mm"),
#           legend.background = element_blank(),
#           legend.text.align = 1,
#           legend.key = element_blank(),
#           panel.background = element_blank(),
#           panel.border = element_blank(),
#           text = element_text(size = 10),
#     )+
#     xlab(" ")
# )
# 
# 
# lcol <- min(thicc_lines$PDE, na.rm=T)/max(shp$PDE, na.rm=T)
# ucol <- max(thicc_lines$PDE, na.rm=T)/max(shp$PDE, na.rm=T)
# upl = 1/ucol # upper limit

# # Conservation hotspot coverage 
# ggplot(shp2) + 
#   geom_sf(aes(col=hotspot_coverage), show.legend=T)+
#   geom_sf(aes(fill=hotspot_coverage), col="grey80", lwd=.1)+
#   # scale_fill_gradient("hotspot \ncoverage", low="white", high="red")+
#   # scale_color_gradient(low="white", high="red")+
#   scale_color_distiller("hotspot \ncoverage", palette="BuGn", direction=1)+
#   scale_fill_distiller("hotspot \ncoverage", palette="BuGn", direction=1)+
#   
#   theme_void()+
#   theme(legend.position = c(0.2, 0.3),
#         legend.key.height = unit(6,"mm"),
#         legend.background = element_blank(),
#         legend.key = element_blank(),
#         legend.title=element_text(size=9))+
#   coord_sf(expand=F)
# ggsave("figures/hotspot_coverage.png", units="in", dpi=300, width=7, height=4.1)


