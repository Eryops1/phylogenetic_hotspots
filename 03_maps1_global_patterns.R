
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
gc()
library(sf)
library(ggplot2)
library(cowplot)
#library(beepr)
library(biscale)
library(spdep)
library(scico)
#library(ggpattern)
library(data.table)
library(terra)
#library(treemapify)
if(!dir.exists("figures"))dir.create("figures")
source("99_functions.R")
library(extrafont)
font_import() # takes 5 minutes if loading all, select wisely
loadfonts()
# DejaVu Sans Condensed?
theme_set(theme_bw()+theme(text=element_text(size=7, family="Helvetica"), 
                           panel.grid=element_blank()))


# Load data ---------------------------------------------------------------

gallpeters_projection <- "+proj=cea +lon_0=0 +x_0=0 +y_0=0 +lat_ts=45 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
# Gall-Peters. Horizontally compressed version of the Lambert equal-area.
# Standard parallels at 45°N/S. Aspect ratio of ~1.6. Similar is Balthasar
# projection with standard parallels at 50°N/S.
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

my_projection <- gallpeters_projection
shp <- readRDS("fin_shp.rds")

names(shp)[grep("SES\\.PD", names(shp))] <- "sesPD"
names(shp)[grep("SES\\.PE", names(shp))] <- "sesPE"


# transform projection
shp <- st_transform(shp, my_projection)

# remove not needed data
shp <- shp[!shp$LEVEL3_COD=="ANT",]
shp <- shp[,-grep("obs_p|obs_rank|reps|LEVEL2|LEVEL1|LEVEL_3_CO|ID|\\.3|_rw|CONTI|REGION|AvTD|TTD|mpd", names(shp))]
names(shp)<- gsub("\\.1", "_mean", names(shp))
names(shp)<- gsub("\\.2", "_sd", names(shp))

# Hotspot shapefile
h <- st_read("hotspots_fixed.gpkg")
h <- st_wrap_dateline(h, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
h <- st_transform(h, my_projection)




source("03_island_definition.R")
#source("03_area_standardisation.R")

# run this to add island + continental unit classification, and chose here if
# you wish to standardize for area. Maps for area standardisation are produced
# and for hotspots with those standardised measures

# Plotting settings
min.area <- 6e+9
thicc_lines <- shp[which(shp$area<min.area),]
shp2 <- shp[!is.na(shp$sesPD),]





# Global patterns  ------------------------------------------------

lcol <- min(thicc_lines$PD_obs)/max(shp$PD_obs)
ucol <- max(thicc_lines$PD_obs)/max(shp$PD_obs)
(pd_map <- ggplot(shp) + 
    geom_sf(data=shp, aes(fill=PD_obs),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=PD_obs), show.legend=F)+
    scale_colour_scico("PD", palette = "batlow", trans = "sqrt", 
                            begin = lcol, end = sqrt(ucol))+
    theme_void()+
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
    theme_void()+
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
# Simple SR 
lcol <- min(thicc_lines$richness)/max(shp$richness)
ucol <- max(thicc_lines$richness)/max(shp$richness)
(sr_map <- ggplot(shp) + 
    geom_sf(aes(fill=richness),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=richness), show.legend=F)+
    scale_colour_scico("SR", palette="batlow", trans = "sqrt", 
                           begin = lcol, end = sqrt(ucol))+
    scale_fill_scico("SR", palette="batlow", trans = "sqrt")+ #, 
    theme_void()+coord_sf(expand=F)+
    theme(legend.position = c(0.2, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )
)


# lcol <- min(thicc_lines$sdensity)/max(shp$sdensity)
# ucol <- max(thicc_lines$sdensity)/max(shp$sdensity)
# transform <- "log"
# (sd_map <- ggplot(shp) + 
#     geom_sf(aes(fill=sdensity),lwd=0, col=NA) + 
#     geom_sf(data=thicc_lines, lwd=1, aes(col=sdensity), show.legend=F)+
#     scale_colour_scico("", palette="batlow", trans = transform, 
#                        begin = lcol, end=ucol)+
#     scale_fill_scico("Species density", palette="batlow", trans=transform)+ #,breaks = c(0, 0.01, 0.1, 1, 10), labels = c("","0.01","0.1","1", "10")
#     theme_void()+
#     theme(legend.position = c(0.22, 0.3),
#           legend.key.height = unit(6,"mm"),
#           legend.background = element_blank(),
#           legend.key = element_blank(),
#           panel.background = element_blank(),
#           panel.border = element_blank(),
#           text = element_text(size = 10)
#     )
# )
# lcol <- min(thicc_lines$PD_SR)/max(shp$PD_SR)
# ucol <- max(thicc_lines$PD_SR)/max(shp$PD_SR)
# transform <- "log"
# (PD_SR_map <- ggplot(shp[shp$richness>1000,]) + 
#     geom_sf(aes(fill=PD_SR),lwd=0, col=NA) + 
#     geom_sf(data=thicc_lines[thicc_lines$richness>1000,], lwd=2, aes(col=PD_SR), show.legend=F)+
#     scale_colour_scico("", palette="batlow", trans = transform, 
#                        begin = lcol, end=ucol)+
#     scale_fill_scico("PD/SR", palette="batlow", trans=transform)+ #,breaks = c(0, 0.01, 0.1, 1, 10), labels = c("","0.01","0.1","1", "10")
#     theme_void()+
#     theme(legend.position = c(0.22, 0.3),
#           legend.key.height = unit(6,"mm"),
#           legend.background = element_blank(),
#           legend.key = element_blank(),
#           panel.background = element_blank(),
#           panel.border = element_blank(),
#           text = element_text(size = 10)
#     )
# )
# (nrow(shp[shp$richness>1000,])/100)*5 # 7 bot countries
# tmp <- shp[shp$richness>1000,]
# tmp$LEVEL_NAME[order(tmp$PD_SR, decreasing=T)][1:14]



# simple WE
thicc_lines <- shp2[which(shp2$area<min.area),]
lcol <- min(thicc_lines$WE)/max(shp$WE)
ucol <- max(thicc_lines$WE)/max(shp$WE)
(we_map <- ggplot(shp2) + 
    geom_sf(aes(fill=WE),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=WE), show.legend=F)+
    scale_colour_scico("WE", palette="batlow", trans = "sqrt", 
                           begin = lcol, end = sqrt(ucol))+
    scale_fill_scico("WE", palette="batlow", trans = "sqrt")+ #, 
    theme_void()+coord_sf(expand=F)+
    theme(legend.position = c(0.2, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10)
    )+
    xlab(" ")
)

# lcol <- min(thicc_lines$WE_s)/max(shp$WE_s)
# ucol <- max(thicc_lines$WE_s)/max(shp$WE_s)
# (we_s_map <- ggplot(shp2) + 
#     geom_sf(aes(fill=WE_s),lwd=0, col=NA) + 
#     geom_sf(data=thicc_lines, lwd=1, aes(col=WE_s), show.legend=F)+
#     scale_colour_scico("WE_s", palette="batlow", trans = "sqrt", 
#                        begin = lcol, end = sqrt(ucol))+
#     scale_fill_scico("WE_s", palette="batlow", trans="sqrt")+ #, 
#     theme_void()+
#     theme(legend.position = c(0.22, 0.3),
#           legend.key.height = unit(6,"mm"),
#           legend.background = element_blank(),
#           legend.key = element_blank(),
#           panel.background = element_blank(),
#           panel.border = element_blank(),
#           text = element_text(size = 10)
#     )+
#     xlab(" ")
# )


lcol <- 1-min(thicc_lines$sesPD)/(min(shp2$sesPD)-max(shp2$sesPD))
ucol <- max(thicc_lines$sesPD)/max(shp2$sesPD)
(pd_ses_map <- ggplot(shp2) + 
    geom_sf(aes(fill=sesPD),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=sesPD), show.legend=F)+
    scale_colour_scico("sesPD", palette="batlow", begin=lcol, end=ucol)+
    scale_fill_scico("sesPD", palette="batlow")+  
    theme_void()+coord_sf(expand=F)+
    theme(legend.position = c(0.2, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")
)

lcol <- min(thicc_lines$sesPE+abs(min(shp2$sesPE)))/diff(range(shp2$sesPE)) 
ucol <- max(thicc_lines$sesPE)/max(shp2$sesPE)
(pe_ses_map <- ggplot(shp2) + 
    geom_sf(aes(fill=sesPE),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=sesPE), show.legend=F)+
    scale_colour_scico("sesPE", palette="batlow",
                           begin = lcol, end = ucol)+
    scale_fill_scico("sesPE", palette="batlow")+ 
    theme_void()+coord_sf(expand=F)+
    theme(legend.position = c(0.20, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")
)

# RAW
plot_grid(pd_map, pe_map, ncol = 2,
          labels=c("A","B"), label_fontface=1)
ggsave("figures/maps_unstandardised.png", width=12.5, height=3.5, units = "in", dpi = 600, bg = "white")

plot_grid(sr_map+ggtitle("Species richness\n")+theme(plot.title = element_text(hjust = 0.5, size=8)), 
          we_map+ggtitle("Weighted endemism\n")+theme(plot.title = element_text(hjust = 0.5, size=8)), 
          pd_ses_map+ggtitle("Phylogenetic diversity, standardized effect size\n")+theme(plot.title = element_text(hjust = 0.5, size=8)), 
          pe_ses_map+ggtitle("Phylogenetic endemism, standardized effect size\n")+theme(plot.title = element_text(hjust = 0.5, size=8)),
          ncol = 2, labels=c("A","B","C","D"), label_fontface=1, label_fontfamily="Helvetica", 
          scale=1)
ggsave("figures/maps.png", width=12.5, height=7.5, units = "in", dpi = 300, bg = "white")

# plot_grid(sd_map, we_s_map, pd_ses_map, pe_ses_map, ncol = 2, 
#           labels=c("A","B","C","D"), label_fontface=1)
# ggsave("figures/maps_s.png", width=12.5, height=7, units = "in", dpi = 600, bg = "white")



ggplot(shp2) + 
  geom_sf(aes(fill=island_parts),lwd=0, col=NA) + 
  geom_sf(data=thicc_lines, lwd=1, aes(col=island_parts), show.legend=F)+
  scale_colour_scico("island_parts", palette="batlow",
                     begin = lcol, end = ucol, trans="log")+
  scale_fill_scico("island_parts", palette="batlow", trans="log")+ 
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

(island_map <- ggplot(shp2) + 
    geom_sf(aes(fill=factor(area_class)),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=factor(area_class)), show.legend=F)+
    scale_colour_scico_d("Land mass type", palette="batlow", begin=0, end=.5, alpha=.7, na.translate=F)+
    scale_fill_scico_d("Land mass type", palette="batlow", begin=0, end=.5, alpha=.7, na.translate=F)+ 
    theme_void()+
    theme(legend.position = c(0.15, 0.12),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10),
    )+
    xlab(" ")
)
ggsave("figures/islands.png", width=6, height=3.5, units = "in", dpi = 300, bg = "white")










# Myer hotspot proportions ------------------------------------

# shp2.tmp <- shp2
# shp2.tmp$hotspot_coverage[shp2.tmp$hotspot_coverage==0] <- NA
ggplot(shp2) + 
  geom_sf(aes(col=hotspot_coverage), show.legend=F)+
  geom_sf(aes(fill=hotspot_coverage), col="grey80", lwd=.1)+
  scale_fill_gradient("hotspot \ncoverage", low="white", high="red")+
  scale_color_gradient(low="white", high="red")+
  theme_void()+
  theme(legend.position = c(0.18, 0.3),
        legend.key.height = unit(6,"mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title=element_text(size=9))+
  coord_sf(expand=F)
ggsave("figures/hotspot_coverage.png", units="in", dpi=300, width=7, height=4.1)


## correlation with PD,PE,SR,WE ----

tmp <- shp2[,grep("LEVEL|hotspot_coverage|sesPE$|sesPD$|WE|richness", names(shp2))]
nb <- spdep::poly2nb(tmp, row.names = tmp$LEVEL3_COD)
col.W <- spdep::nb2listw(nb, style="W", zero.policy = TRUE)
tmp <- st_drop_geometry(tmp)
leeHS <- apply(tmp[,grep("sesPD|sesPE|WE|richness", names(tmp))], 2, spdep::lee.test, y=tmp$hotspot_coverage, 
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











# Choropleth PD + SR / PE + WE  ---------------------------

# Top 2.5% 
shp2$LEVEL_NAME[order(shp2$sesPD, decreasing=T)][1:9]
shp2$LEVEL_NAME[order(shp2$sesPE, decreasing=T)][1:9]


dim <- 4
shp2$PE_hotspot <- bi_class(shp2, x = WE, y = sesPE, style = "jenks", dim = dim)$bi_class
shp2$PD_hotspot <- bi_class(shp2, x = richness, y = sesPD, style = "jenks", dim = dim)$bi_class
thicc_lines <- shp2[which(shp2$area<min.area),]


df <- shp2[, grepl("PE_hotspot|PD_hotspot|SES\\.PD$|SES\\.PE$|area$", names(shp2))]
df.mlt <- tidyr::pivot_longer(df, cols=grep("hotspot", names(df)), names_to="group")
thicc.mlt <- df.mlt[which(df.mlt$area<min.area),]

my_pal <- my_pal1

hs <- ggplot() +
  geom_sf(df.mlt, mapping = aes(fill=value), color = NA, size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal=my_pal, dim=dim, na.value="white") +
  geom_sf(data=thicc.mlt, lwd=1, aes(col=value), show.legend=F)+
  bi_scale_color(pal=my_pal, dim=dim, na.value="white")+
  coord_sf(expand=T)+theme_void()+
  facet_wrap(~group, ncol=2)+
  theme(strip.background=element_blank())
  

#### create legends 
legs <- unique(df.mlt$group)
legs.names <- unique(df.mlt$group)
na.x <- c("WE", "SR")
na.y <- c("sesPE", "sesPD")
lab.color <- rep("white", 16)
lab.color[11] <- "red"
for(i in 1:length(legs)){
  tmp <- bi_legend(pal=my_pal, dim=dim, xlab=na.x[i], ylab=na.y[i], size=5)+
    theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text = element_blank(), legend.background=element_blank())+
    annotate("text", x=rep(1:4,each=4), y=rep(1:4,4), 
             label=class_col(df.mlt$value[df.mlt$group==legs[i]]), col=rep("white", 16), size=3, alpha=.7) # no 11 is a hotspot
  assign(legs.names[i], tmp)
}

## draw maps 
ggdraw() + draw_plot(hs, 0, 0, 1, 1) + 
  draw_plot(PD_hotspot, -0.1, .1, .35, .35)+
  draw_plot(PE_hotspot, 0.4, .1, .35, .35)+
  draw_figure_label(position="top.left", "A", size=11)+
  draw_figure_label(position="top", "B", size=11)
ggsave(paste0("figures/choropleth_hotspots.png"),width=10, height=3, units = "in", dpi = 600, bg = "white")




### islands vs continental ----------------------------------

dim <- 4
table(bi_class(shp2[shp2$area_class=="continental",], x = WE, y = sesPE, style = "jenks", dim = dim)$bi_class)
table(bi_class(shp2[shp2$area_class=="continental",], x = richness, y = sesPD, style = "jenks", dim = dim)$bi_class)
table(bi_class(shp2[shp2$area_class=="island",], x = WE, y = sesPE, style = "jenks", dim = dim)$bi_class)
table(bi_class(shp2[shp2$area_class=="island",], x = richness, y = sesPD, style = "jenks", dim = dim)$bi_class)


tmp <- st_drop_geometry(bi_class(shp2[shp2$area_class=="continental",], x = WE, y = sesPE, style = "jenks", dim = dim))
shp2 <- merge(shp2, tmp[,c("LEVEL3_COD", "bi_class")], all.x=T)
names(shp2)[grep("bi_class$", names(shp2))] <- "PE_hotspot_continental"

tmp <- st_drop_geometry(bi_class(shp2[shp2$area_class=="continental",], x = richness, y = sesPD, style = "jenks", dim = dim))
shp2 <- merge(shp2, tmp[,c("LEVEL3_COD", "bi_class")], all.x=T)
names(shp2)[grep("bi_class$", names(shp2))] <- "PD_hotspot_continental"

tmp <- st_drop_geometry(bi_class(shp2[shp2$area_class=="island",], x = WE, y = sesPE, style = "jenks", dim = dim))
shp2 <- merge(shp2, tmp[,c("LEVEL3_COD", "bi_class")], all.x=T)
names(shp2)[grep("bi_class$", names(shp2))] <- "PE_hotspot_island"

tmp <- st_drop_geometry(bi_class(shp2[shp2$area_class=="island",], x = richness, y = sesPD, style = "jenks", dim = dim))
shp2 <- merge(shp2, tmp[,c("LEVEL3_COD", "bi_class")], all.x=T)
names(shp2)[grep("bi_class$", names(shp2))] <- "PD_hotspot_island"


thicc_lines <- shp2[which(shp2$area<min.area),]

df <- shp2[, grepl("continental|island|area", names(shp2))]
df.mlt <- tidyr::pivot_longer(df, cols=grep("hotspot", names(df)), names_to="group")
thicc.mlt <- df.mlt[which(df.mlt$area_class=="island"),]

my_pal <- my_pal1

group.new <- c("continental", "continental", "islands", "islands")
names(group.new) <- unique(df.mlt$group)
hs_islands <- ggplot() +
  geom_sf(df.mlt, mapping = aes(fill=value), color = NA, size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal=my_pal, dim=dim, na.value="grey95") +
  geom_sf(data=thicc.mlt, lwd=2, aes(col=value), show.legend=F)+
  bi_scale_color(pal=my_pal, dim=dim, na.value="grey95")+
  coord_sf(expand=F)+
  facet_wrap(~group, ncol=2, labeller=labeller(group = group.new))+
  theme(strip.background=element_blank())


#### create legends 
legs <- unique(df.mlt$group)
legs.names <- unique(df.mlt$group)
na.x <- c("WE", "SR")
na.y <- c("sesPE", "sesPD")
for(i in 1:length(legs)){
  tmp <- bi_legend(pal=my_pal, dim=dim, xlab=na.x[i], ylab=na.y[i], size=5)+ # 
    theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text = element_blank(), legend.background=element_blank())
#    annotate("text", x=rep(1:4,each=4), y=rep(1:4,4), 
#             label=class_col(df.mlt$value[df.mlt$group==legs[i]]), col="white", size=2, alpha=.7)
  assign(legs.names[i], tmp)
}

## draw maps 
ggdraw() + draw_plot(hs_islands, 0, 0, 1, 1) + 
  draw_plot(PD_hotspot_continental, 0.02, .55, .15, .12)+
  draw_plot(PE_hotspot_continental, 0.02, .065, .15, .12)

ggsave(paste0("figures/choropleth_hotspots_islands_vs_cont.png"), height=6, width=9.545, unit="in", dpi=300)





# get number myers hotspots covered:
tmp <- shp2[shp2$PD_hotspot=="3-3" | shp2$PE_hotspot %in% c("4-3", "3-4", "4-4"), ]
s <- as(st_geometry(tmp), "Spatial")
m3 <- as(st_geometry(h), "Spatial")

res <- intersect(m3, s) 
# first = myer hotspot, second=evol hotspot
myer_number <- gsub(" [0-9]{1,2}", "", names(res))
evol_number <- gsub("[0-9]{1,2} ", "", names(res))
unique(h$NAME[as.numeric(myer_number)])
unique(tmp$LEVEL_NAME[as.numeric(evol_number)])



# Scaling effects area + SR -----------------------------------------------

## standardize ----
stand_fun <- function(x){(x-mean(x))/sd(x)}

tmp <- st_drop_geometry(shp2[,c("area", "PD_obs", "richness", "sesPD", "LEVEL_NAME")])
tmp$PD_richness <- tmp$PD/tmp$richness
names(tmp) <- c("area", "PD", "SR", "sesPD", "LEVEL_NAME", "PD_richness")

tmp[,c(2,3)] <- apply(tmp[,c(2,3)], 2, normalized)
  # reshape
tmp <- data.table::melt(tmp, id.var=c("area", "LEVEL_NAME"))

## plot ----
library(scales)
bc <- c("#52548D", "#C57391", "#EFB984") # PD, PD+PE, PE
(p1 <- ggplot(tmp[grep("sesPD", tmp$variable),], aes(x=area, y=value, col=variable))+
  geom_point(alpha=0.1)+
  geom_smooth(se=F, method="lm", lwd=0.7)+
  scale_x_continuous(trans="log", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = math_format(
                       format = function(x){number(log10(x), accuracy = 1)}))+
  ylab("sesPD")+
  scale_color_manual(values=c(bc[1], "grey20"))+
  theme(legend.position="none"))
(p2 <- ggplot(tmp[grep("PD_richness", tmp$variable),], aes(x=area, y=value, col=variable))+
    geom_point(alpha=0.1)+
    geom_smooth(se=F, method="lm", lwd=0.7)+
    scale_x_continuous(trans="log", breaks = trans_breaks("log10", function(x) 10^x),
                       labels = math_format(
                         format = function(x){number(log10(x), accuracy = 1)}))+
    ylab("PD/species richness")+
    scale_color_manual(values=c("grey20"))+
    theme(legend.position="none"))
(pempty <- ggplot(tmp[grep("PD_richness", tmp$variable),], aes(x=area, y=value, col=variable))+
    geom_point(alpha=0)+ theme_void()+ theme(legend.position="none"))
(p3 <- ggplot(tmp[grep("^PD$|^SR", tmp$variable),], aes(x=area, y=value, col=variable))+
  geom_point(alpha=0.1)+
  scale_x_continuous(trans="log", breaks = trans_breaks("log10", function(x) 10^x),
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
plot_grid(plot_grid(p3, pempty, labels=c("A",""),label_y=1,  label_size=11, label_fontface="plain"),
          plot_grid(p1, p2, labels=c("B","C"),label_y=1.02,  label_size=11, label_fontface="plain"), 
          ncol=1, labels=c("A","", ""), label_size=11, label_fontface="plain")

ggsave("figures/scaling_log.pdf", height=3.4, width=3.4, unit="in", dpi=300)

## Table regression coefficients ----
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
huxtable::quick_html(tab) 


# schematic figure random tree ----
# dev.off()
# tree <- rtree(12, rooted = T, tip.label = NULL, br = runif, equiprob=T)
# tree <- phytools::force.ultrametric(tree)
# library(ape)
# ape::plot.phylo(tree, direction="upwards", show.tip.label=F, type="ti", use.edge.length=T)
# #ape::plot.phylo(rtopology(10, rooted = FALSE, tip.label = NULL, br = runif))
# #("figures/random_tree.pdf", height=2, width=3, unit="in")
# plot(x=c(1,2,3,4,5,6,7), y=c(10,26,37,47,55,61,66), ylim=c(0,70), type="b")
# points(x=c(1,2,3,4,5,6,7), y=c(4,8,12,16, 20,24,28), col="red", type="b")



# TABLE 1 -------------------------------------------------
# check
table(shp2$PD_hotspot)
table(shp2$PE_hotspot)
table(shp2$PD_hotspot_island)
table(shp2$PE_hotspot_island)
table(shp2$PD_hotspot_continental)
table(shp2$PE_hotspot_continental)

# define
pdspots <- "3-3"
pespots <- c("3-4", "4-3", "4-4")
pdspots_island <- c("4-4", "3-4", "4-3")
pespots_island <- c("3-4", "4-3", "4-4", "3-3")
pdspots_cont <- c("3-3", "3-4")
pespots_cont <- pespots

#assign
shp2$hotspot_type <- NA
shp2$hotspot_type[shp2$PD_hotspot==pdspots] <- "PD"
shp2$hotspot_type[shp2$PE_hotspot%in%pespots] <- "PE"
shp2$hotspot_type[shp2$PE_hotspot%in%pespots&shp2$PD_hotspot==pdspots] <- "PD & PE"

shp2$hotspot_type_island <- NA
shp2$hotspot_type_island[shp2$PD_hotspot_island %in% pdspots_island] <- "PD"
shp2$hotspot_type_island[shp2$PE_hotspot_island%in%pespots_island] <- "PE"
shp2$hotspot_type_island[shp2$PE_hotspot_island%in%pespots_island&shp2$PD_hotspot_island%in%pdspots_island] <- "PD & PE"

shp2$hotspot_type_continental <- NA
shp2$hotspot_type_continental[shp2$PD_hotspot_continental %in% pdspots_cont] <- "PD"
shp2$hotspot_type_continental[shp2$PE_hotspot_continental%in%pespots_cont] <- "PE"
shp2$hotspot_type_continental[shp2$PE_hotspot_continental%in%pespots_cont&shp2$PD_hotspot_cont%in%pdspots_cont] <- "PD & PE"



###### All -----
tabs <- st_drop_geometry(shp2)
tabs <- tabs[!is.na(tabs$hotspot_type), ]
tmp <- apply(tabs[,c("deforestation2", "hfp_mean", "mat_change", "pre_change")], 2, function(x){rank(x)/length(x)})
colnames(tmp) <- paste0(colnames(tmp), "_rank")
tabs <- cbind(tabs, tmp)

tab1 <- tabs[,c("LEVEL_NAME", "hotspot_type", "richness", "deforestation2",
                "hfp_mean", "mat_change", "pre_change",
                "hotspot_coverage")]
tab1[,-c(1,2)] <- round(tab1[,-c(1,2)], 2)
tab1 <- tab1[order(tab1$LEVEL_NAME),]
knitr::kable(tab1, digits=2, format="simple", row.names=F)


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

# merge biomes for clarity?
biomes$'(Sub)Tropical forests' <- rowSums(biomes[,c("(Sub)Tropical moist broadleaf forests",
          "(Sub)Tropical dry broadleaf forests",
          "(Sub)Tropical coniferous forests")])
biomes$'Temperate forests' <- rowSums(biomes[,c("Temperate broadleaf and mixed forests",
                                                "Temperate coniferous forests")])
biomes <- biomes[,-grep("Tropical moist broadleaf forests|Tropical dry broadleaf forests|Tropical coniferous forests|Temperate broadleaf and mixed forests|Temperate coniferous forests", names(biomes))]

# remove biomes with less than 3% for clarity
for(i in 1:nrow(biomes)){
  tmp <- biomes[i,-1]
  tmp[which(tmp<0.03)] <- 0
  biomes[i,-1] <- tmp
}

rowSums(biomes[,-1], na.rm=T)
# scale to 100% to account for minor boundary inaccuracies
biomes[,-1] <- t(apply(biomes[,-1], 1, function(x){x/sum(x, na.rm=T)}) )
rowSums(biomes[,-1], na.rm=T)

hab.df <- tabs[!is.na(tabs$hotspot_type), c("LEVEL_NAME", "LEVEL3_COD")]
hab.df <- merge(hab.df, biomes, all.x=TRUE, by.x="LEVEL3_COD", by.y="country")
hab.mlt <- melt(setDT(hab.df))
hab.mlt <- hab.mlt[order(LEVEL_NAME, value),]
hab.mlt$sort <- rep(length(unique(hab.mlt$variable)):1, length(unique(hab.mlt$LEVEL_NAME)))
# remove empty biomes
hab.mlt <- hab.mlt[hab.mlt$value!=0,]

hab.mlt <- droplevels(hab.mlt)
hab.mlt <- as.data.frame(hab.mlt)
hab.mlt$variable <- factor(hab.mlt$variable, levels=levels(hab.mlt$variable)[c(8,1,9,3,6,2,5,4,7)])

my_pal <- c(scico(8, alpha = .7, begin = 0.1, end = .9, direction = 1, palette = "batlow"), "#dbdbdb")
# ggplot(hab.mlt, aes(area=value, fill=variable)) +
#   geom_treemap()+
#   scale_fill_manual("biome", values=my_pal)+
#   facet_wrap(~LEVEL_NAME, ncol=1, strip.position="left")+
#   theme_void() + theme(strip.background=element_blank(), strip.text.y.left=element_text(angle=0, hjust=1))
# ALT
ggplot(hab.mlt, aes(x=value, y=LEVEL_NAME, fill=variable)) +
  geom_bar(stat="identity")+
  scale_y_discrete(limits=rev)+
  scale_fill_manual("Biome", values=my_pal)+
  guides(fill=guide_legend(nrow=3,byrow=TRUE))
ggsave("figures/table1_biomes.pdf", width=12, height=6, units="in")


###### Continental -----

tabs <- st_drop_geometry(shp2)
tabs <- tabs[!is.na(tabs$hotspot_type_continental), ]
tmp <- apply(tabs[,c("deforestation2", "hfp_mean", "mat_change", "pre_change")], 2, function(x){rank(x)/length(x)})
colnames(tmp) <- paste0(colnames(tmp), "_rank")
tabs <- cbind(tabs, tmp)

tab1 <- tabs[,c("LEVEL_NAME", "hotspot_type_continental", "richness", "deforestation2",
                "hfp_mean", "mat_change", "pre_change",
                "hotspot_coverage")]
tab1[,-c(1,2)] <- round(tab1[,-c(1,2)], 2)
tab1 <- tab1[order(tab1$LEVEL_NAME),]
knitr::kable(tab1, digits=2, format="simple", row.names=F)

hab.df.continental <- tabs[, c("LEVEL_NAME", "LEVEL3_COD")]
hab.df.continental <- merge(hab.df.continental, biomes, all.x=TRUE, by.x="LEVEL3_COD", by.y="country")
hab.mlt.continental <- melt(setDT(hab.df.continental))
hab.mlt.continental <- hab.mlt.continental[order(LEVEL_NAME, value),]
hab.mlt.continental$sort <- rep(length(unique(hab.mlt.continental$variable)):1, length(unique(hab.mlt.continental$LEVEL_NAME)))
# remove empty biomes
hab.mlt.continental <- hab.mlt.continental[hab.mlt.continental$value!=0,]

hab.mlt.continental <- droplevels(hab.mlt.continental)
hab.mlt.continental <- as.data.frame(hab.mlt.continental)
hab.mlt.continental$variable <- factor(hab.mlt.continental$variable, levels=levels(hab.mlt.continental$variable)[c(7,8,1,3,6,2,5,4)])

# ggplot(hab.mlt.continental, aes(area=value, fill=variable)) +
#   geom_treemap()+
#   scale_fill_manual("biome", values=my_pal)+
#   facet_wrap(~LEVEL_NAME, ncol=1, strip.position="left")+
#   theme_void() + theme(strip.background=element_blank(), strip.text.y.left=element_text(angle=0, hjust=1))
ggplot(hab.mlt.continental, aes(x=value, y=LEVEL_NAME, fill=variable)) +
  geom_bar(stat="identity")+
  scale_y_discrete(limits=rev)+
  scale_fill_manual("Biome", values=my_pal)+
  guides(fill=guide_legend(nrow=3,byrow=TRUE))


###### Islands -----
tabs <- st_drop_geometry(shp2)
tabs <- tabs[!is.na(tabs$hotspot_type_island), ]
tmp <- apply(tabs[,c("deforestation2", "hfp_mean", "mat_change", "pre_change")], 2, function(x){rank(x)/length(x)})
colnames(tmp) <- paste0(colnames(tmp), "_rank")
tabs <- cbind(tabs, tmp)

tab1 <- tabs[,c("LEVEL_NAME", "hotspot_type_island", "richness", "deforestation2",
                "hfp_mean", "mat_change", "pre_change",
                "hotspot_coverage")]
tab1[,-c(1,2)] <- round(tab1[,-c(1,2)], 2)
tab1 <- tab1[order(tab1$LEVEL_NAME),]
knitr::kable(tab1, digits=2, format="simple", row.names=F)

# no need to merge here, load data again
biomes <- readRDS("../DATA/PDiv/biomes_olson_ALL.rds")
names(biomes)[grepl("^X", names(biomes))] <- biome_names
# remove biomes with less than 3% for clarity
for(i in 1:nrow(biomes)){
  tmp <- biomes[i,-1]
  tmp[which(tmp<0.03)] <- 0
  biomes[i,-1] <- tmp
}
rowSums(biomes[,-1], na.rm=T)
# scale to 100% to account for minor boundary inaccuracies
biomes[,-1] <- t(apply(biomes[,-1], 1, function(x){x/sum(x, na.rm=T)}) )
rowSums(biomes[,-1], na.rm=T)

hab.df.island <- tabs[, c("LEVEL_NAME", "LEVEL3_COD")]
hab.df.island <- merge(hab.df.island, biomes, all.x=TRUE, by.x="LEVEL3_COD", by.y="country")
hab.mlt.island <- melt(setDT(hab.df.island))
hab.mlt.island <- hab.mlt.island[order(LEVEL_NAME, value),]
hab.mlt.island$sort <- rep(16:1, length(unique(hab.mlt.island$LEVEL_NAME)))
# remove empty biomes
hab.mlt.island <- hab.mlt.island[hab.mlt.island$value!=0,]

hab.mlt.island <- droplevels(hab.mlt.island)
hab.mlt.island <- as.data.frame(hab.mlt.island)

# ggplot(hab.mlt.island, aes(area=value, fill=variable)) +
#   geom_treemap()+
#   scale_fill_scico_d("biome", palette="batlow", begin=.2, end=.5, alpha=.7)+
#   facet_wrap(~LEVEL_NAME, ncol=1, strip.position="left")+
#   theme_void() + theme(strip.background=element_blank(), strip.text.y.left=element_text(angle=0, hjust=1))
ggplot(hab.mlt.island, aes(x=value, y=LEVEL_NAME, fill=variable)) +
  geom_bar(stat="identity")+
  scale_y_discrete(limits=rev)+
  scale_fill_manual("Biome", values=my_pal)+
  guides(fill=guide_legend(nrow=3,byrow=TRUE))



save.image("workspace_point1.RData")






# Hotspot maps  ------------------------------------------
# SAVEPOINT ------------------------
load("workspace_point1.RData")


## load simplified hotspot map
# ha <- st_read("hotspot_area.gpkg")
# ha_united <- st_union(ha)
# st_write(ha_united, "hotspot_area_united.gpkg")
#ha_united <- st_read("hotspot_area_united.gpkg")
#pm <- st_read("Poly-Micronesia.gpkg")
gc()
ha_united <- st_read("hotspots_2016_1/hotspots_merged_cleaned.gpkg")


# change projection + solve dateline issue
ha_united <- st_wrap_dateline(ha_united, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
#pm <- st_wrap_dateline(pm, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
ha_united <- st_transform(ha_united, my_projection)
#pm <- st_transform(pm, my_projection)


## Hotspot map ----------------------------------------
shp2$hotspot_type <- NA
shp2$hotspot_type[shp2$PD_hotspot==pdspots] <- "PD"
shp2$hotspot_type[shp2$PE_hotspot%in%pespots] <- "PE"
shp2$hotspot_type[shp2$PE_hotspot%in%pespots&shp2$PD_hotspot==pdspots] <- "PD & PE"

diss <- st_read("../DATA/shapefile_bot_countries/level3_dissolved.gpkg")
diss <- st_simplify(diss)

thicc_lines <- shp2[which(shp2$area<6e+08),]

bc <- c("#52548D", "#C57391", "#EFB984")
 shp3 <- st_simplify(shp2)
 thicc_lines3 <- st_simplify(thicc_lines)
 ha_united3 <- st_simplify(ha_united)
 
# ggplot() +
#     geom_sf(shp3, mapping = aes(fill=hotspot_type), color=NA, size=0.1) +
#     geom_sf(data=thicc_lines3, lwd=1, col="grey98", show.legend=F)+
#     geom_sf(data=diss, lwd=0.1, col="black", fill=NA, show.legend=F)+
#     scale_fill_manual("Hotspot\ntype", values=bc, na.value="white")+
#     geom_sf(data=ha_united3, col="grey", fill="grey", alpha=.4, size=0.1)+
#     geom_sf_text(data=ha_united3, aes(label=NAME), size=2)+
#     theme_void()+ coord_sf(expand = F)+
#     theme(legend.position = c(0.1, 0.15),
#       legend.key.height = unit(5,"mm"),
#       legend.background = element_blank(),
#       legend.key = element_blank(),
#       panel.background = element_blank(), 
#       legend.text=element_text(size=7),
#       legend.title=element_text(size=8))#, 
# ggsave("figures/one_hotspots_map.pdf", width=150, height=100, units="mm", dpi=300, bg="white")  
ggplot() +
  geom_sf(shp3, mapping = aes(fill=hotspot_type), color=NA, size=0.1) +
#  geom_sf(data=thicc_lines3, lwd=1, col="grey40", show.legend=F)+
  geom_sf(data=ha_united3, col="grey99", fill="grey", alpha=.3, size=0.05)+
    geom_sf(data=diss, lwd=0.1, col="grey40", fill=NA, show.legend=F)+
  scale_fill_manual("Hotspot type", values=bc, na.value=NA, na.translate=F)+

  geom_sf_text(data=shp3[shp3$hotspot_type!="no hotspot",], aes(label=LEVEL_NAME), size=2, fontface="italic", col="grey30")+
  theme_void()+ coord_sf(ylim=c(-7001217, 8928740))+
  theme(legend.position = c(0.15, 0.15),
        legend.key.height = unit(5,"mm"),
        legend.background = element_rect(fill="white", colour="white", size=3),
#        legend.key = element_blank(),
        panel.background = element_blank(), 
        legend.text=element_text(size=7),
        legend.title=element_text(size=8))
#ggsave("figures/one_hotspots_map_v2.pdf", width=160, height=90, units="mm", dpi=300, bg="white")  



### alt panel layout ----
diss <- st_read("../DATA/shapefile_bot_countries/level3_dissolved.gpkg")
diss <- st_simplify(diss)
labsize <- 2

p1 <- ggplot() +
    geom_sf(shp3, mapping = aes(fill=hotspot_type), color=NA, size=0.1) +
    geom_sf(data=thicc_lines3, lwd=1, col=NA, show.legend=F)+
    geom_sf(data=diss, lwd=0.1, col="black", fill=NA, show.legend=F)+
    scale_fill_manual("Hotspot\ntype", values=bc, na.value="white")+
    geom_sf(data=ha_united3, col=NA, fill="grey80", alpha=0.5, lwd=0)+
  geom_sf_text(data=ha_united3, aes(label=NAME), size=labsize)+
#    geom_sf_label(data=ha_united, aes(label=NAME), size=labsize)+
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


# ggplot() +
#   geom_sf(shp3, mapping = aes(fill=hotspot_type), color=NA, size=0.1) +
#   geom_sf_text(data=shp3, aes(label=LEVEL_NAME), size=4)+
#   theme_void()+ coord_sf(xlim=c(4802430, 10502430), ylim=c(-1328740, 5028740))



## Hotspot type map - ISLANDS ------------------------------------------
pdspots_island <- c("4-4", "3-4", "4-3")
pespots_island <- c("3-4", "4-3", "4-4", "3-3")
pdspots_cont <- c("3-3", "3-4")
pespots_cont <- pespots

shp2$hotspot_type_island <- NA
shp2$hotspot_type_island[shp2$PD_hotspot_island %in% pdspots_island] <- "PD"
shp2$hotspot_type_island[shp2$PE_hotspot_island%in%pespots_island] <- "PE"
shp2$hotspot_type_island[shp2$PE_hotspot_island%in%pespots_island&shp2$PD_hotspot_island%in%pdspots_island] <- "PD & PE"

shp2$hotspot_type_continental <- NA
shp2$hotspot_type_continental[shp2$PD_hotspot_continental %in% pdspots_cont] <- "PD"
shp2$hotspot_type_continental[shp2$PE_hotspot_continental%in%pespots_cont] <- "PE"
shp2$hotspot_type_continental[shp2$PE_hotspot_continental%in%pespots_cont&shp2$PD_hotspot_cont%in%pdspots_cont] <- "PD & PE"

bc <- c("#52548D", "#C57391", "#EFB984")
# geom_sf(shp3, mapping = aes(fill=hotspot_type), color=NA, size=0.1) +
#   #  geom_sf(data=thicc_lines3, lwd=1, col="grey40", show.legend=F)+
#   geom_sf(data=ha_united3, col="grey99", fill="grey", alpha=.3, size=0.05)+
#   geom_sf(data=diss, lwd=0.1, col="grey40", fill=NA, show.legend=F)+
#   scale_fill_manual("Hotspot type", values=bc, na.value=NA, na.translate=F)+
#   
#   geom_sf_text(data=shp3[shp3$hotspot_type!="no hotspot",], aes(label=LEVEL_NAME), size=2, fontface="italic", col="grey30")+

shp2 <- st_simplify(shp2)
hotspots_continental <- ggplot() +
    geom_sf(shp2, mapping = aes(fill=hotspot_type_continental), color=NA, size=0.1) +
    geom_sf(data=thicc_lines, lwd=1, col="white", show.legend=F)+
    scale_fill_manual("Hotspot\ntype", values=bc, na.value="white", na.translate=F)+
    geom_sf(data=ha_united3, col="grey99", fill="grey", alpha=.3, size=0.05)+
  geom_sf(data=diss, lwd=0.1, col="grey40", fill=NA, show.legend=F)+
    theme_void()+ coord_sf(expand = F)+
  geom_sf_text(data=shp2[shp2$hotspot_type_continental!="no hotspot",], aes(label=LEVEL_NAME), size=2, fontface="italic", col="grey30")+
  theme(legend.position = c(0.15, 0.15),
        legend.key.height = unit(5,"mm"),
        legend.background = element_rect(fill="white", colour="white", size=3),
        #        legend.key = element_blank(),
        panel.background = element_blank(), 
        legend.text=element_text(size=7),
        legend.title=element_text(size=8))
hotspots_continental

# recenter map for island spots
# gallpeters_pacific <- "+proj=cea +lon_0=150 +x_0=0 +y_0=0 +lat_ts=45 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
# shp3 <- st_transform(shp2, gallpeters_pacific) 
# thicc_lines3 <- shp3[shp3$area_class=="island",]
# thicc_lines2 <- shp2[shp2$area_class=="island",]

hotspots_islands_main <-ggplot() +  #
    geom_sf(shp2, mapping = aes(fill=hotspot_type_island), color=NA, size=0.1, show.legend=F) +
    geom_sf(thicc_lines2, mapping=aes(col=hotspot_type_island, fill=hotspot_type_island), lwd=1, show.legend=F)+
    scale_fill_manual("Hotspot\ntype", values=bc[-1], na.value="white")+
    scale_color_manual("Hotspot\ntype", values=bc[-1], na.value="white")+
  geom_sf(data=ha_united3, col="grey99", fill="grey", alpha=.3, size=0.05)+
  geom_sf(data=diss, lwd=0.1, col="grey40", fill=NA, show.legend=F)+
  theme_void()+ coord_sf(expand = F)+
  geom_sf_text(data=shp2[shp2$hotspot_type_island!="no hotspot",], aes(label=LEVEL_NAME), size=2, fontface="italic", col="grey30")+
  theme(legend.position = c(0.15, 0.15),
        legend.key.height = unit(5,"mm"),
        legend.background = element_rect(fill="white", colour="white", size=3),
        #        legend.key = element_blank(),
        panel.background = element_blank(), 
        legend.text=element_text(size=7),
        legend.title=element_text(size=8))
hotspots_islands_main

# inset_map <- ggplot() +  
#   geom_sf(shp3, mapping = aes(fill=hotspot_type_island), color=NA, size=0.1, show.legend=F) +
#   geom_sf(thicc_lines3, mapping=aes(col=hotspot_type_island, fill=hotspot_type_island), lwd=1, show.legend=F)+
#   scale_fill_manual("Hotspot\ntype", values=bc[-1], na.value="grey95")+
#   scale_color_manual("Hotspot\ntype", values=bc[-1], na.value="grey95")+
#   geom_sf(data=ha_united, col=NA, fill="grey", alpha=.5)+
#   coord_sf(xlim=c(-4092430, 3092211), ylim=c(-3500000, 4000000))+
#   theme(legend.position = c(0.9, 0.8),
#         legend.key.height = unit(5,"mm"),
#         legend.background = element_blank(),
#         legend.key = element_blank(),
#         panel.background = element_blank(),
#         axis.text=element_blank(), axis.ticks=element_blank(), 
#         plot.margin=unit(c(0,0,0,0),"mm"))
# 
# hotspots_islands <- ggdraw() + draw_plot(hotspots_islands_main, 0, 0, 1, 1) + 
#   draw_plot(inset_map, 0, .1, .7, .8)

#manual tweaking of map
ggsave(plot=hotspots_islands_main, "figures/island_hotspots.pdf")

tmp <- plot_grid(labels="AUTO", hotspots_continental, hotspots_islands_main, label_fontface=1, scale=0.99)
tmp
ggsave(plot=tmp, "figures/one_hotspots_map_islands.pdf", width=10, height=3.2, dpi=300, bg="white")  

#View(st_drop_geometry(shp2[,grep("hotspot|LEVEL_NAME", names(shp2))]))


## Build tiny barplots -----

# DO NOT make the new subset the max values, draws a wrong picture (eg
# Cape povinces having hfp of 100%)
plot.df <- st_drop_geometry(shp2[,c("LEVEL_NAME", "deforestation2", "hfp_mean", 
                                    "pre_change", "mat_change", "PD_hotspot", "PE_hotspot")])
plot.df$PD_hotspot[plot.df$PD_hotspot!=pdspots] <- "n"
plot.df$PD_hotspot[plot.df$PD_hotspot==pdspots] <- "y"
plot.df$PE_hotspot[!plot.df$PE_hotspot%in%pespots] <- "n"
plot.df$PE_hotspot[plot.df$PE_hotspot%in%pespots] <- "y"
plot.hp <- plot.df[plot.df$PD_hotspot=="y" | plot.df$PE_hotspot=="y" ,]

plot.hp$deforestation2 <- plot.hp$deforestation2/max(shp2$deforestation2, na.rm=T)
plot.hp$hfp_mean <- plot.hp$hfp_mean/max(shp2$hfp_mean, na.rm=T)
plot.hp$mat_change <- plot.hp$mat_change/max(shp2$mat_change, na.rm=T)
plot.hp$pre_change <- plot.hp$pre_change/max(shp2$pre_change, na.rm=T)

# plot a facet for each country showing each variable
test <- melt(plot.hp)

bat <- scico(palette="hawaii", n=4, begin=0, alpha=1)
# arrange colors (bat[2] = defore)
bat2 <- c(bat[3], bat[2], bat[4], bat[1])
ggplot(test, aes(fill=variable, x=1, y=value))+
  geom_bar(position=position_dodge(), stat="identity")+
  scale_fill_manual(values=bat2)+
  facet_wrap(~LEVEL_NAME, scales="free")+
  theme_void()+coord_cartesian(expand=F, ylim=c(NA, 1))+
  theme(panel.border=element_rect(fill=NA), 
        strip.text=element_text(size=8),
        plot.background=element_rect("white"))
ggsave("figures/little_bars_alpha09.pdf", width=5, height=5)

save.image("workspace_maps.RData")


# Mittermeier hotspots not touched by our hotspots
ha_united
#Maybe ranking botanical regions by PD/PE (as in fig.5) and marked them by % of
#Myers' hotspots coverage? Would this make sense?
shp6 <- st_drop_geometry(shp2[,c("sesPD", "sesPE", "hotspot_coverage", "LEVEL_NAME")])
plot.df <- shp6[order(shp6$sesPD, shp2$hotspot_coverage),]
plot.df$x <- c(1:352)
#plot.df$lab <- plot.df$LEVEL_NAME
#plot.df$lab[which(plot.df$PD_hotspot!="y"|plot.df$PE_hotspot!="y")] <- NA
#plot.df.sub <- plot.df[which(plot.df$PD_hotspot=="y"|plot.df$PE_hotspot=="y"),]

ggplot(plot.df, aes(x=hotspot_coverage, y=sesPD))+
  geom_point()
ggplot(plot.df, aes(x=x, y=sesPD))+
  geom_area(fill="grey85")+
  geom_bar(data=plot.df, aes(x=x, y=sesPD, fill=hotspot_coverage, col=hotspot_coverage), 
            stat="identity", width=.9)


# Hotspot threats -----------------------------------------------------
# SAVEPOINT ------------------------
## Choropleth threats + PE / PD -------------------------------

load("workspace_maps.RData")

# dim <- 4
# # PE_deforest_bi <- bi_wrap(dat=shp2, x=deforestation2, y=sesPE, style="jenks", dim=4)
# 
# shp2$bi_class_PE_deforest <- bi_class(shp2, x = deforestation2, y = sesPE, style = "jenks", dim = 4)$bi_class
# shp2$bi_class_PD_deforest <- bi_class(shp2, x = deforestation2, y = sesPD, style = "jenks", dim = 4)$bi_class
# shp2$bi_class_PE_hfp <- bi_class(shp2, x = hfp_mean, y = sesPE, style = "jenks", dim = 4)$bi_class
# shp2$bi_class_PD_hfp <- bi_class(shp2, x = hfp_mean, y = sesPD, style = "jenks", dim = 4)$bi_class
# # pre_change fix: that one neg value fucks up clusters, move closer to the rest
# shp2$pre_change[which.min(shp2$pre_change)] <- -2000
# shp2$bi_class_PE_pre_change <- bi_class(shp2, x=pre_change, y=sesPE, style="jenks", dim=4)$bi_class
# shp2$bi_class_PD_pre_change <- bi_class(shp2, x=pre_change, y=sesPD, style="jenks", dim=4)$bi_class
# shp2$bi_class_PE_mat_change <- bi_class(shp2, x=mat_change, y=sesPE, style="jenks", dim=4)$bi_class
# shp2$bi_class_PD_mat_change <- bi_class(shp2, x=mat_change, y=sesPD, style="jenks", dim=4)$bi_class
# thicc_lines <- shp2[which(shp2$area<min.area),]
# 
# 
# #### facet wrap 
# df <- shp2[, grepl("bi_class|SES\\.PD$|SES\\.PE$|area$", names(shp2))]
# df <- df[,!grepl("bi_class$|bi_class_ses$", names(df))]  # remove unwanted columns
# df.mlt <- tidyr::pivot_longer(df, cols=grep("bi_class", names(df)), names_to="group")
# thicc.mlt <- df.mlt[which(df.mlt$area<min.area),]
# 
# my_pal <- my_pal1
# 
# fw <- ggplot() +
#   geom_sf(df.mlt, mapping = aes(fill=value), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal=my_pal, dim=dim, na.value="white") +
#   geom_sf(data=thicc.mlt, lwd=1, aes(col=value), show.legend=F)+
#   bi_scale_color(pal=my_pal, dim=dim, na.value="white")+
#   coord_sf(expand=F)+
#   facet_wrap(~group, ncol=2)+
#   theme(strip.background=element_blank())
# 
# #### legends 
# legs <- unique(df.mlt$group)
# legs.names <- unique(df.mlt$group)
# na <- gsub("bi_class_", "", legs)
# na.x <- gsub("^.*?_", "", na)
# na.y <- paste0("ses", gsub("_.*", "", na))
# for(i in 1:8){
#   tmp <- bi_legend(pal=my_pal, dim=dim, xlab=na.x[i], ylab=na.y[i], size=5)+
#     theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text = element_blank(), legend.background=element_blank())+
#     annotate("text", x=rep(1:4,each=4), y=rep(1:4,4), 
#              label=class_col(df.mlt$value[df.mlt$group==legs[i]]), col="white", size=2, alpha=.7)
#   assign(legs.names[i], tmp)
# }
# 





# Hotspots on the threat spectrum -----------------------
### reshape dataframe for plotting
plot.df <- st_drop_geometry(shp2[,c("LEVEL_NAME", "deforestation2", "hfp_mean", 
                                    "pre_change", "mat_change", "hotspot_type")])
plot.df <- reshape::melt(plot.df)
ypos <- rep(rev(seq(0.2,0.80,0.04)), 4)

library(ggtext)
# alt visualisation
plot.df <- plot.df[order(plot.df$variable, plot.df$value),]
plot.df$x <- rep(1:(nrow(plot.df)/lunique(plot.df$variable)), lunique(plot.df$variable))
plot.df$lab <- plot.df$LEVEL_NAME
#plot.df$lab[which(plot.df$PD_hotspot!="y"|plot.df$PE_hotspot!="y")] <- NA
plot.df.sub <- plot.df[which(!is.na(plot.df$hotspot_type)),]

plot.df$value[plot.df$variable=="pre_change" & plot.df$value <= (-2000) & !is.na(plot.df$value)] <- (-2000)
# cosmetics!!!



# theme_set(theme_bw()+theme(text=element_text(size=10, family="sans"), 
#                            panel.grid=element_blank()))
theme_set(theme_bw()+theme(text=element_text(size=8, family="Helvetica"), 
                           panel.grid=element_blank()))
#font_import(pattern="helvn")
#loadfonts(device="pdf")

p.spectrum <- ggplot(plot.df, aes(y=value, x=x))+
  geom_area(fill="grey85")+
  geom_bar(data=plot.df.sub, aes(x=x, y=value, fill=hotspot_type, col=hotspot_type), 
           stat="identity", width=.9)+
  xlab("rank")+
  facet_wrap(~variable, scales="free")+
  scale_color_manual("Hotspot type", values=bc)+scale_fill_manual("Hotspot type", values=bc)+
  scale_y_continuous("Value")+
  facet_wrap("variable", scales="free", ncol=2, 
             labeller=as_labeller(c("deforestation2"="deforestation",
                                    "hfp_mean"="human footprint", "pre_change"="precipitation change",
                                    "mat_change"="temperature change")))+
  coord_cartesian(expand=F)+
  geom_richtext(data=plot.df.sub, aes(x=x, y=value, label=LEVEL_NAME, col=hotspot_type),
                size=3, angle=90, hjust=-0.05, label.padding=unit(.05,"mm"),
                show.legend=F, label.color=NA, fill=NA, 
                family="Helvetica")+ 
  theme(strip.background=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks=element_blank(), 
        panel.grid=element_blank(),
        legend.position=c(0.08, 0.90), 
        legend.background = element_blank(), 
        plot.background = element_blank())
p.spectrum

# the Realm barplots
plot.df$value[plot.df$variable=="pre_change"] <- abs(plot.df$value[plot.df$variable=="pre_change"])
# scale to max
tst <- tapply(plot.df$value, plot.df$variable, function(x){x/max(x, na.rm=T)})
plot.df$value_scale <- as.numeric(unlist(tst))
# subset to hotspots
plot.hp <- plot.df[!is.na(plot.df$hotspot_type) ,]
## check max of accumulated threat values for each country to scale axes later
sort(tapply(plot.df$value_scale, plot.df$LEVEL_NAME, sum)) # Hello South Carolina!

# add Realm 
realm.l <- list(Neotropical = c("Mexico Gulf", "Guatemala"),
              Afrotropical = c("Cape Provinces"),
              Indomalayan = c("India", "China South-Central", "East Himalaya", "Philippines",
                              "Thailand","China Southeast","Myanmar","Vietnam","Borneo",
                              "Sumatera","Malaya"),
              Australasia = c("Western Australia", "Queensland"))
realm <- data.frame(realm = c(rep("Neotropical", 2), rep("Afrotropical", 1), rep("Australasia", 2),
                        rep("Indomalayan", 11)),
              LEVEL_NAME = c("Mexico Gulf", "Guatemala", "Cape Provinces",
                          "Western Australia", "Queensland", "India", "China South-Central",
                          "East Himalaya", "Philippines", "Thailand","China Southeast",
                          "Myanmar","Vietnam","Borneo","Sumatera","Malaya"))
plot.hp <- merge(plot.hp, realm, all.x=T)
plot.hp$LEVEL_NAME <- factor(plot.hp$LEVEL_NAME, levels=unlist(realm.l))



bat <- scico(palette="hawaii", n=4, begin=0, alpha=1)
# arrange colors (bat[2] = defore)
bat2 <- c(bat[3], bat[2], bat[4], bat[1])
p.realm <- ggplot(plot.hp, aes(fill=variable, y=value_scale, x=LEVEL_NAME))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual("Threat", values=bat2, labels=as_labeller(c("deforestation2"="deforestation",
                                                            "hfp_mean"="human foot print", 
                                                          "pre_change"="precipitation change",
                                                            "mat_change"="temperature change")))+
  ylab("Scaled, added threat values")+
  coord_cartesian(ylim=c(NA,2), expand=F)+
  xlab("")+
  geom_hline(yintercept=1.8165854, lty=2)+
  annotate("text", x = 16, y = 1.88, label = "max", size=2.5)+
  # geom_rect(aes(xmin = -.12, xmax = -.03, ymin = .55, ymax = 2.45), fill="grey")+
  # geom_rect(aes(xmin = -.12, xmax = -.03, ymin = 2.55, ymax = 3.45), fill="grey")+
  # geom_rect(aes(xmin = -.12, xmax = -.03, ymin = 3.55, ymax = 14.45), fill="grey")+
  # geom_rect(aes(xmin = -.12, xmax = -.03, ymin = 14.55, ymax = 16.45), fill="grey")+
  theme(panel.border=element_rect(fill=NA, color=NA),         
        legend.position="right", 
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=40, hjust=1, size=8))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))
p.realm

plot_grid(p.spectrum, p.realm, ncol=1, labels="AUTO", label_fontface="plain",
          label_size=11, rel_heights=c(1,.5))
ggsave("figures/threats_per_realm.pdf", width=17, height=20, units = "cm", dpi = 300, bg="white")





-----------------------------------------------# (spectrum_PD <- ggplot(plot.df, aes(x=value, y=..scaled..))+
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





# ## more facets - currently not working
# # library(tidytext
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
# ggsave("figures/hotspots_on_the_spectrum.pdf", width=8, height=7, units = "in", dpi = 300, bg="white")
# ggsave("figures/hotspots_on_the_spectrum.png", width=8, height=7, units = "in", dpi = 300, bg="white")





# Hotspots VS non-hotspots - threats ---------------------------

#bc <- c("grey", "#52548D", "#C57391", "#EFB984")
bat <- c(bat[2], bat[3],bat[4],bat[1])
shp2$hotspot_type[is.na(shp2$hotspot_type)] <- "no hotspot"
shp2.tmp <- st_drop_geometry(shp2[,c("hotspot_type", "deforestation2", "hfp_mean", "mat_change", "pre_change")])
shp2.tmp[,-1] <- apply(shp2.tmp[,-1], 2, normalized)

shp2.tmp <- data.table::melt(shp2.tmp)

ggplot(shp2.tmp, aes(x=hotspot_type, y=value, fill=variable))+
  geom_boxplot(varwidth=F)+
  scale_fill_manual("Threat type", values=bat, labels=c("deforestation", "HFP", "MAT change", "PRE change"))+
  scale_y_continuous("Scaled threat values", lim=c(0,1.03))+
  xlab("")+
  theme(legend.position=c(.8,.89), legend.background=element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE, keyheight=unit(6,"mm")))+
  geom_text(aes(x=1.2,y=1.02), label="*", size=10)+
#  geom_text(aes(x=1.7,y=1.01), label="*", size=10, col="red")+
  geom_segment(aes(x=0.7,y=1.015, xend=1.7,yend=1.015), lwd=0.1)
ggsave("figures/hotspots_boxplots.png", width=5, height=4, units = "in", dpi = 300, bg="white")

kruskal.test(shp2$hfp_mean, shp2$hotspot_type)
kruskal.test(shp2$mat_change, shp2$hotspot_type)
kruskal.test(shp2$pre_change, shp2$hotspot_type)
kruskal.test(shp2$deforestation2, shp2$hotspot_type)
pairwise.wilcox.test(shp2$deforestation2, shp2$hotspot_type, p.adjust.method="fdr")

# quantify
tapply(shp2$deforestation2, shp2$hotspot_type, psych::describe)





# Hotspots VS equal-area non-HS ----
# super hotspots + guatemala + Mexico Gulf

# procedure:
# 1. select 5 countries that are most similar in area to some selected hotspots.
# 1.1. the selected country would ideally be located within a Myers hotspot
# 2. check PD, SR, sesPD ratios.
# 3. check biome composition

# search for equal area pairs
a1 <- shp2$area[shp2$hotspot_type %in% c("PD & PE") | shp2$LEVEL_NAME %in% c("Guatemala", "Mexico Gulf")]
a2 <- list(NULL)
for(i in 1:length(a1)){
  tm <- order(abs(shp2$area - a1[i]), decreasing=F)
  pairorder <- shp2$LEVEL_NAME[tm]
  tm <- shp2$LEVEL3_COD[sort(abs(shp2$area - a1[i]))]
  a2[[i]] <- pairorder[1:15]
}

# prep data frame with values to compare
centroids <- st_centroid(shp2) %>% 
  st_coordinates()
shp2$lng <- centroids[,1] # x=lng
shp2$lat <- centroids[,2]
names(a2) <- shp2$LEVEL_NAME[shp2$hotspot_type %in% c("PD & PE") | shp2$LEVEL_NAME %in% c("Guatemala", "Mexico Gulf")]
a2.df <- unlist(a2)
a2.df <- reshape2::melt(a2.df, value.name="twin")
a2.df$position <- as.numeric(regmatches(row.names(a2.df), regexpr("[0-9]{1,2}", row.names(a2.df))))
a2.df$hotspot <- gsub("[0-9]{1,2}", "", row.names(a2.df))
a2.df <- merge(a2.df, 
               st_drop_geometry(shp2[,grep("NAME|area$|richness$|PD_obs$|SES\\.PD$|WE|SES\\.PE|mrd|lat|trop|temp|mont|desert|coverage|soil|mat_mean|pre_mean|tra_mean|prs_mean", names(shp2))]), 
               by.x="twin", by.y="LEVEL_NAME", all.x=T) # 
#View(a2.df)#[a2.df$position %in% c(1,2),])
# show area ranges in groups
ggplot(a2.df, aes(x=hotspot, y=area))+
  geom_boxplot()+
  scale_y_log10()+
  geom_point(data=a2.df[a2.df$position==1,], aes(x=hotspot, y=area), col="red")
# show lat ranges in groups
ggplot(a2.df, aes(x=hotspot, y=abs(lat/100000), label=twin))+
  geom_boxplot()+
  geom_text(size=3)+
  geom_point(data=a2.df[a2.df$position==1,], aes(x=hotspot, y=abs(lat/100000)), col="red")

# Pairs chosen according to similar area, latitude and biome affiliation:
a2.sub <- a2.df
a2.sub$pair <- paste(a2.sub$hotspot, "-", a2.sub$twin)

a2.sub <- a2.sub[grep("China Southeast - Brazil Northeast|China Southeast - China Southeast|Guatemala - Benin|Guatemala - Guatemala|Mexico Gulf - Sierra Leone|Mexico Gulf - Mexico Gulf|Queensland - Saudi Arabia|Queensland - Queensland|Thailand - Brazil South|Thailand - Thailand", a2.sub$pair),]

a2.sub$position[which(a2.sub$position!=1)] <- 2



# extract position 1+2 + hotspot identical, turn into longformat for plotting
et <- a2.sub[a2.sub$position %in% c(1,2),-grep("geometry|pair", names(a2.sub))]
# abs latitude 
et$lat <- abs(et$lat/100000)
#et$type <- c("hotspot", "twin")[et$position]
#format: hotspot_group | position | variables
pdf <- reshape::melt(et, id.vars=c("hotspot", "position", "twin"))
pdf$pair <- paste0(pdf$hotspot, "\n", pdf$twin)
#View(pdf)
ggplot(pdf[-grep("temp_cf|mont_gs|area|lat", pdf$variable),], 
       aes(x=factor(position), y=value, group=hotspot, col=hotspot))+
  geom_point(alpha=.5, size=3)+
  geom_line(alpha=.5, size=3)+
  scale_color_scico_d(labels=sort(unique(pdf$pair))[c(1,3,6,8,9)], end=.9)+
  xlab("1=hotspot, 2=similar non-hotspot")+
  facet_wrap("variable", scales="free_y", ncol=4)+
  theme(strip.background=element_blank(),
        legend.position=c(.6,.1), 
        legend.direction="horizontal")
ggsave("figures/hotspot_pairs.png", width=8, height=6, unit="in", bg="white", dpi=300)








# SI PLOTS #####
load("workspace_maps.RData")




## island + continental comparison equal areas ----
# . Does it make sense to treat
# islands as special case? what to compare: SR, sesPD, etc
# 
# The idea: comparing absolute values is difficult due to unknown area dynamics.
# Here I will compare the difference in equal-area botanical country pairs for
# several parameters

fams <- readRDS("families_per_level.rds")
tmp <- data.frame(LEVEL3_COD=names(fams), fams = lengths(fams))
shp2 <- merge(shp2, tmp, all.x=T)


# build pairs
shp2$area_rank <- rank(shp2$area) # rank areas
tmp <- st_drop_geometry(shp2[order(shp2$area_rank),])
pairs <- vector(length=nrow(tmp))
for(i in seq(1, length(tmp$area_rank), by=2)){
  pairs[c(i,i+1)] <- i
}
tmp$pairs <- pairs
# View(tmp[,c(1:10,ncol(tmp))]) # check it out


## sesPD 
ggplot(tmp, aes(x=area_class, y=sesPD, fill=area_class))+
  geom_boxplot()
tmp1 <- tmp[tmp$area>=min(tmp$area[tmp$area_class=="continental"], na.rm=T) & 
               tmp$area<=max(tmp$area[tmp$area_class=="island"], na.rm=T), ]
kruskal.test(tmp$sesPD, tmp$area_class)
# higher on islands

ggplot(tmp1, aes(x=area_class, y=sesPD, fill=area_class))+
  geom_boxplot()
kruskal.test(tmp1$sesPD, tmp1$area_class)
# just the overlap zone: significant differences between islands and continental

ggplot(tmp, aes(x=area, y=sesPD, col=area_class))+
  geom_point(col="black")+
  geom_point(data=tmp1, size=2)+
  geom_smooth(data=tmp1, aes(group=area_class), method="lm")+
  scale_x_continuous(trans="log")
# in the intermittend area zone, islands have on average a higher sesPD than
# continental regions. Looks like sesPD drops faster in continental units than
# on islands. There are more continental units in the upper end of this overlap
# zone though, this might explain the larger sesPD values between islands and
# continents. Now check area vs sesPD for actual EQUAL AREA PAIRS:

# calc pair stats
tmp2 <- data.frame(pairs=unique(tmp$pairs))
tmp2$richness_diff <- tapply(tmp$richness, tmp$pairs, diff)
tmp2$sesPD_diff <- tapply(tmp$sesPD, tmp$pairs, diff)
tmp2$sesPE_diff <- tapply(tmp$sesPE, tmp$pairs, diff)
tmp2$area_diff <- tapply(tmp$area, tmp$pairs, diff)
tmp2$area_pair_mean <- tapply(tmp$area, tmp$pairs, mean)
tmp2$Family_density_diff <- tapply(tmp$fams_richness, tmp$pairs, diff)

# area type mix
tmp2$area_class_diff <- tapply(tmp$area_class, tmp$pairs, function(x){x[1]!=x[2]})

# visualize pairs
(eaps <- ggplot(tmp, aes(x=area, y=sesPD, col=area_class, group=pairs))+
  geom_point(col="black", pch=1)+
  geom_point(data=tmp1, size=2)+
  scale_color_scico_d("Area type", palette="batlow", begin=0, end=.5, alpha=.7)+
  scale_x_continuous(trans="log")+
  geom_line(data=tmp1, size=.5, col="grey40")+
  theme(legend.position=c(.2,.2)))

ggplot(na.omit(tmp2), aes(x=area_pair_mean, y=abs(sesPD_diff), col=area_class_diff))+
  geom_point()+
  scale_x_continuous(trans="log")
ggplot(na.omit(tmp2), aes(x=area_class_diff, y=abs(sesPD_diff), fill=area_class_diff))+
  geom_boxplot(varwidth=T)
kruskal.test(abs(tmp2$sesPD_diff), tmp2$area_class_diff)

# restrict to area range that has both mixed and non mixed pairs
tmp3 <- tmp2[tmp2$area_pair_mean>=min(tmp2$area_pair_mean[tmp2$area_class_diff==T], na.rm=T) & 
              tmp2$area_pair_mean<=max(tmp2$area_pair_mean[tmp2$area_class_diff==T], na.rm=T), ]
(eaps_ses.pd_box <- ggplot(na.omit(tmp3), aes(x=area_class_diff, y=abs(sesPD_diff), fill=area_class_diff))+
  geom_boxplot(varwidth=T)+
    ylab("Absolute sesPD difference")+
    xlab("")+
    scale_fill_scico_d(palette="lapaz", begin=0.2, end=.8, alpha=.8)+
    theme(legend.position="none"))
kruskal.test(abs(tmp3$sesPD_diff), tmp3$area_class_diff)
# We expect that if there are differences between islands and continental units
# in sesPD which are independent of area, then the differences in mixed equal
# area pairs would be larger than in equal pairs. 
# Comparing the differences in sesPD between equal-area botanical country
# pairs, we find no significant difference between mixed and equal pairs.
# p=0.08). 


## Species richness 
ggplot(na.omit(tmp2), aes(x=area_pair_mean, y=abs(richness_diff), col=area_class_diff))+
  geom_point()+
  geom_hline(yintercept=0, lty=2)+
  scale_x_continuous(trans="log")
ggplot(na.omit(tmp2), aes(x=area_class_diff, y=abs(richness_diff), fill=area_class_diff))+
  geom_boxplot(varwidth=T)
kruskal.test(abs(tmp2$richness_diff), tmp2$area_class_diff)
# restrict to area range that has both mixed and non mixed pairs
(eaps_sr_box <- ggplot(na.omit(tmp3), aes(x=area_class_diff, y=abs(richness_diff), fill=area_class_diff))+
  geom_boxplot(varwidth=T)+
    ylab("Absolute SR difference")+
    xlab("")+
    scale_fill_scico_d("Mixed pair", palette="lapaz", begin=0.2, end=.8, alpha=.8)+
    theme(legend.position=c(.8, .8)))
kruskal.test(abs(tmp3$richness_diff), tmp3$area_class_diff)
# No significant difference in SR between mixed and equal pairs

## sesPE
ggplot(na.omit(tmp2), aes(x=area_pair_mean, y=abs(sesPE_diff), col=area_class_diff))+
  geom_point()+
  geom_hline(yintercept=0, lty=2)+
  scale_x_continuous(trans="log")
ggplot(na.omit(tmp2), aes(x=area_class_diff, y=abs(sesPE_diff), fill=area_class_diff))+
  geom_boxplot(varwidth=T)
kruskal.test(abs(tmp2$sesPE_diff), tmp2$area_class_diff)
# restrict to area range that has both mixed and non mixed pairs
(eaps_ses.pe_box <- ggplot(na.omit(tmp3), aes(x=area_class_diff, y=abs(sesPE_diff), fill=area_class_diff))+
  geom_boxplot(varwidth=T)+
    ylab("Absolute sesPE difference")+
    xlab("")+
    scale_fill_scico_d(palette="lapaz", begin=0.2, end=.8, alpha=.8)+
    theme(legend.position="none"))
kruskal.test(abs(tmp3$sesPE_diff), tmp3$area_class_diff)
# No significant difference in sesPE between mixed and equal pairs

## Family density
ggplot(na.omit(tmp2), aes(x=area_pair_mean, y=abs(Family_density_diff), col=area_class_diff))+
  geom_point()+
  geom_hline(yintercept=0, lty=2)+
  scale_x_continuous(trans="log")
ggplot(na.omit(tmp2), aes(x=area_class_diff, y=abs(Family_density_diff), fill=area_class_diff))+
  geom_boxplot(varwidth=T)
kruskal.test(abs(tmp2$Family_density_diff), tmp2$area_class_diff)
# restrict to area range that has both mixed and non mixed pairs
(eaps_fams_box <- ggplot(na.omit(tmp3), aes(x=area_class_diff, y=abs(Family_density_diff), fill=area_class_diff))+
    geom_boxplot(varwidth=T)+
    ylab("Absolute family density difference")+
    xlab("")+
    scale_fill_scico_d(palette="lapaz", begin=0.2, end=.8, alpha=.8)+
    theme(legend.position="none"))
kruskal.test(abs(tmp3$Family_density_diff), tmp3$area_class_diff)
# No significant difference in sesPE between mixed and equal pairs

# Summary: there are hints that sesPD and sesPE are higher on islands than on
# continental units of equal size, however these are not significant
plot_grid(eaps, eaps_sr_box, eaps_ses.pd_box, eaps_ses.pe_box, 
          ncol=2, label_size=9, label_fontface="plain", labels="AUTO")
ggsave("figures/equal_area_pairs.png", width=7, height=6)

nrow(tmp3)

### number of families per island ----

# plot
min.area <- 6e+9
thicc_lines <- shp2[which(shp2$area<min.area),]
lcol <- min(thicc_lines$fams_richness)/max(fams_richness)
ucol <- 1 #max(thicc_lines$fams_richness)/max(shp$fams_richness)
ggplot(shp2) + 
    geom_sf(aes(fill=fams_richness),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=fams_richness), show.legend=F)+
    scale_colour_scico("", palette = "batlow", trans = "sqrt", 
                       begin = lcol, end = sqrt(ucol))+
    scale_fill_scico("Families/SR", palette = "batlow", trans="sqrt")+ 
  theme_void()+
    theme(legend.position = c(0.22, 0.3),
          legend.key.height = unit(6,"mm"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 10))+
    xlab(" ")
shp2$LEVEL_NAME[order(shp2$fams/shp2$richness, decreasing=T)]

ggplot(na.omit(shp2[,c("area_class", "fams_richness")]), aes(x=area_class, y=fams_richness, fill=area_class))+
  geom_boxplot()+
  ylab("Family density")+
  xlab("")+
  scale_fill_scico_d(palette="lapaz", begin=0.2, end=.8, alpha=.8)+
  theme(legend.position="none")




## plot bot countries sesPE vs sesPD ----
ggplot(shp2, aes(x=sesPD, y=sesPE, label=LEVEL_NAME))+
  geom_point()+
  geom_label(nudge_x=1, hjust=0, label.size=0, label.padding=unit(0.5,"mm"), 
             alpha=0, color = alpha('black', .5), size=3)+
  scale_x_continuous(limits=c(-80,20))
ggsave("figures/SES_PE_vs_SES_PD.png", width=5, height=5)



### kruskal-test for area ----
names(shp2)
shp2$yn <- as.character(!is.na(shp2$hotspot_type)) # TRUE=HS
kruskal.test(shp2$area, shp2$yn)
ggplot(shp2)+
  geom_boxplot(aes(y=area, x=yn))


# a better way to display phylogenetic diversity: 
# get barplot figure, but with phylogenies instead of threats to show where diversity comes from


+# C A C H E 
  
  # ggplot(plot.df, aes(y=value, x=x))+
  #   geom_area(fill="grey85")+
  #   geom_bar(data=plot.df.sub, aes(x=x, y=value, fill=hotspot_type, col=hotspot_type), 
  #            stat="identity", width=.9)+
  #   xlab("rank")+
  #   facet_wrap(~variable, scales="free")+
  # #  scale_color_scico_d(end=0.9)+scale_fill_scico_d(end=0.9)+
  #   scale_color_manual("Hotspot type", values=bc)+scale_fill_manual("Hotspot type", values=bc)+
  #   scale_y_continuous("Value")+
  #   facet_wrap("variable", scales="free", 
#              labeller=as_labeller(c("deforestation2"="deforestation",
#                                     "hfp_mean"="HFP", "pre_change"="future PRE change",
#                                     "mat_change"="future MAT change")))+
#   coord_cartesian(expand=F)+
#   geom_richtext(data=plot.df.sub, aes(x=x, y=value, label=LEVEL_NAME, col=hotspot_type),
#                 size=2.5, angle=90, hjust=-0.05, label.padding=unit(.05,"mm"),
#                 show.legend=F, label.color=NA, family="sans", fill=NA)+
#   theme(strip.background=element_blank(), 
#         axis.text.x=element_blank(), 
#         axis.ticks=element_blank(), 
#         panel.grid=element_blank(),
#         legend.position=c(0.07, 0.85), 
#         legend.background = element_blank(), 
#         plot.background = element_blank())

#ggsave("figures/hotspots_on_the_spectrum_rank_order.pdf", width=17, height=10, units = "cm", dpi = 300, bg="white")


# NEW NEW NEW, version %&%W#¤ECFD¤¤"#"#33232345699!!!111elf


# 
# 
  # # no log on response (crappy model)
  # l1 <- lm(stand_fun(shp2$richness)~log(shp2$area))
  # l2 <- lm(stand_fun(shp2$PD_obs)~log(shp2$area))
  # summary(l1)
  # summary(l2)
  # # --> PD increases faster
  # 
  # # log on response BEFORE standardizing
  # l3 <- lm(stand_fun(log(shp2$richness))~log(shp2$area))
# l4 <- lm(stand_fun(log(shp2$PD_obs))~log(shp2$area))
# summary(l3)
# summary(l4)
# # --> SR increases faster
# 
# # log on response AFTER standardizing - doesnt work, NAs
# # l5 <- lm(log(stand_fun(shp2$richness))~log(shp2$area))
# # l6 <- lm(log(stand_fun(shp2$PD_obs))~log(shp2$area))
# # summary(l5)
# # summary(l6)
# 
# # sswitch to normalize
# l7 <- lm(normalized(shp2$richness)~log(shp2$area))
# l8 <- lm(normalized(shp2$PD_obs)~log(shp2$area))
# summary(l7)
# summary(l8)
# 
# shp2$richness_norm_log <- log(normalized(shp2$richness))
# shp2$PD_obs_norm_log <- log(normalized(shp2$PD_obs))
# 
# tmp <- st_drop_geometry(shp2[,c("richness_norm_log", "PD_obs_norm_log", "area")])
# tmp <- tmp[!is.infinite(rowSums(tmp)),]
# l9 <- lm(data=tmp, richness_norm_log~log(area))
# l10 <- lm(data=tmp, PD_obs_norm_log~log(area))
# summary(l9)
# summary(l10)

# # Ternary plots 
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
#   geom_sf(data=shp2[shp2$sesPE_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
#   geom_sf(data=thicc_lines[thicc_lines$sesPE_spot==1, ], col="#6575B1", show.legend=F, lwd=2) +
#   # PD layer
#   geom_sf(data=shp2[shp2$sesPD_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
#   geom_sf(data=thicc_lines[thicc_lines$sesPD_spot==1, ], col=alpha("#6575B1",0.3), show.legend=F, lwd=2) +
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
# x <- st_drop_geometry(shp2[,c("SR_spot","sesPD_spot", "sesPE_spot")])
# x$SR_spot <- ifelse(x$SR_spot<0, NA,x$SR_spot)
# x$sesPE_spot <- ifelse(x$sesPE_spot<0, NA,x$sesPE_spot)
# x$sesPD_spot <- ifelse(x$sesPD_spot<0, NA,x$sesPD_spot)
# shp2$hotspot_sum <- rowSums(x, na.rm = T)
# 
# ggplot()+
#   # world layer
#   geom_sf(data=shp2, aes(fill=factor(hotspot_sum)),col=NA, lwd=0) +
#   geom_sf(data=thicc_lines[which(thicc_lines$LEVEL3_COD %in% shp2$LEVEL3_COD[shp2$hotspot_sum!=0]), ],
#           col="#6575B1", show.legend=F, lwd=2) +
#   # hotspot sum laye...
#   # geom_sf(aes(fill=hotspot_sum), col=NA, lwd=0, show.legend = F, alpha=al) + 
#   # geom_sf(data=thicc_lines[thicc_lines$sesPE_spot==1, ], col="#6575B1", show.legend=F, lwd=2) +
#   # # PD layer
#   # geom_sf(data=shp2[shp2$sesPD_spot==1, ], fill="#6575B1", col=NA, lwd=0, show.legend = F, alpha=al) + 
#   # geom_sf(data=thicc_lines[thicc_lines$sesPD_spot==1, ], col=alpha("#6575B1",0.3), show.legend=F, lwd=2) +
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
# # Hotspot definition 
# # Myer2000 uses endemism (at least x number of endemics) and habitat loss(at least 70% habitat loss)
# 
# shp2 <- readRDS("including_hotspot_coverage.rds")
# 
# # get a reasonable percent:
# kmean_calc <- function(df, ...){
#   kmeans(df, scaled = ..., nstart = 30)
# }
# km2 <- kmean_calc(shp2$sesPE, 2)
# km3 <- kmean_calc(shp2$sesPE, 3)
# km4 <- kmeans(shp2$sesPE, 4)
# km5 <- kmeans(shp2$sesPE, 5)
# km6 <- kmeans(shp2$sesPE, 6)
# km7 <- kmeans(shp2$sesPE, 7)
# km8 <- kmeans(shp2$sesPE, 8)
# km9 <- kmeans(shp2$sesPE, 9)
# km_perc <- function(x)x$betweenss/x$totss
# plot(unlist(lapply(list(km2,km3,km4,km5,km6,km7,km8,km9),km_perc)))
# plot(shp2$sesPE~shp2$richness, col=km6$cluster) # top 8 look good, matches 2.5% well
# 
# # highest 2.5% of everything:
# ## Endemism hotspots 
# C <- coldspots(shp2$sesPE) # coldspots
# H <- hotspots(shp2$sesPE) # hotspots
# ## Merge endemism values to shapefile of grid cells.
# DF <- data.frame(LEVEL3_COD=shp2$LEVEL3_COD, cold=C, hot=H)
# DF$sesPE_spot <- 0
# DF$sesPE_spot[DF$hot==1] <- 1
# DF$sesPE_spot[DF$cold==1] <- -1
# shp2 <- merge(shp2, DF[,c("LEVEL3_COD", "sesPE_spot")], by = "LEVEL3_COD", all = TRUE)
# shp2$sesPE_spot
# 
# ggplot(shp2, aes(x=factor(sesPE_spot), y=hotspot_coverage, label=LEVEL3_COD))+
#   #geom_point() +
#   geom_text(position="jitter")
# 
# 
# ## PD hotspots 
# C <- coldspots(shp2$sesPD) # coldspots
# H <- hotspots(shp2$sesPD) # hotspots
# DF <- data.frame(LEVEL3_COD=shp2$LEVEL3_COD, cold=C, hot=H)
# DF$sesPD_spot <- 0
# DF$sesPD_spot[DF$hot==1] <- 1
# DF$sesPD_spot[DF$cold==1] <- -1
# shp2 <- merge(shp2, DF[,c("LEVEL3_COD", "sesPD_spot")], by = "LEVEL3_COD", all = TRUE)
# 
# 
# 
# ## Richness hotspots 
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
# ## Deforestation hotspots 
# H <- hotspots(shp2$deforestation2) # hotspots
# shp2$deforest_spot <- H
# 
# # Myers & PD hotspot matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage>0 & shp2$sesPD_spot==1]
# # Myers & PD hotspot no matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage==0 & shp2$sesPD_spot==1]
# 
# # Myers & PE hotspot matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage>0 & shp2$sesPE_spot==1]
# # Myers & PE hotspot no matches:
# shp2$LEVEL_NAME[shp2$hotspot_coverage==0 & shp2$sesPE_spot==1]
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
# ## Maps hotspots 
# min.area <- 1.5e+9
# thicc_lines <- shp2[which(shp2$area<min.area),]
# 
# # extract overlapping parts for another color
# hf <- st_read("hotspots_fixed.gpkg")
# wrld_wrap <- st_wrap_dateline(hf, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
# m <- st_transform(wrld_wrap, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
# # fix sf trouble...
# #sf::sf_use_s2(FALSE)
# #PE.int <- st_intersection(shp2[shp2$sesPE_spot==1,], m)
# # subtract countries from hotspot areas to not cover them
# #hs.int <- st_(m, shp2[shp2$sesPE_spot==1,])
# 
# 
# # put myers map on top
# thicc_lines$sesPE_spot[thicc_lines$sesPE_spot!=0] # no hot/cold zones in small areas
# (PE_hotspot_map <- ggplot(shp2)+
#     geom_sf(aes(fill=factor(sesPE_spot==1)),lwd=.1, col=NA, show.legend = F) + 
#     geom_sf(data=thicc_lines, aes(col=factor(sesPE_spot==1)), show.legend=F, lwd=2)+
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
#     ggtitle("Top sesPE 2.5% with Myers biodiv hotspots"))
# 
# 
# PD_hotspot_map <- ggplot(shp2)+
#   geom_sf(aes(fill=factor(sesPD_spot==1)),lwd=0, col=NA, show.legend = F) + 
#   geom_sf(data=thicc_lines, aes(col=factor(sesPD_spot==1)), show.legend=F, lwd=2)+
#   scale_fill_manual(values = c("grey90", "red"))+
#   scale_color_manual(values = c("grey90", "red"))+
#   theme(panel.border = element_blank())+
#   ggtitle("sesPD Hotspots")+
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
# #breaks=bi_class_breaks(shp2, x=sesPD, y=sesPE, style="jenks", dim=dim)
# shp2$bi_class <- bi_class(shp2, x = richness, y = WE, style = "jenks", dim = dim)$bi_class
# shp2$bi_class_ses <- bi_class(shp2, x = sesPD, y = sesPE, style = "jenks", dim = dim)$bi_class
# thicc_lines <- shp2[which(shp2$area<min.area),]

# (choropl1 <- ggplot() +
#   geom_sf(na.omit(shp2), mapping = aes(fill = bi_class_ses), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = my_pal, dim=dim, na.value="white") +
#   geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class_ses), show.legend=F)+
#   bi_scale_color(pal = my_pal, dim=dim, na.value="white") +
#   bi_theme()+ theme(plot.margin = margin(0.1, -0.2, 0.1, -0.2, "cm"))+
#   coord_sf(expand = F)) 
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="sesPD", ylab="sesPE ", size=8,
#                     breaks=bi_class_breaks(shp2, x=sesPD, y=sesPE, style="jenks", dim=dim))+
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
# ## Threat boxplots for groups 
# # threats: HFP, future clim change, deforestation
# 
# # boxplot for deforest vs biclass groups -a more quantitative plot
# #coropleth_palette <- read.csv("figures/coropleth_palette.txt", header = F)$V1 # actual order
# 
# (biclass_vs_deforest <- ggplot(shp2)+
#   geom_boxplot(aes(x=bi_class_ses, y=deforestation2), varwidth = T,  show.legend=F)+ #, fill=bi_class_ses
#   #scale_fill_manual(values=coropleth_palette)+
#   xlab("sesPD-sesPE group"))
# # it does not look like there is a strong pattern in deforestation. might be single hotspots in it though. gotta inspect visually. check for hfp and future climate change as well:
# (biclass_vs_hfp <- ggplot(shp2)+
#     geom_boxplot(aes(x=bi_class_ses, y=hfp_mean), varwidth = T,  show.legend=F)+ 
#     xlab("sesPD-sesPE group"))
# (biclass_vs_fut_prec <- ggplot(shp2)+
#     geom_boxplot(aes(x=bi_class_ses, y=pre_change), varwidth = T,  show.legend=F)+
#     xlab("sesPD-sesPE group"))
# (biclass_vs_fut_mat <- ggplot(shp2)+
#     geom_boxplot(aes(x=bi_class_ses, y=mat_change), varwidth = T,  show.legend=F)+
#     xlab("sesPD-sesPE group"))
# plot_grid(biclass_vs_deforest, biclass_vs_hfp, biclass_vs_fut_prec, biclass_vs_fut_mat, ncol=2)
# ggsave("figures/boxplot_PD_PE_threats.png", width=7, height=7, units = "in", dpi = 600, bg = "white")
# 
# ## continuous plot
# y <- c("deforestation2", "hfp_mean", "pre_change", "mat_change")
# for(i in 1:4){
#   response <- y[i] 
#   g <- ggplot(shp2, aes_string(x="sesPD", y=response))+
#     geom_point()+
#     geom_smooth()
#   assign(response, g) #generate an object for each plot
# }
# for(i in 1:4){
#   response <- y[i] 
#   g <- ggplot(shp2, aes_string(x="sesPE", y=response))+
#     geom_point()+
#     geom_smooth()
#   assign(paste0(response, "_PE"),g) #generate an object for each plot
# }
# plot_grid(deforestation2, hfp_mean, pre_change, mat_change, 
#           deforestation2_PE, hfp_mean_PE, pre_change_PE, mat_change_PE, ncol=4)
# 
# 
# 
# ## adding Myer hotspots 
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
#                     breaks=bi_class_breaks(shp2, x=sesPD, y=sesPE, style="jenks", dim=dim))+
#   theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_text(angle=45, hjust=1))
# pd_pe_myer <- ggdraw() + draw_plot(choropl3, 0, 0, 1, 1) + draw_plot(legend, 0.05, 0, 0.3, 0.3)
# ggsave("figures/PD_PE_Myer.png", width=10, height=5.74, units = "in", dpi = 600, bg = "white")

# OLD SINGLE CHOROLPLETH PLOTS THREATS 
# # legend position and size
# l.x=0.07; l.y=0.01; l.width=0.25; l.height=0.25
# 
# ## PLOTS 
# PE_deforest <- ggplot() +
#   geom_sf(shp2, mapping = aes(fill = bi_class_PE_deforest), color = NA, size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = my_pal, dim=dim, na.value="white") +
#   geom_sf(data=thicc_lines, lwd=1, aes(col=bi_class_PE_deforest), show.legend=F)+
#   bi_scale_color(pal = my_pal, dim=dim, na.value="white") +
#   bi_theme()+
#   coord_sf(expand = F)
# (legend <- bi_legend(pal=my_pal, dim=dim, xlab="deforestation", ylab="sesPE ", size=8,
#                      breaks=bi_class_breaks(shp2, x=deforestation2, y=sesPE, style="jenks", dim=dim))+
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
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="deforestation", ylab="sesPD ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=deforestation2, y=sesPD, style="jenks", dim=dim))+
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
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="HFP", ylab="sesPE ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=hfp_mean, y=sesPE, style="jenks", dim=dim))+
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
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="HFP", ylab="sesPD ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=hfp_mean, y=sesPD, style="jenks", dim=dim))+
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
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="pre_change", ylab="sesPE ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=pre_change, y=sesPE, style="jenks", dim=dim))+
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
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="pre_change", ylab="sesPD ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=pre_change, y=sesPD, style="jenks", dim=dim))+
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
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="mat_change", ylab="sesPE ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=mat_change, y=sesPE, style="jenks", dim=dim))+
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
# legend <- bi_legend(pal = my_pal, dim=dim, xlab="mat_change", ylab="sesPD ", size=8, 
#                     breaks=bi_class_breaks(shp2, x=mat_change, y=sesPD, style="jenks", dim=dim))+
#   theme(plot.margin=margin(0, 0, 0, 0, "cm"), axis.text.x = element_text(angle=45, hjust=1))
# ggdraw() + draw_plot(matchange_pd, 0, 0, 1, 1) + draw_plot(legend, l.x, l.y, l.width, l.height)
# ggsave("figures/choropleth_PD_matchange.png", width=10, height=5.74, units="in", dpi=300, bg="white")

# 
# ## Ordination of threats 
# dat <- st_drop_geometry(shp2[,c("deforestation2", "hfp_mean", "pre_change", "mat_change")])
# dat <- na.omit(dat)
# tmp <- apply(dat, 2, normalized)
# par(mfrow=c(2,2))
# apply(tmp, 2, hist)
# dat$deforestation2 <- sqrt(dat$deforestation2)
# dat$hfp_mean <- sqrt(dat$hfp_mean)
# tmp <- apply(dat, 2, normalized)
# 
# pca1 <- prcomp(tmp, center = TRUE, scale. = TRUE)
# summary(pca1)
# biplot(pca1, cex=0.9)
# 
# ## Kernel PCA (non-linear pca)
# library(kernlab)
# kpca1 <- kpca(~., data = as.data.frame(tmp), kernel = 'rbfdot')
# pcv(kpca1)
# plot(rotated(kpca1),
#      xlab="1st Principal Component",ylab="2nd Principal Component")


# plot_grid(
# ggplot(shp2, aes(x=hotspot_coverage, y=richness))+
#   geom_density_2d()+
#   geom_smooth()+
#   annotate("text", x=0.9,y=7500, label="rho=ns"),
# ggplot(shp2, aes(x=hotspot_coverage, y=WE))+
#   geom_density_2d_filled(show.legend=F)+
#   scale_fill_scico_d(palette="lajolla")+
#   geom_smooth() + coord_cartesian(ylim=c(0,3500))+
#   annotate("text", x=0.9,y=3000, label="rho=0.17"),
# ggplot(shp2, aes(x=hotspot_coverage, y=sesPD))+
#   geom_rug()+
#   geom_smooth()+
#   annotate("text", x=0.9,y=-10, label="rho=0.26"),
# ggplot(shp2, aes(x=hotspot_coverage, y=sesPE))+
#   geom_point(aes(col=sesPE), show.legend=F)+  scale_color_scico(palette="lajolla")+
#   geom_smooth()+
#   annotate("text", x=0.9,y=25, label="rho=0.34"), 
# nrow=2)

# tmp <- shp2[,grep("LEVEL|hotspot_coverage|sesPE$|sesPD$|WE_s|sdensity", names(shp2))]
# nb <- spdep::poly2nb(tmp, row.names = tmp$LEVEL3_COD)
# col.W <- spdep::nb2listw(nb, style="W", zero.policy = TRUE)
# tmp <- st_drop_geometry(tmp)
# leeHS <- apply(tmp[,grep("sesPD|sesPE|WE_s|sdensity", names(tmp))], 2, spdep::lee.test, y=tmp$hotspot_coverage, 
#                listw=col.W, zero.policy = TRUE, alternative="two.sided", na.action=na.omit)
# leeHS.df <- sapply(leeHS, "[[", "estimate")
# leeHS.df <- rbind(leeHS.df, sapply(leeHS, "[[", "p.value"))
# row.names(leeHS.df) <- c("Lee", "expect", "var", "pvalue")
# leeHS.df <- as.data.frame(t(leeHS.df))
# leeHS.df$twoSD <- 2*sqrt(leeHS.df$var)
# 
# ggplot(leeHS.df, aes(y=expect, x=row.names(leeHS.df)))+
#   geom_linerange(aes(ymin=expect-twoSD, ymax=expect+twoSD), col="grey70")+
#   geom_point(aes(y=Lee, x=row.names(leeHS.df), col=factor(pvalue<0.05)), show.legend = F)+
#   scale_color_manual("p<0.05",values=c("#e4974a"))+
#   ylab("Lee's L with hotspot coverage")+
#   xlab("")+  coord_flip()
# ggsave("figures/LeesL_HS_coverage_and_PD_vars_s.png", width=3, height=2, units = "in", dpi = 300)

