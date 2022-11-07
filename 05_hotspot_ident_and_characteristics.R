
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
library(ggtext)
if(!dir.exists("figures"))dir.create("figures")
source("99_functions.R")


load(file="data/workspace_point0.RData")
min.area <- 6e+9
gallpeters_projection <- "+proj=cea +lon_0=0 +x_0=0 +y_0=0 +lat_ts=45 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
my_projection <- gallpeters_projection
shp <- readRDS("data/fin_shp.rds")

# Choropleth PD + SR / PE + WE  ---------------------------

# number of clusters per variable
dim <- 4
shp2$PE_hotspot <- bi_class(shp2, x = WE, y = sesPE, style = "jenks", dim = dim)$bi_class
shp2$PD_hotspot <- bi_class(shp2, x = richness, y = sesPD, style = "jenks", dim = dim)$bi_class
thicc_lines <- shp2[which(shp2$area<min.area),]


df <- shp2[, grepl("PE_hotspot|PD_hotspot|sesPD$|SESPE$|area$", names(shp2))]
df.mlt <- tidyr::pivot_longer(df, cols=grep("hotspot", names(df)), names_to="group")
thicc.mlt <- df.mlt[which(df.mlt$area<min.area),]

my_pal <- my_pal1

# the plots
hs <- ggplot() +
  geom_sf(df.mlt, mapping = aes(fill=value), color = NA, size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal=my_pal, dim=dim, na.value="white") +
  geom_sf(data=thicc.mlt, lwd=1, aes(col=value), show.legend=F)+
  bi_scale_color(pal=my_pal, dim=dim, na.value="white")+
  coord_sf(expand=T)+theme_void()+
  facet_wrap(~group, ncol=2)+
  theme(strip.background=element_blank())

# create legends 
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

# draw maps 
ggdraw() + draw_plot(hs, 0, 0, 1, 1) + 
  draw_plot(PD_hotspot, -0.1, .1, .35, .35)+
  draw_plot(PE_hotspot, 0.4, .1, .35, .35)+
  draw_figure_label(position="top.left", "A", size=11)+
  draw_figure_label(position="top", "B", size=11)

ggsave(paste0("figures/choropleth_hotspots.png"),width=10, height=3, units = "in", dpi = 600, bg = "white")




# TABLE 2 -------------------------------------------------

# check
table(shp2$PD_hotspot)
table(shp2$PE_hotspot)

# define
pdspots <- "3-3"
pespots <- c("3-4", "4-3", "4-4")

#assign
shp2$hotspot_type <- NA
shp2$hotspot_type[shp2$PD_hotspot==pdspots] <- "PD"
shp2$hotspot_type[shp2$PE_hotspot%in%pespots] <- "PE"
shp2$hotspot_type[shp2$PE_hotspot%in%pespots&shp2$PD_hotspot==pdspots] <- "PD & PE"


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

# biome column
biomes <- readRDS("data/biomes_olson_ALL.rds")
biome_names <- c("(Sub)Tropical moist broadleaf forests", 
                 "(Sub)Tropical dry broadleaf forests", 
                 "(Sub)Tropical coniferous forests",
                 "Temperate broadleaf and mixed forests",
                 "Temperate coniferous forests", "Boreal forests/taiga",
                 "(Sub)Tropical grasslands, savannas, shrublands",
                 "Temperate grasslands, savannas, shrublands",
                 "Flooded_grasslands and savannas",
                 "Montane grasslands and shrublands", "Tundra",
                 "Mediterranean forests, woodlands, scrub",
                 "Deserts and xeric shrublands", "Mangroves", "Lakes", "Rock and ice")

names(biomes)[grepl("^X", names(biomes))] <- biome_names

# remove biomes that make up < 3% for clarity
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
rowSums(hab.df[,-c(1,2)], na.rm=T)

hab.mlt <- melt(setDT(hab.df))
hab.mlt <- hab.mlt[order(LEVEL_NAME, value),]
hab.mlt$sort <- rep(length(unique(hab.mlt$variable)):1, length(unique(hab.mlt$LEVEL_NAME)))
# remove empty biomes
hab.mlt <- hab.mlt[hab.mlt$value!=0,]

hab.mlt <- droplevels(hab.mlt)
hab.mlt <- as.data.frame(hab.mlt)

my_pal <- c(scico(12, alpha = .7, begin = 0.1, end = .9, direction = 1, palette = "batlow"), "#dbdbdb")
ggplot(hab.mlt, aes(x=value, y=LEVEL_NAME, fill=variable)) +
  geom_bar(stat="identity", col="black", size=0.2)+
  scale_y_discrete(limits=rev)+
  scale_fill_manual("Biome", values=my_pal)+
  guides(fill=guide_legend(nrow=4,byrow=TRUE))+
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"), 
        legend.spacing.y = unit(.5,"mm"))
ggsave("figures/table1_biomes.pdf", width=6.5, height=5, units="in")
ggsave("figures/table1_biomes2.pdf", width=2, height=5, units="in")

# biome stats
tmp = st_drop_geometry(shp2)
setDF(tmp)
tmp <- merge(tmp, biomes, all.x=TRUE, by.x="LEVEL3_COD", by.y="country")
tmp$hotspot = !is.na(shp2$hotspot_type)
tmp = tmp[,c(66:ncol(tmp))]
tmp2 = melt(tmp, value.name = 'coverage', variable.name='biome')
ggplot(tmp2)+
  geom_boxplot(aes(x=biome, y=coverage, fill=hotspot), outlier.alpha = .3)+
  #  scale_y_continuous(trans='sqrt')+
  scale_fill_scico_d("Hotspot", alpha = .5)+
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        plot.margin = margin(c(5,3,1,25, "mm")))
ggsave('figures/biome_proportions.png', width=7, height=4, dpi=300)

# single tests
for(i in 1:length(unique(tmp2$biome))){
  s = tmp2[tmp2$biome==unique(tmp2$biome)[i],]
  print(kruskal.test(s$coverage, s$hotspot))
  print(tapply(s$coverage, s$hotspot, psych::describe))
  
}

# Div rates
shp2$hotspot_type[is.na(shp2$hotspot_type)] <- "no hotspot"
kruskal.test(shp2$mrd, shp2$hotspot_type)
pairwise.wilcox.test(shp2$mrd, shp2$hotspot_type, p.adjust.method="fdr")

bc <- c("#52548D", "#C57391", "#EFB984")
ggplot(shp2)+
  geom_boxplot(aes(x=hotspot_type, y=mrd, fill=hotspot_type), show.legend = F)+
  labs(y="Net diversification rate (mean root distance)",
       x="")+
  scale_fill_manual(values=c("grey", bc))+
  annotate("text", x=1.5, y=120.7, label="P=0.038", size=2)+
  scale_x_discrete(labels=c("no hotspot", expression(PD[std]), expression(PD[std]+PE[std]), expression(PE[std])))+
  geom_segment(aes(x=1,y=120, xend=2,yend=120), lwd=0.1) #, arrow = arrow(length = unit(0.1,"cm"))+
ggsave("figures/hotspots_boxplots_MRD.png", width=3, height=3, units = "in", dpi = 300, bg="white")





# Hotspots on the threat spectrum -----------------------
### reshape dataframe for plotting
plot.df <- st_drop_geometry(shp2[,c("LEVEL_NAME", "deforestation2", "hfp_mean", 
                                    "pre_change", "mat_change", "hotspot_type")])
plot.df <- reshape::melt(plot.df)
ypos <- rep(rev(seq(0.2,0.80,0.04)), 4)

plot.df <- plot.df[order(plot.df$variable, plot.df$value),]
plot.df$x <- rep(1:(nrow(plot.df)/lunique(plot.df$variable)), lunique(plot.df$variable))
plot.df$lab <- plot.df$LEVEL_NAME
#plot.df$lab[which(plot.df$PD_hotspot!="y"|plot.df$PE_hotspot!="y")] <- NA
plot.df$hotspot_type[plot.df$hotspot_type=="no hotspot"] <- NA
plot.df.sub <- plot.df[which(!is.na(plot.df$hotspot_type)),]

plot.df$value[plot.df$variable=="pre_change" & plot.df$value <= (-2000) & !is.na(plot.df$value)] <- (-2000)
# cosmetics!!!
theme_set(theme_bw()+theme(text=element_text(size=8, family="Helvetica"), 
                           panel.grid=element_blank()))

bc <- c("#52548D", "#C57391", "#EFB984")
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









# Hotspots VS non-hotspots - threats ---------------------------

#bc <- c("grey", "#52548D", "#C57391", "#EFB984")
shp2$hotspot_type[is.na(shp2$hotspot_type)] <- "no hotspot"
shp2.tmp <- st_drop_geometry(shp2[,c("hotspot_type", "deforestation2", "hfp_mean", "mat_change", "pre_change")])
shp2.tmp[,-1] <- apply(shp2.tmp[,-1], 2, normalized)

shp2.tmp <- data.table::melt(shp2.tmp)

ggplot(shp2.tmp, aes(x=hotspot_type, y=value, fill=variable))+
  geom_boxplot(varwidth=F)+
  scale_fill_manual("Threat type", values=bat2, labels=c("deforestation", "HFP", "MAT change", "PRE change"))+
  scale_y_continuous("Scaled threat values", lim=c(0,1.03))+
  xlab("")+
  theme(legend.position=c(.8,.89), legend.background=element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE, keyheight=unit(6,"mm")))+
  annotate("text", x=1.2,y=1.025, label="p=0.016", size=2)+
  geom_segment(aes(x=0.7,y=1.015, xend=1.7,yend=1.015), lwd=0.1)+
  scale_x_discrete(labels=c("no hotspot", expression(PD[std]), expression(PD[std]+PE[std]), expression(PE[std])))
ggsave("figures/hotspots_boxplots.png", width=5, height=4, units = "in", dpi = 300, bg="white")

kruskal.test(shp2$hfp_mean, shp2$hotspot_type)
kruskal.test(shp2$mat_change, shp2$hotspot_type)
kruskal.test(shp2$pre_change, shp2$hotspot_type)
kruskal.test(shp2$deforestation2, shp2$hotspot_type)
pairwise.wilcox.test(shp2$deforestation2, shp2$hotspot_type, p.adjust.method="fdr")

# quantify
tapply(shp2$deforestation2, shp2$hotspot_type, psych::describe)


save(list=c("shp2", "pdspots", "pespots"), file="data/workspace_point1.RData")





