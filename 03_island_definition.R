# Islands ------------------------------------------------------------------


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
library(data.table) # fast csv reading
library(sf)
library(terra)
library(ggplot2)
theme_set(theme_bw())

is <- vect("data/shapefile_bot_countries/level3_and_islands_merged.gpkg")
g <- as.data.frame(geom(is))
g <- g[-which(g$hole==1),] # remove holes
is$island_parts <- as.numeric(tapply(g$part, g$geom, function(x)length(unique(x))))

# # area measure
ar <- expanse(is, unit="km", transform=TRUE)
is$area <- ar

df <- as.data.frame(is)
ggplot(df, aes(area, group=island, fill=island))+
  geom_density()+
  scale_x_continuous(trans="log")+
  scale_fill_scico_d(palette="batlow", alpha=.5, begin=0, end=.5)

ggplot(df, aes(x=area, fill=island)) +
  geom_histogram(position="dodge")+
  scale_x_continuous(trans="log")


ggplot(df, aes(area))+
  geom_histogram(data=subset(df,island == T), fill = "#808133", alpha = 0.5) +
  geom_histogram(data=subset(df,island == F), fill = "#001959", alpha = 0.5) +
  #  geom_text(aes(x = area, y = 1, label=LEVEL_NAME), angle = 45)+
  scale_x_continuous(trans="log")+
  scale_fill_scico_d(palette="batlow", alpha=0.5, begin=0, end=.5)+
  geom_vline(xintercept=quantile(df$area[df$island==T], c(.05, .5, .95)), col="#808133", alpha=.5)+
  geom_vline(xintercept=quantile(df$area[df$island==F], c(.05, .5, .95)), col="#001959", alpha=.5)
# A meaningful threshold

## islands bigger than smallest non-island 
df$LEVEL_NAME[df$area >= min(df$area[df$island==F]) & df$island==T]
# Islands area sizes, decreasing
df$LEVEL_NAME[df$island==T][order(df$area[df$island==T], decreasing=T)]
# "New Guinea" "Madagascar" "Sumatera" "Japan" "Philippines"  "Sulawesi"  "New
# Zealand South" "Jawa" "New Zealand North" "Newfoundland" "Cuba"  "Lesser Sunda
# Is." "Iceland" "Dominican Republic_Haiti" "Sakhalin"  "Tasmania" "Sri Lanka"
km2 <- kmeans(df$area[df$island==T], 3)
df$cluster <- NA
df$cluster[df$island==T] <- km2$cluster

ggplot(df, aes(x=area, fill=factor(cluster))) +
  geom_histogram(position="dodge")+
  scale_x_continuous(trans="log")+
  theme(legend.position=c(.2,.7))
#ggsave("figures/island_yes_no.png", width=4, height=3, units = "in", dpi = 300, bg = "white")


## non-islands smaller than biggest island
df$LEVEL_NAME[df$area <= max(df$area[df$island==T]) & df$island==F]

# Non-islands area sizes, increasing
df$LEVEL_NAME[df$island==F][order(df$area[df$island==F])]
#"District of Columbia" "Rhode I." "Delaware" "Cabinda" "Gambia, The"
#"Connecticut"  "Kuwait" "Swaziland"

df$area_class <- "continental"
# Sssssmallest cluster
big_clust <- na.omit(df$cluster[order(df$area, decreasing=F)])[1]
df$area_class[df$cluster == big_clust] <- "island"


# # merge into the main shapefile
shp <- readRDS("data/fin_shp.rds")
shp <- merge(shp, df[,c("LEVEL_3_CO", "island", "island_parts", "area_class")],
             by.x="LEVEL3_COD", by.y="LEVEL_3_CO", all.x=T)

# Dom Rep + Haiti are one big unit and will be treated as continental
shp$area_class[shp$LEVEL3_COD=="DOM"] <- "continental"
shp$area_class[shp$LEVEL3_COD=="HAI"] <- "continental"

# # Remove District of Columbia
# shp <- shp[-which(shp$LEVEL3_COD=="WDC"),]



## island + continental comparison equal areas ----

# fams <- readRDS("data/families_per_level.rds")
# tmp <- data.frame(LEVEL3_COD=names(fams), fams = lengths(fams))
# shp2 <- merge(shp, tmp, all.x=T)


# build pairs
names(shp)[grep("SES\\.PD", names(shp))] <- "sesPD"
names(shp)[grep("SES\\.PE", names(shp))] <- "sesPE"
shp$area_rank <- rank(shp$area) # rank areas
tmp <- st_drop_geometry(shp[order(shp$area_rank),])
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
#tmp2$Family_density_diff <- tapply(tmp$fams_richness, tmp$pairs, diff)

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

# ## Family density
# ggplot(na.omit(tmp2), aes(x=area_pair_mean, y=abs(Family_density_diff), col=area_class_diff))+
#   geom_point()+
#   geom_hline(yintercept=0, lty=2)+
#   scale_x_continuous(trans="log")
# ggplot(na.omit(tmp2), aes(x=area_class_diff, y=abs(Family_density_diff), fill=area_class_diff))+
#   geom_boxplot(varwidth=T)
# kruskal.test(abs(tmp2$Family_density_diff), tmp2$area_class_diff)
# # restrict to area range that has both mixed and non mixed pairs
# (eaps_fams_box <- ggplot(na.omit(tmp3), aes(x=area_class_diff, y=abs(Family_density_diff), fill=area_class_diff))+
#     geom_boxplot(varwidth=T)+
#     ylab("Absolute family density difference")+
#     xlab("")+
#     scale_fill_scico_d(palette="lapaz", begin=0.2, end=.8, alpha=.8)+
#     theme(legend.position="none"))
# kruskal.test(abs(tmp3$Family_density_diff), tmp3$area_class_diff)
# # No significant difference in sesPE between mixed and equal pairs

# Summary: there are hints that sesPD and sesPE are higher on islands than on
# continental units of equal size, however these are not significant
plot_grid(eaps, eaps_sr_box, eaps_ses.pd_box, eaps_ses.pe_box, 
          ncol=2, label_size=9, label_fontface="plain", labels="AUTO")
ggsave("figures/equal_area_pairs.png", width=7, height=6)

nrow(tmp3)

