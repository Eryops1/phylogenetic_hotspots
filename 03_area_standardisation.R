# Area standardization --------------------------------------------------------


## CHOSE (N)ONE /// all are explored in more detail in file species_area_standardisation.qmd
# 
# olson <- readRDS("../DATA/PDiv/biomes_olson_ALL.rds")
# shp <- merge(shp, olson[,c("country", "X15", "X16")], by.x="LEVEL3_COD", by.y="country", all.x=TRUE)
# shp$area_ni <- shp$area * (1-shp$X16)
# 
# 
# 
# ## standardize by SAC -----
# library(sars)
# # Uses a dataset in the form of a dataframe with two columns: the first with
# # island/site areas, and the second with the species richness of each
# # island/site.
# 
# fit_mult <- sar_multi(data.frame(area=shp$area_ni[!shp$LEVEL3_COD=="ANT"], shp$richness[!shp$LEVEL3_COD=="ANT"]))
# sort(sapply(fit_mult, "[[", "R2"))
# 
# tmp <- sar_pred(fit_mult$weibull4, shp$area_ni) # gets the theoretical SR
# shp$sdensity <- shp$richness/tmp$Prediction
# 
# ggplot(shp, aes(fill=sdensity))+
#   geom_sf()+
#   scale_fill_scico(palette="batlow")
# 
# 
# plot(shp$area_ni, shp$richness)
# plot(shp$area_ni, shp$sdensity)
# plot(shp$richness, shp$sdensity)
# 
# # same for weighted endemism
# fit_mult <- sar_multi(data.frame(area=shp$area_ni[!shp$LEVEL3_COD=="ANT"], shp$WE[!shp$LEVEL3_COD=="ANT"]))
# sort(sapply(fit_mult, "[[", "R2"))
# tmp <- sar_pred(fit_mult$logistic, shp$area_ni)
# shp$WE_s <- shp$WE/tmp$Prediction
# 
# 
# 
# 
# 
# 
# ## standardize by inhabitable area -------
# ## inhabitable = area - icesheets / glaciated areas
# # 15= lake
# # 16=ice + rock
# 
# shp$sdensity <- shp$richness/(shp$area_ni/1000000) # from sqm -> sqkm
# ggplot(shp, aes(fill=PD_obs, col=PD_obs))+
#   geom_sf()
# View(st_drop_geometry(shp[,c("LEVEL3_COD", "LEVEL_NAME", "area", "area_ni", "sdensity", "richness")]))
# 
# # Selvagens as extreme island with 90 species on 2sqkm
# hist(log(shp$sdensity), breaks=20)
# ggplot(shp, aes(fill=sdensity, col=sdensity))+
#   geom_sf()+
#   scale_fill_scico(trans="log")+
#   scale_color_scico(trans="log")
# 
# shp$WE_s <- shp$WE/(shp$area_ni/1000000) # from sqm -> sqkm
# 
# 
# 
# 
# 
# 
# 
# ## standardize by inhabitable area, min 1500 species -------
# 
# shp$sdensity <- shp$richness/(shp$area_ni/1000000) # from sqm -> sqkm
# shp <- shp[shp$richness>=1000,]
# ggplot(shp, aes(fill=sqrt(sdensity), col=sdensity))+
#   geom_sf()
# View(st_drop_geometry(shp[,c("LEVEL3_COD", "LEVEL_NAME", "area", "area_ni", "sdensity", "richness")]))
# 
# shp$WE_s <- shp$WE/(shp$area_ni/1000000) # from sqm -> sqkm
# 
# 
# 
# 
# ## standardize by SAC > 1500 -----
# shp <- shp[shp$richness>=1000,]
# library(sars)
# # Uses a dataset in the form of a dataframe with two columns: the first with
# # island/site areas, and the second with the species richness of each
# # island/site.
# 
# fit_mult <- sar_multi(data.frame(area=shp$area_ni[!shp$LEVEL3_COD=="ANT"], shp$richness[!shp$LEVEL3_COD=="ANT"]))
# sort(sapply(fit_mult, "[[", "R2"))
# 
# tmp <- sar_pred(fit_mult$weibull4, shp$area_ni) # gets the theoretical SR
# shp$sdensity <- shp$richness/tmp$Prediction
# 
# ggplot(shp, aes(fill=sdensity))+
#   geom_sf(col=NA)+
#   scale_fill_scico(palette="batlow")
# 
# 
# plot(shp$area_ni, shp$richness)
# plot(shp$area_ni, shp$sdensity)
# plot(shp$richness, shp$sdensity)
# 
# # same for weighted endemism
# fit_mult <- sar_multi(data.frame(area=shp$area_ni[!shp$LEVEL3_COD=="ANT"], shp$WE[!shp$LEVEL3_COD=="ANT"]))
# sort(sapply(fit_mult, "[[", "R2"))
# tmp <- sar_pred(fit_mult$logistic, shp$area_ni)
# shp$WE_s <- shp$WE/tmp$Prediction
# 
# 
# 
# 
# ## standardise PD_obs by SR -----
# shp <- shp[!shp$LEVEL3_COD=="ANT",]
# shp$PD_SR <- shp$PD_obs / shp$richness
# hist(log(shp$PD_SR))
# ggplot(shp, aes(fill=PD_SR))+
#   geom_sf(col=NA)+
#   scale_fill_scico(palette="batlow", trans="log")
# 
# 
# 
# # plot rank changes (standardization effects)
# 
# df <- st_drop_geometry(shp)
# df$SR_rank <- rank(df$richness)
# df$SES.PD_rank <- rank(df$SES.PD)
# df$SR_SAC_rank <- rank(df$SR_SAC)
# 
# library(tidyr); library(ggbump); library(dplyr)
# df <- df[,c("LEVEL_NAME", "SR_rank", "SES.PD_rank", "SR_SAC_rank", "area")]
# df2 <- pivot_longer(df, cols = c("SR_rank", "SES.PD_rank", "SR_SAC_rank"))
# 
# alevel = 0.4
# sm <- 0
# df2$name <- factor(df2$name, levels=c("SR_rank", "SES.PD_rank", "SR_SAC_rank"))
# 
# ggplot(df2, aes(x=name, y=value, group=LEVEL_NAME, color=area))+
#   geom_point()+
#   geom_line(alpha=.4, size=1)+
#   coord_cartesian(xlim=c(1.5,2.5))+
#   scale_color_scico("Area", trans="sqrt", alpha = alevel, palette="batlow")


# Islands ------------------------------------------------------------------

is <- vect("../DATA/shapefile_bot_countries/level3_and_islands_merged.gpkg")
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
ggsave("figures/island_yes_no.png", width=4, height=3, units = "in", dpi = 300, bg = "white")


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


# merge into the main shapefile
shp <- merge(shp, df[,c("LEVEL_3_CO", "island", "island_parts", "area_class")], by.x="LEVEL3_COD", by.y="LEVEL_3_CO", all.x=T)

# Dom Rep + Haiti are one big unit and will be treated as continental
shp$area_class[shp$LEVEL3_COD=="DOM"] <- "continental"
shp$area_class[shp$LEVEL3_COD=="HAI"] <- "continental"

# Remove District of Columbia
shp <- shp[-which(shp$LEVEL3_COD=="WDC"),]











## Maps for different species richness standardisations

#### Standardize by inhabitable area

#Divide species richness by area of botanical country - icesheets.

## inhabitable = area - icesheets / glaciated areas
# 15= lake
# 16=ice + rock
shp$SR_area <- shp$richness/(shp$area_ni/1000000) # from sqm -> sqkm
shp$WE_area <- shp$WE/(shp$area_ni/1000000) # from sqm -> sqkm



#### Species-area-curve

#Fit species-area-curve, then divide species richness by theoretical species richness of botanical country.

## standardize by SAC -----
library(sars)
# Uses a dataset in the form of a dataframe with two columns: the first with
# island/site areas, and the second with the species richness of each
# island/site.

fit_mult <- sar_weibull4(data.frame(area=shp$area_ni[!shp$LEVEL3_COD=="ANT"],
                                    shp$richness[!shp$LEVEL3_COD=="ANT"]))
#sort(sapply(fit_mult, "[[", "R2"))

tmp <- sar_pred(fit_mult, shp$area_ni) # gets the theoretical SR
shp$SR_SAC <- shp$richness/tmp$Prediction

# same for weighted endemism
fit_mult <- sar_logistic(data.frame(area=shp$area_ni[!shp$LEVEL3_COD=="ANT"], shp$WE[!shp$LEVEL3_COD=="ANT"]))
tmp <- sar_pred(fit_mult, shp$area_ni)
shp$WE_SAC <- shp$WE/tmp$Prediction


#### Species-area-curve, \>1000 species

#Fit species-area-curve, but for countries with \> 1000 species only. Then divide species richness by theoretical species richness of botanical country.

shp_sub <- shp[shp$richness>=1000,]
library(sars)
# Uses a dataset in the form of a dataframe with two columns: the first with
# island/site areas, and the second with the species richness of each
# island/site.

fit <- sar_weibull4(data.frame(area=shp_sub$area_ni[!shp_sub$LEVEL3_COD=="ANT"], shp_sub$richness[!shp_sub$LEVEL3_COD=="ANT"]))
#sort(sapply(fit_mult, "[[", "R2"))

tmp <- sar_pred(fit, shp_sub$area_ni) # gets the theoretical SR
shp_sub$SR_SAC1000 <- shp_sub$richness/tmp$Prediction

# same for weighted endemism
fit <- sar_logistic(data.frame(area=shp_sub$area_ni[!shp_sub$LEVEL3_COD=="ANT"],
                               shp_sub$WE[!shp_sub$LEVEL3_COD=="ANT"]))
#sort(sapply(fit_mult, "[[", "R2"))
tmp <- sar_pred(fit, shp_sub$area_ni)
shp_sub$WE_SAC1000 <- shp_sub$WE/tmp$Prediction

# merge into shp dataframe
shp <- merge(shp, st_drop_geometry(shp_sub[,c("LEVEL_NAME", "SR_SAC1000", "WE_SAC1000")]), all.x=TRUE)

#### Divide PD by SR

#Alternative to SES.PD: simply divide PD by SR.

## standardise PD_obs by SR -----
#shp <- shp[!shp$LEVEL3_COD=="ANT",]
shp$PD_SR <- shp$PD_obs / shp$richness


### Plot maps

# transform projection
shp <- st_transform(shp, my_projection)

min.area <- 6e+9
shp <- shp[!shp$LEVEL3_COD=="ANT",]
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
# Simple SR 
lcol <- min(thicc_lines$richness)/max(shp$richness)
ucol <- max(thicc_lines$richness)/max(shp$richness)
(sr_map <- ggplot(shp) + 
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


lcol <- min(thicc_lines$SR_area)/max(shp$SR_area)
ucol <- max(thicc_lines$SR_area)/max(shp$SR_area)
transform <- "log"
(SR_area_map <- ggplot(shp) + 
    geom_sf(aes(fill=SR_area),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=SR_area), show.legend=F)+
    scale_colour_scico("", palette="batlow", trans = transform, 
                       begin = lcol, end=ucol)+
    scale_fill_scico("Species/Area", palette="batlow", trans=transform)+ #,breaks = c(0, 0.01, 0.1, 1, 10), labels = c("","0.01","0.1","1", "10")
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

lcol <- min(thicc_lines$SR_SAC)/max(shp$SR_SAC)
ucol <- max(thicc_lines$SR_SAC)/max(shp$SR_SAC)
transform <- "sqrt"
(SR_SAC_map <- ggplot(shp) + 
    geom_sf(aes(fill=SR_SAC),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=SR_SAC), show.legend=F)+
    scale_colour_scico("", palette="batlow", trans = transform, 
                       begin = lcol, end=ucol)+
    scale_fill_scico("Species/SAC", palette="batlow", trans=transform)+ #,breaks = c(0, 0.01, 0.1, 1, 10), labels = c("","0.01","0.1","1", "10")
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

lcol <- min(thicc_lines$SR_SAC1000, na.rm=T)/max(shp$SR_SAC1000, na.rm=T)
ucol <- max(thicc_lines$SR_SAC1000, na.rm=T)/max(shp$SR_SAC1000, na.rm=T)
transform <- "sqrt"
(SR_SAC1000_map <- ggplot(shp) + 
    geom_sf(aes(fill=SR_SAC1000),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=SR_SAC1000), show.legend=F)+
    scale_colour_scico("", palette="batlow", trans = transform, 
                       begin = lcol, end=ucol)+
    scale_fill_scico("Species/SAC1000", palette="batlow", trans=transform)+ #,breaks = c(0, 0.01, 0.1, 1, 10), labels = c("","0.01","0.1","1", "10")
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



lcol <- min(thicc_lines$PD_SR)/max(shp$PD_SR)
ucol <- max(thicc_lines$PD_SR)/max(shp$PD_SR)
transform <- "log"
(PD_SR_map <- ggplot(shp[shp$richness>1000,]) + 
    geom_sf(aes(fill=PD_SR),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines[thicc_lines$richness>1000,], lwd=2, aes(col=PD_SR), show.legend=F)+
    scale_colour_scico("", palette="batlow", trans = transform, 
                       begin = lcol, end=ucol)+
    scale_fill_scico("PD/SR", palette="batlow", trans=transform)+ #,breaks = c(0, 0.01, 0.1, 1, 10), labels = c("","0.01","0.1","1", "10")
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
(nrow(shp[shp$richness>1000,])/100)*5 # 7 bot countries
tmp <- shp[shp$richness>1000,]
tmp$LEVEL_NAME[order(tmp$PD_SR, decreasing=T)][1:14]


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

lcol <- min(thicc_lines$WE_SAC)/max(shp$WE_SAC)
ucol <- max(thicc_lines$WE_SAC)/max(shp$WE_SAC)
(WE_SAC_map <- ggplot(shp2) + 
    geom_sf(aes(fill=WE_SAC),lwd=0, col=NA) + 
    geom_sf(data=thicc_lines, lwd=1, aes(col=WE_SAC), show.legend=F)+
    scale_colour_scico("WE_SAC", palette="batlow", trans = "sqrt", 
                       begin = lcol, end = sqrt(ucol))+
    scale_fill_scico("WE_SAC", palette="batlow", trans="sqrt")+ #, 
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


## Species richness maps

plot_grid(sr_map, SR_area_map, SR_SAC_map, SR_SAC1000_map, ncol = 2, 
          labels=c("A","B","C","D"), label_fontface=1)

## Phylodiversity maps

plot_grid(pd_map, pd_ses_map, PD_SR_map, ncol = 2, 
          labels=c("A","B","C"), label_fontface=1)




## Hotspots
table(bi_class(shp2, x = SR_area, y = SES.PD, style="jenks", dim = 4)$bi_class) # 7 hotspots
table(bi_class(shp2, x = SR_SAC, y = SES.PD, style = "jenks", dim = 4)$bi_class) # 19 hotspots
table(bi_class(shp2, x = SR_SAC1000, y = SES.PD, style = "jenks", dim = 4)$bi_class) # 15 hotspots
table(bi_class(shp2[shp2$richness>1000,], x = richness, y = SES.PD, style = "jenks", dim = 4)$bi_class) # 8 hotspots


# Top 2.5% (= 9 or 7 bot countries) instead:
(nrow(shp2)/100)*2.5
shp2$LEVEL_NAME[order(shp2$SR_area, decreasing=T)][1:9]
shp2$LEVEL_NAME[order(shp2$SR_SAC, decreasing=T)][1:9]
shp2$LEVEL_NAME[order(shp2$SR_SAC1000, decreasing=T)][1:9]

shp2$LEVEL_NAME[order(shp2$PD_obs, decreasing=T)][1:9]
shp2$LEVEL_NAME[order(shp2$SES.PD, decreasing=T)][1:9]
shp2$LEVEL_NAME[order(shp2$PD_SR, decreasing=T)][1:9]
shp2$LEVEL_NAME[order(shp2$SES.PE, decreasing=T)][1:9]


dim <- 4
# regular
shp2$PE_hotspot <- bi_class(shp2, x = WE, y = SES.PE, style = "jenks", dim = dim)$bi_class
shp2$PD_hotspot <- bi_class(shp2, x = richness, y = SES.PD, style = "jenks", dim = dim)$bi_class

# area
shp2$PE_hotspot_area <- bi_class(shp2, x = WE_area, y = SES.PE, style = "jenks", dim = dim)$bi_class
shp2$PD_hotspot_area <- bi_class(shp2, x = SR_area, y = SES.PD, style = "jenks", dim = dim)$bi_class

# SAC
shp2$PE_hotspot_SAC <- bi_class(shp2, x = WE_SAC, y = SES.PE, style = "jenks", dim = dim)$bi_class
shp2$PD_hotspot_SAC <- bi_class(shp2, x = SR_SAC, y = SES.PD, style = "jenks", dim = dim)$bi_class

# SAC >1000
shp2$PE_hotspot_SAC1000 <- bi_class(shp2, x = WE_SAC1000, y = SES.PE, style = "jenks", dim = dim)$bi_class
shp2$PD_hotspot_SAC1000 <- bi_class(shp2, x = SR_SAC1000, y = SES.PD, style = "jenks", dim = dim)$bi_class

thicc_lines <- shp2[which(shp2$area<min.area),]


#Getting the top 2.5% (or 5%) of each variable produces no overlap between a species richness measure and a phylodiversity measure.

df <- shp2[, grepl("PE_hotspot|PD_hotspot|SES\\.PD$|SES\\.PE$|area$", names(shp2))]
df.mlt <- tidyr::pivot_longer(df, cols=grep("hotspot", names(df)), names_to="group")
thicc.mlt <- df.mlt[which(df.mlt$area<min.area),]

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

hs <- ggplot() +
  geom_sf(df.mlt, mapping = aes(fill=value), color = NA, size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal=my_pal, dim=dim, na.value="white") +
  geom_sf(data=thicc.mlt, lwd=1, aes(col=value), show.legend=F)+
  bi_scale_color(pal=my_pal, dim=dim, na.value="white")+
  coord_sf(expand=F)+
  facet_wrap(~group, ncol=2)+
  theme(strip.background=element_blank())

hs

#table(shp2$PD_hotspot)
hsyes <- c("3-3", "3-4", "4-3", "4-4") # this is more inclusive than currently in the manuscript (PE 3-3) 

ggplot() +
  geom_sf(df.mlt, mapping = aes(fill=value%in%hsyes), color = NA, size = 0.1, show.legend = FALSE) +
  #  bi_scale_fill(pal=my_pal, dim=dim, na.value="white") +
  geom_sf(data=thicc.mlt, lwd=1, aes(col=value%in%hsyes), show.legend=F)+
  #  bi_scale_color(pal=my_pal, dim=dim, na.value="white")+
  coord_sf(expand=F)+
  facet_wrap(~group, ncol=2)+
  theme(strip.background=element_blank())

# Clustering = "jenks", 4 groups. Hotspots are in blue.
# 
# -   strictly dividing by area: only island hotspots. This patterns is stable for up to 5000 species (there are enough islands)
# 
# -   dividing by SAR-curve SR yields more balanced results, with small differences (islands)

