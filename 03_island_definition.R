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
