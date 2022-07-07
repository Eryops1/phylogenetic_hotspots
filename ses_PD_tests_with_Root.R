# SES.PD test in and excluding root


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
#library(arrow) # fast parquet file reading
library(data.table) # fast csv reading
library(castor) # fast tree reading
library(phyloregion) # PD calculations
library(raster)
library(sf)
library(ggplot2)
theme_set(theme_bw())
library(cowplot)
library(rgdal)
library(beepr)



# load("PD_nullmodel/comm_and_phy.RData")
# phy <- subphy[[1]]
# rm(subphy)
# gc()
# save(phy, submat, file="tmp.RData")
load("tmp.RData")

tmp <- PD_ses(submat, phy, "tipshuffle", reps=10)

plot(tmp$PD_obs, tmp$pd_rand_mean)
plot(tmp$richness, tmp$zscore)
abline(a=0,b=1)

tmp2 <- PD_ses(submat, phy, "rowwise", reps=10)

plot(tmp2$PD_obs, tmp2$pd_rand_mean, xlim=c(0,150000), ylim=c(0,70000))
abline(a=0,b=1)
plot(tmp2$richness, tmp2$zscore)

plot(tmp$richness, tmp$zscore/tmp$richness)

# diff between tipshuffle + rowwise: 
# “tipshuffle”: shuffles tip labels multiple times.
# “rowwise”: shuffles sites (i.e., varying richness) and keeping species occurrence frequency constant.
# “colwise”: shuffles species occurrence frequency and keeping site richness constant.

# x=matrix, x <- x[, intersect(p$tip.label, colnames(x))] #just makes sure only phylogeny taxa are included, same for phylo 

#p=phylogeny
# rt <- function (phy) {
#   phy$tip.label <- phy$tip.label[sample(length(phy$tip.label))]
#   return(phy)
# } # this shuffles the phylogeny
pd.rand <- switch(model, 
                  tipshuffle = lapply(seq_len(reps), function(i) PD(x, rt(p))), 
                  rowwise = lapply(seq_len(reps), function(i) PD(x[sample(nrow(x)), ], p)), 
                  colwise = lapply(seq_len(reps), function(i) PD(x[, sample(ncol(x))], p)))

test <- submat[sample(nrow(submat)), ]
b <- rowSums(test)
b2 <- rowSums(test)

one <- PD(submat[sample(nrow(submat)), ], phy)
two <- PD(submat[sample(nrow(submat)), ], phy)
three <- PD(submat[sample(nrow(submat)), ], phy)
y <- rbind(as.numeric(one), as.numeric(two))
y <- rbind(y, as.numeric(three))

pd_rand_mean <- apply(X = y, MARGIN = 2, FUN = mean, na.rm = TRUE)
pd_rand_sd <- apply(X = y, MARGIN = 2, FUN = sd, na.rm = TRUE)
pd_rand_mean[1:10]
pd_rand_sd[1:10]
y[,1:10]

# What rowwise shuffle does: Instead of mixing the phylogeny, we are mixing locations: The composition of species stays the same, but we compare the observed PD to PDs that actually occur in other locations for different SRs . This keeps species occurrence frequency (i.e. in how many countries a species occurs) constant. Rare species do not show up in more places than they should. 
# Instead of comparing to made up PDs that do not exist, we compare to real PDs from other locations.
# yeah this does not make sense, produces SR



# Check TUV (Tuvalu) who comes up with random SES.PD (positive, ns):
tuv <- names(which(submat["TUV",]==1))
wcp <- readRDS("../DATA/wcp_apg.rds")
w <- wcp[wcp$plant_name_id %in% tuv,]
length(unique(w$family)) # 30 unique families
tp <- get_subtree_with_tips(phy, only_tips= tuv, collapse_monofurcations = F, force_keep_root = T)$subtree
library(phytools)
plot.phylo(tp, show.node.label = T, type="fan")
sum(tp$edge.length)

# now check cape provinces
cpp <- names(which(submat["CPP",]==1))
w <- wcp[wcp$plant_name_id %in% cpp,]
length(unique(w$family)) # 176 unique families
tp <- get_subtree_with_tips(phy, only_tips= cpp, collapse_monofurcations = F, force_keep_root = T)$subtree
library(phytools)
plot.phylo(tp, type="fan", show.node.label =TRUE)
sum(tp$edge.length)

#compare SR - PD relations between countries:

46/3970 # tuv
15362/88523 # cpp 
# we would expect either less SR or more PD in south africa
plot(tmp$richness, tmp$richness/tmp$PD_obs)
# the ratio gets bigger, meaning PD compared to richness gets smaller - which is to be expected, just not that much. 
points(tmp$richness, tmp$richness/tmp$pd_rand_mean, col="blue")

plot(tmp$richness, tmp$PD_obs)
# the ratio gets bigger, meaning PD compared to richness gets smaller - which is to be expected, just not that much. 
points(tmp$richness, tmp$pd_rand_mean, col="blue")

## A more likely null model:
# How to make the null distribution more likely: We do not expect species to be randomly distributed across the globe, we expect some degree of phylogenetic clustering
# use tipshuffle but only draw null distribution from phylogenetically neighbouring species!!
# -- get MRCA of all taxa in the country, sample from that.

rt2 <- function (phy, x) {
  # subsample, but from the clade and not the entire tree
  presents <- names(which(x[i,]==1))
  nod <- castor::get_mrca_of_set(phy, presents)
  p2 <- extract.clade(phy, nod)
  #p2 <- castor::get_subtree_at_node(phy, nod)
  phy$tip.label <- p2$tip.label[sample(length(phy$tip.label))] # this is wronf, should be length presents
  res <- list(phy, presents, ext.tips=p2$tip.label)
  return(res)
} 

# original=SR, extended= tip label of entire clade
# randomize differently for each country (=row) or for the entire set. for each row! PD uses the entire matrix. need a different phylogeny for each country....???? defies the purpose of null model
ext.tipn <- c()
for(i in 1:nrow(submat)){
  presents <- names(which(submat[i,]==1))
  fams <- unique(wcp$family.apg[wcp$plant_name_id %in% presents])
  specs <- unique(wcp$accepted_plant_name_id[wcp$family.apg %in% fams])
  p2 <- get_subtree_with_tips(phy, specs)$subtree
  #nod <- castor::get_mrca_of_set(phy, presents)
  #p2 <- extract.clade(phy, nod)
  ext.tipn <- c(ext.tipn, length(p2$tip.label))
  cat(i,"\r")
}
plot(ext.tipn) # looks fine for family sampling
plot(tmp$richness, ext.tipn)
hist(ext.tipn)


rt_family <- function (phy, x, wcp) {
  # randomize tip labels, but from the families included only
  presents <- names(which(x[i,]==1))
  fams <- unique(wcp$family.apg[wcp$plant_name_id %in% presents]) # get present families
  specs <- unique(wcp$accepted_plant_name_id[wcp$family.apg %in% fams & wcp$taxon_rank=="Species"]) # all species from present families
  specs <- specs[specs %in% phy$tip.label] # subset all species from families to those present in the phylogeny
  p2 <- get_subtree_with_tips(phy, specs)$subtree # get the tree for present species
  phy$tip.label <- p2$tip.label[sample(length(presents))]
  res <- list(phy, presents, ext.tips=p2$tip.label)
  return(res.int)
} 
environment(rt_family) <- asNamespace('phyloregion')
PD_ses_custom <- 
function (x, phy, wcp, model = c("tipshuffle", "rowwise", "colwise"), reps = 1000, ...) {
  colnames(x) <- gsub(" ", "_", colnames(x))
  p <- keep.tip(phy, intersect(phy$tip.label, colnames(x)))
  x <- x[, intersect(p$tip.label, colnames(x))]
  PD_obs <- PD(x, p)
  pd.rand <- lapply(seq_len(reps), function(i) PD(x, rt_family(p, x, wcp)[[1]])) # check out function above
  y <- do.call(rbind, pd.rand)
  pd_rand_mean <- apply(X = y, MARGIN = 2, FUN = mean, na.rm = TRUE)
  pd_rand_sd <- apply(X = y, MARGIN = 2, FUN = var, na.rm = TRUE)
  zscore <- (PD_obs - pd_rand_mean)/sqrt(pd_rand_sd)
  pd_obs_rank <- apply(X = rbind(PD_obs, y), MARGIN = 2, FUN = rank)[1, ]
  pd_obs_rank <- ifelse(is.na(pd_rand_mean), NA, pd_obs_rank)
  m <- data.frame(grids = rownames(x), PD_obs, pd_rand_mean, 
                  pd_rand_sd, pd_obs_rank, zscore, pd_obs_p = pd_obs_rank/(reps + 
                                                                             1), reps = reps, row.names = row.names(x),
                  org.tips = length(res.int[[2]]))
  z <- data.frame(table(sparse2long(x)$grids))
  names(z) <- c("grids", "richness")
  res <- Reduce(function(x, y) merge(x, y, by = "grids", all = TRUE), 
                list(z, m))
  res
}
environment(PD_ses_custom) <- asNamespace('phyloregion')



tmp3 <- PD_ses_custom(x=submat, phy=phy, wcp=wcp, "tipshuffle", reps=10)

plot(tmp$richness, tmp$zscore)
plot(tmp3$richness, tmp3$zscore)
plot(tmp$zscore, tmp3$zscore)
abline(a=1,b=1)




library(picante)

dat <- as.matrix(submat[1:2,])

pic <- ses.pd(dat, phy, null.model = "independentswap",
       runs = 2, iterations = 2, include.root=TRUE)

# read results from GDK








# whats the average PD per sample size?
plot(tmp$richness, tmp$pd_rand_mean)
points(tmp$richness, tmp$PD_obs, col="blue")


## what if you go beyond max species richness? 
library(castor)
res <- c()
for(i in seq(2,length(phy$tip.label),1000)){
  tp <- get_subtree_with_tips(phy, only_tips= sample(phy$tip, i))$subtree
#  tm <- submat[,which(colnames(submat) %in% tp$tip.label)]
  res <- c(res, sum(tp$edge.length))
  cat(i,"\r")
}
plot(seq(2,length(phy$tip.label),1000), res, ylab="sum branch lengths with increasing SR")
points(tmp$richness, tmp$PD_obs, col="blue")

