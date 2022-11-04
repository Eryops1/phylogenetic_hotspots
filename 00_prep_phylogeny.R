
# Prep tree with World Checklist Vascular Plant IDs ------------------------


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(phytools)
library(data.table)
library(phytools)
library(castor)
library(stringr)



# read phylogeny ---------------------------------------------------------
phylo <- read_tree(file="data/GBMB.tre") # genbank tree
# clean labels
phylo$tip.label <- gsub("_", " ", phylo$tip.label)

# all tip label names are from Genbank (NCBI); to get more taxanomy information
# for matching names with the WCVP, download ncbi taxonomy taxdump files here:
# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/



# get ID + authority from taxonomy.dmp --------------------------

text <- readChar("data/names.dmp", 
                 nchars = file.info( "data/names.dmp" )$size, 
                 useBytes = TRUE)
#replace column and row separator
text <- gsub( "\t\\|\t", "|", text)
text <- gsub( "\t\\|\n", "\n", text)
authority <- data.table::fread( text, 
                                header = FALSE, 
                                sep = "|",
                                nThread = 2, 
                                quote="",
                                strip.white = TRUE,
                                fill = TRUE,
                                col.names = c("tax_id", "name_txt", "unique_name", "name_class"))
rm(text)
dat <- data.frame(tips=phylo$tip.label,
                  tax_id=authority$tax_id[match(phylo$tip.label, authority$name_txt)])
long_name <- authority[authority$name_class=="authority",]
dat$name_full <- long_name$name_txt[match(dat$tax_id, long_name$tax_id)]
# if name_full is missing, then there is no author information

dat$tips[!(dat$tips %in% authority$name_txt)] # 98 missing: no match in NCBI

dat$name_clean <- gsub(" \\(.*", "", dat$name_full)
dat$name_clean <- gsub(" [A-Z].*?$", "", dat$name_clean)
dat$author <- NA
for(i in 1:nrow(dat)){
  dat$author[i] <- gsub(dat$name_clean[i], "", dat$name_full[i])
  if(!i%%10)cat(i,"\r")
}
# remove any spaces
dat$author <- str_trim(dat$author)


# get higher tax ranks from taxonomy.dmp ----------------------------------

text <- readChar("data/rankedlineage.dmp", 
                 nchars = file.info( "data/rankedlineage.dmp" )$size, 
                 useBytes = TRUE)
text <- gsub( "\t\\|\t", "|", text)
text <- gsub( "\t\\|\n", "\n", text)
ncbi <- data.table::fread( text, 
                             header = FALSE, 
                             sep = "|",
                             nThread = 2, 
                             quote="",
                             strip.white = TRUE,
                             fill = TRUE,
                             col.names = c("tax_id", "tax_name", "species", "genus", "family", "order",
                                           "class", "phylum", "kingdom", "superkingdom"))
rm(text)
ncbi


# merge higher ranks into dat
dat <- merge(dat, ncbi, by="tax_id", all.x=TRUE)
dat <- dat[-which(is.na(dat$tax_id)),]
fwrite(dat, "data/ncbi.csv")




# create input common format for taxonomy merger ----------------------------

## get authors from authority data frame, rest from lineage data frame
dat <- fread("data/ncbi.csv")
table(phylo$tip.label %in% dat$tips) # 98 who could not be matched

# create unique ID (some tax ids are >1 due to synonyms)
(dups <- which(duplicated(dat$tax_id)))
tmp <- dat[dat$tax_id %in% dat$tax_id[dups],] # 157 rows
tmp <- tmp[!tmp$tips != tmp$tax_name,] # discarded 104

# no duplicate IDs, please
dat <- dat[!dat$tax_id %in% dat$tax_id[dups],]
dat <- rbind(dat, tmp)

which(duplicated(dat$tax_id)) # zero?
table(phylo$tip.label %in% dat$tips) # 202 missing due to removed duplicate IDs

# extract species name ----------------------------------------------------

dat$species <- gsub("^.*? ", "", dat$tax_name)
dat$species <- gsub(" .*", "", dat$species)



# common format file --------------------------------------------------

split_length <- unlist(lapply(strsplit(as.character(dat$tax_name), split = " "), length))
table(split_length)

ncbi_input <- data.frame(taxonID = dat$tax_id, 
                      scientificName= dat$tax_name,
                      tip_label = dat$tips,
                      family = dat$family,
                      genus = dat$genus,
                      species = dat$species,
                      split_length=split_length,
                      author = dat$author, 
                      genus_hybrid = rep(NA, nrow(dat)),
                      species_hybrid = rep(NA, nrow(dat)),
                      taxon_rank = rep(NA, nrow(dat)),
                      infra_name = rep(NA, nrow(dat)),
                      comment = rep(NA, nrow(dat)),
                      usable = rep(NA, nrow(dat)))

# order the dataframe by split length
ncbi_input <- ncbi_input[order(ncbi_input$split_length),]
ncbi_input$id <- c(1:nrow(ncbi_input))

split_list <- strsplit(as.character(ncbi_input$scientificName), split = " ")
names(split_list) <- ncbi_input$id

# Vectorize + loop mix

## split length == 2
ind <- which(ncbi_input$split_length==2)
ncbi_input$taxon_rank[ind] <- "species"

## split length == 3
ind <- which(ncbi_input$split_length==3)
for(i in 1:length(ind)){
  if(split_list[[ind[i]]][1]=="x"){
    ncbi_input$genus_hybrid[ind[i]] <- "x"
    ncbi_input$genus[ind[i]] <- split_list[[ind[i]]][2]
    ncbi_input$species[ind[i]] <- split_list[[ind[i]]][3]
    ncbi_input$taxon_rank[ind[i]] <- "species"
  }
  if(split_list[[ind[i]]][2]=="x"){ # does not occur but you never know
    if(grepl("[A-Z]",split_list[[ind[i]]][3])){
      ncbi_input$genus_hybrid[ind[i]] <- "x"
      ncbi_input$usable[ind[i]] <- "no"
    }else{
      ncbi_input$species_hybrid[ind[i]] <- "x"
      ncbi_input$genus[ind[i]] <- split_list[[ind[i]]][[1]]
      ncbi_input$species[ind[i]] <- split_list[[ind[i]]][[3]]
      ncbi_input$taxon_rank[ind[i]] <- "species"
    }
  }
  if(split_list[[ind[i]]][2]=="hybrid" & split_list[[ind[i]]][3]=="cultivar"){
    ncbi_input$species_hybrid[ind[i]] <- "x"
    ncbi_input$usable[ind[i]] <- "no"
    
  }
  if(split_list[[ind[i]]][1]!="x" & split_list[[ind[i]]][2]!="x" & !grepl("[A-Z]",split_list[[ind[i]]][3])){ 
    ncbi_input$genus[ind[i]] <- split_list[[ind[i]]][[1]]
    ncbi_input$species[ind[i]] <- split_list[[ind[i]]][[2]]
    ncbi_input$taxon_rank[ind[i]] <- NA
    ncbi_input$infra_name[ind[i]] <- split_list[[ind[i]]][[3]]
  }
  if(split_list[[ind[i]]][2]!="x" & grepl("[A-Z]",split_list[[ind[i]]][3])){
    ncbi_input$genus[ind[i]] <- split_list[[ind[i]]][[1]]
    ncbi_input$taxon_rank[ind[i]] <- "genus"
  }
}

## split length == 4
ind <- which(ncbi_input$split_length==4)
for(i in 1:length(ind)){
  if(split_list[[ind[i]]][[1]]=="x" & split_list[[ind[i]]][[3]]=="x"){
    ncbi_input$genus_hybrid[ind[i]] <- "x"
    ncbi_input$species_hybrid[ind[i]] <- "x"
    #ncbi_input$genus[ind[i]] <- split_list[[ind[i]]][[2]]
    ncbi_input$species[ind[i]] <- split_list[[ind[i]]][[4]]
    ncbi_input$taxon_rank[ind[i]] <- "species"
  }
  if(split_list[[ind[i]]][[1]]!="x" & split_list[[ind[i]]][[3]]=="x"){
    ncbi_input$species_hybrid[ind[i]] <- "x"
    #ncbi_input$genus[ind[i]] <- split_list[[ind[i]]][[1]]
    ncbi_input$species[ind[i]] <- split_list[[ind[i]]][[2]]
    #    ncbi_input$infra_name[ind[i]] <- split_list[[ind[i]]][[4]]
    ncbi_input$taxon_rank[ind[i]] <- "species"
  }
  if(split_list[[ind[i]]][[1]]!="x" & split_list[[ind[i]]][[3]]!="x"){
    #ncbi_input$genus[ind[i]] <- split_list[[ind[i]]][[1]]
    ncbi_input$species[ind[i]] <- split_list[[ind[i]]][[2]]
    ncbi_input$taxon_rank[ind[i]] <- split_list[[ind[i]]][[3]]
    ncbi_input$infra_name[ind[i]] <- split_list[[ind[i]]][[4]]
    if(split_list[[ind[i]]][[3]]==""){ncbi_input$taxon_rank[ind[i]] <- NA}
  }
  if(split_list[[ind[i]]][[2]]=="sp."){
    #ncbi_input$species_hybrid[ind[i]] <- "x"
    #ncbi_input$genus[ind[i]] <- split_list[[ind[i]]][[1]]
    ncbi_input$species[ind[i]] <- split_list[[ind[i]]][[2]]
    #    ncbi_input$infra_name[ind[i]] <- split_list[[ind[i]]][[4]]
    ncbi_input$taxon_rank[ind[i]] <- "species"
  }
  if(!i%%100)cat(i,"\r")
}

## split length == 5
ind <- which(ncbi_input$split_length==5)
if(length(ind!=0)){
  for(i in 1:length(ind)){
    if(split_list[[ind[i]]][2]=="x"){
      ncbi_input$genus[ind[i]] <- split_list[[ind[i]]][1]
      ncbi_input$species_hybrid[ind[i]] <- split_list[[ind[i]]][2]
      ncbi_input$species[ind[i]] <- split_list[[ind[i]]][3]
      ncbi_input$taxon_rank[ind[i]] <- split_list[[ind[i]]][4]
      ncbi_input$infra_name[ind[i]] <- split_list[[ind[i]]][5]
    }
    if(split_list[[ind[i]]][4]=="x"){
      ncbi_input$genus[ind[i]] <- split_list[[ind[i]]][1]
      ncbi_input$species_hybrid[ind[i]] <- split_list[[ind[i]]][4]
      ncbi_input$species[ind[i]] <- split_list[[ind[i]]][2]
      ncbi_input$taxon_rank[ind[i]] <- split_list[[ind[i]]][3]
      ncbi_input$infra_name[ind[i]] <- split_list[[ind[i]]][5]
    }
    if(split_list[[ind[i]]][3]=="x"){  # e.g. Dieffenbachia nitidipetiolada x d. oerstedii
      ncbi_input$usable[ind[i]] <- "no"
    }
  }
}

## split length == 6 
ind <- which(ncbi_input$split_length==6)
ncbi_input$usable[ind] <- "no"
ncbi_input$usable[ncbi_input$split_length %in% c(7,8)] <- "no"
ncbi_input$usable[grep("cf\\.|aff\\.", ncbi_input$scientificName)] <- "no"
ncbi_input$usable[ncbi_input$species=="sp."] <- "no"

table(phylo$tip.label %in% ncbi_input$tip_label)

# remove unusable entries?
ncbi_input <- ncbi_input[-which(ncbi_input$usable=="no"),] 
table(phylo$tip.label %in% ncbi_input$tip_label) # 1167 tip labels which will not be matched

# save
saveRDS(ncbi_input, "data/ncbi_common_format.rds")


# check family names APG IV -----------------------------------------------

## NCBI
ncbi_input <- readRDS("ncbi_common_format.rds")
apg <- fread("data/apgweb_parsed.csv")
unique(ncbi_input$family)
ncbi_input$family[!ncbi_input$family %in% apg$Syn_Fam]
# Ripogonaceae is a new family, used to be included in Smilacaceae
# Crambidae should be Cactaceae
ncbi_input$family[ncbi_input$family=="Crambidae"] <- "Cactaceae"
ncbi_input$family <- apg$Acc_Fam[match(ncbi_input$family, apg$Syn_Fam)]
ncbi_input$family.apg <- ncbi_input$family
saveRDS(ncbi_input, "data/ncbi_input.rds")


## WCVP
wcp <- fread("data/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", quote="")
unique(wcp$family)
unique(wcp$family[!wcp$family %in% apg$Syn_Fam])
# Pseudotubulare = not a family name, Synonym entry. Remove
# Gigaspermaceae = moss family
# Incertae_sedis: remove
# Tiganophytaceae: includes monotypic genus, but alright
# Osmundaceae: a fern family
wcp <- wcp[!wcp$family=="Pseudotubulare",]

wcp$family <- apg$Acc_Fam[match(wcp$family, apg$Syn_Fam)]
wcp$family.apg <- wcp$family
saveRDS(wcp, "data/wcp_apg.rds")





# run taxonomy matcher to match with WCVP IDs -----------------------------

source("00_taxonomy_matcher_for_PDiv.R", print.eval = TRUE, chdir = TRUE)




# Replace tip labels with WCVP IDs ----------------------------------------
gc()

ncbi <- readRDS("data/fin_species_match_NCBI_wcvp2022.rds")
table(is.na(ncbi$accepted_plant_name_id))
names(ncbi)

phylo <- read_tree(file="data/GBMB.tre") # genbank tree
phylo$tip.label <- gsub("_", " ", phylo$tip.label)
phylo_org <- phylo
table(phylo$tip.label %in% ncbi$tip_label) 

# FALSE = tip labels that have been marked unusable etc while building the
# common format, eg cf. aff. etc). These tip labels did not enter the taxonomy
# matching procedure

# use species level
ncbi$accepted_plant_name_id[!is.na(ncbi$elevated_to_species_id)] <- ncbi$elevated_to_species_id[!is.na(ncbi$elevated_to_species_id)]

# actual replacement
table(is.na(match(phylo$tip.label, ncbi$tip_label))) # 1167, this is correct
phylo$tip.label[!phylo$tip.label %in% ncbi$tip_label]

# mark NAs in ncbi replacement file with "no matching found"
ncbi$accepted_plant_name_id[which(is.na(ncbi$accepted_plant_name_id))] <- "no match with WCVP found"

phylo$tip.label <- ncbi$accepted_plant_name_id[match(phylo$tip.label, ncbi$tip_label)]
table(is.na(phylo$tip.label)) 
table(phylo$tip.label == "no match with WCVP found") 

phylo <- keep.tip(phylo, phylo$tip.label[!is.na(phylo$tip.label)]) # remove NA tips
phylo <- keep.tip(phylo, phylo$tip.label[phylo$tip.label!="no match with WCVP found"]) # remove unmatched tips

# duplicates?
length(which(table(phylo$tip.label)>1)) 


write.tree(phylo, "data/gbmb_matched_wcvp2022.tre")




# Fix some non monophyletic clades ----------------------------------------
# for TACT to work clades in the tree have to be monophyletic. A single taxon
# outside causes missing species to be TACTed all over the tree

rm(list = setdiff(ls(), lsf.str())) 
tree <- read_tree(file="data/gbmb_matched_wcvp2022.tre")
wcp <- readRDS("data/wcp_apg.rds")
ferns <- as.vector(read.csv("data/fern_list.txt")[,1])
wcp <- wcp %>% 
  dplyr::select(-family) %>%
  dplyr::rename("family" = "family.apg")

# remove ferns
wcp <- wcp[!wcp$family %in% ferns,]

# add orders. Do this here because its specific for the tree (changes in APG system)
apg <- fread("data/apgweb_parsed.csv")
apg <- unique(apg[,c("Acc_Fam", "Clade")])
apg$Acc_Fam[duplicated(apg$Acc_Fam)]
tree$node.label[!tree$node.label==""]

# check which Order names are being used in the tree
apg[apg$Acc_Fam %in% apg$Acc_Fam[duplicated(apg$Acc_Fam)],]
grep("Asparagales", tree$node.label[!tree$node.label==""])
grep("Orchidales", tree$node.label[!tree$node.label==""])
# Orchidales: no
grep("Sabiales", tree$node.label[!tree$node.label==""])
# Sabiales: no

# remove clade names not included in tree from apg file
apg <- apg[-grep("Orchidales|Sabiales", apg$Clade),]
names(apg)[grep("Clade", names(apg))] <- "order"
table(unique(wcp$family) %in% apg$Acc_Fam)

wcp <- merge(wcp, apg, by.x="family", by.y="Acc_Fam", all.x=TRUE)


# get "good species" list to add with TACT -------------------
goodspp <- wcp[wcp$genus_hybrid=="" &
              ""==wcp$species_hybrid &
              ""==wcp$infraspecific_rank &
              ""!=wcp$species &
                 wcp$taxon_status == "Accepted",]
fwrite(goodspp, "data/goodsp.csv")

goodsp <- fread("data/goodsp.csv")
str(goodsp)
goodsp.sub <- goodsp[which(goodsp$accepted_plant_name_id %in% tree$tip.label),]


# check tree attributes (binary, ultrametric etc) ----------------------------
any(tree$edge.length==0)
#tree$edge.length[tree$edge.length==0] <- 1.6e-05
out <- castor::extend_tree_to_height(tree)
out$max_extension
t <- out$tree
is.ultrametric(t, tol=1e-15) # T
is.binary(t)
t <- multi2di(t)


# check monophylies
t0 <-t
goodsp <- fread("data/goodsp.csv")
goodsp.t0 <- goodsp[goodsp$accepted_plant_name_id %in% t0$tip.label,]

orders <- sort(unique(goodsp.t0$order))
families <- sort(unique(goodsp.t0$family))
genera <- sort(unique(goodsp.t0$genus))
res <- data.frame(unique(goodsp.t0[,c("order", "family", "genus")], ord.mono=NA, fam.mono=NA, gen.mono=NA))
res$ord.mono <- NA
res$fam.mono <- NA
res$gen.mono <- NA
for(i in 1:length(orders)){
  temp <- goodsp.t0[goodsp.t0$order %in% orders[i],]
  res[which(res$order==orders[i]), 4] <- is_monophyletic(t0, temp$accepted_plant_name_id)
  if(!i%%1)cat(i,"\r")
}
for(i in 1:length(families)){
  temp <- goodsp.t0[goodsp.t0$family %in% families[i],]
  res[which(res$family==families[i]), 5] <- is_monophyletic(t0, temp$accepted_plant_name_id)
  if(!i%%1)cat(i,"\r")
}
for(i in 1:length(genera)){
  temp <- goodsp.t0[goodsp.t0$genus %in% genera[i],]
  res[which(res$genus==genera[i]), 6] <- is_monophyletic(t0, temp$accepted_plant_name_id)
  if(!i%%1)cat(i,"\r")
}

nm <- res
nm.orders <- unique(nm$order[!nm$ord.mono])

res <- data.frame(nm.orders, all.in=NA, others.in=NA)
too.much <- list()
missing <- list()
for(i in 1:length(nm.orders)){
  wcvp.taxa <- goodsp.sub$accepted_plant_name_id[goodsp.sub$order%in%nm.orders[i]]
  # get the tree order clade
  node.name <- t$node.label[grep(nm.orders[i], t$node.label)]
  t.ord.clade <- extract.clade(t, node.name)

  # all tips in the actual order node?
  res[i, "all.in"] <- all(wcvp.taxa%in%t.ord.clade$tip.label)
  if(!all(wcvp.taxa%in%t.ord.clade$tip.label)){
    missing[[nm.orders[i]]] <- wcvp.taxa[!wcvp.taxa%in%t.ord.clade$tip.label]
  }

  # extra species in the order node?
  res[i, "others.in"] <- any(!t.ord.clade$tip.label%in%wcvp.taxa)
  if(any(!t.ord.clade$tip.label%in%wcvp.taxa)){
    too.much[[nm.orders[i]]] <- t.ord.clade$tip.label[!t.ord.clade$tip.label%in%wcvp.taxa]
  }

  if(!i%%1)cat(i,"\r")
}
res
all(unlist(too.much) %in% unlist(missing))
(rogue <- unique(c(unlist(too.much),unlist(missing))))

# fix, check + save
tp <- drop.tip(t, rogue)
goodsp.tp <- goodsp[goodsp$accepted_plant_name_id%in%tp$tip.label,]

for(i in 1:length(unique(goodsp.tp$order))){
  temp <- goodsp.tp[goodsp.tp$order %in% unique(goodsp.tp$order)[i],]
  print(is_monophyletic(tp, temp$accepted_plant_name_id))
}
# all TRUE



write.tree(tp, "data/gbmb_matched_monophyletic_orders.tre")




# *** Write Goodsp as clean csv -----------------------------------------------


goodsp <- fread("data/goodsp.csv")
# remove those where accepted != plant id
goodsp <- goodsp[!goodsp$accepted_plant_name_id != goodsp$plant_name_id,]
names(goodsp)
mod <- goodsp[,c("order", "family", "genus", "accepted_plant_name_id")]
any(duplicated(mod$accepted_plant_name_id))

# manual adjustments
unique(mod$genus[mod$family==""])
# Isoetes = fern / fern  ally (Isoetaceae)
# Leptopteris = another fern genus
# Osmunda"      "Osmundastrum" "Todea" --> more ferns
mod$family[mod$genus=="Ripogonum"] <- "Rhipogonaceae"
mod$order[mod$genus=="Ripogonum"] <- "Liliales"

mod$family[mod$genus=="Tiganophyton"] <- "Tiganophytaceae"
mod$order[mod$genus=="Tiganophyton"] <- "Brassicales"

mod <- mod[!mod$genus %in% unique(mod$genus[mod$family==""]),]


# all tips in the taxonomy file?
table(tp$tip.label %in% mod$accepted_plant_name_id)
wcp[wcp$plant_name_id %in% tp$tip.label[which(!tp$tip.label %in% mod$accepted_plant_name_id)],]



fwrite(mod, sep=",", file="data/goodsp_wcp_2022_forTACT.csv")


