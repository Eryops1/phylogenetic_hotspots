
# Prep tree with IDs --------------------------------------------------
# This is using Miaos TACTed phylogeny, but the tip label matching is quick and
# dirty for now, and also using the WCVP 2019, since this is the version he used
# building the tree. Will be changed if the tree is updated to the latest WCVP

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(phytools)
library(data.table)
library(phytools)
library(castor)
library(stringr)

phylo <- read_tree(file="../DATA/phylos/Seed_plant_TOL.full.Sep.22.2021.tre", interpret_quotes = T) # our tree


# match tip labels with WCVP IDs ------------------------------------------

wcp.match <- fread("../DATA/phylos/WCSP_NCBI_merge_Nov182020.ed.txt") # the file with matches
wcp.tips <- paste(wcp.match$order, wcp.match$taxon_name)
wcp.tips <- gsub(" ", "_", wcp.tips)

table(phylo$tip.label %in% wcp.tips) # 234 not matched
phylo$tip.label[!phylo$tip.label %in% wcp.tips] # remove those
phylo <- keep.tip(phylo, phylo$tip.label[phylo$tip.label%in%wcp.tips])

# replace tips labels
new.tips <- wcp.match$accepted_plant_name_id[match(phylo$tip.label, wcp.tips)]
wcp <- readRDS("../BIEN/data/wcp_dec_19.rds")
all(new.tips %in% wcp$accepted_plant_name_id)

phylo$tip.label <- new.tips
phylo$tip.label[1:10]
length(phylo$tip.label)

# drop duplicates
table(duplicated(phylo$tip.label)) # 25 duplicates!
phylo <- drop.tip(phylo, which(duplicated(phylo$tip.label)))
table(duplicated(phylo$tip.label)) # no more!

# drop sub-species levels
table(wcp$taxon_rank, useNA = "ifany")
wcp <- wcp[wcp$taxon_rank=="Species",]

phylo.sp <- keep.tip(phylo, phylo$tip.label[phylo$tip.label %in% wcp$accepted_plant_name_id])
write_tree(phylo.sp, "../DATA/phylos/seed_plant_TOL_WCVP_IDs.tre")


is.binary(phylo.sp)
is.rooted(phylo.sp)
is.ultrametric(phylo.sp)
