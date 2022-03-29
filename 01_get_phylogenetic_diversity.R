# phylogenetic hotspots


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
library(arrow) # fast parquet file reading
#library(data.table) # fast csv reading
library(castor) # fast tree reading
library(phyloregion) # PD calculations



# dat <- fread("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", quote = "")
# dist <- fread("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt")
# write_parquet(dat, "../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_names.parquet")
# write_parquet(dist, "../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.parquet")


dat <- read_parquet("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_names.parquet")
dist <- read_parquet("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.parquet")



# Get phylogeny -----------------------------------------------------------

allmb <- read_tree(file="../DATA/phylos/ALLMB.tre")








