# phylogenetic hotspots


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
library(arrow)
library(data.table)

# dat <- fread("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", quote = "")
# dist <- fread("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt")
# write_parquet(dat, "../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_names.parquet")
# write_parquet(dist, "../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.parquet")


dat <- read_parquet("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_names.parquet")
dist <- read_parquet("../DATA/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.parquet")


