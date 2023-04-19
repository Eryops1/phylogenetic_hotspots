
# get rep number (=tree number)
args <- commandArgs()
print(args)
rep <- as.numeric(args[6])
# Get PD complementarity -----------------------------------------------------

pd_greedy <- function(species_matrix, phylo, m=0, n = nrow(species_matrix)) {
  # rows = area, columns = species. binary.
  # m = missing species tolerated? provide number of species
  # n = max number of botanical countries allowed
  # phylo = phylo object
  
  # set to store the cells
  cell_set <- logical(nrow(species_matrix))
  names(cell_set) <- row.names(species_matrix)
  # store the PD in each cell
  pd <- phyloregion::PD(species_matrix, phylo)
  # store species count in each iteration
  pd_number <- c()
  area <- c()
  
  # Loop until all PD is represented (sum(pd)=0)
  while (sum(pd) > m & sum(cell_set) < n) {
    # Find the cell with the maximum total PD and add it to set
    best_cell <- which.max(pd)
    cell_set[best_cell] <- TRUE
    # store pd with country
    pd_number <- c(pd_number, max(pd))
    area <- c(area, names(best_cell))
    # Update the pd for the remaining cells
    ## set species (=columns) that are represented in best cell(=row) to absent
    counted_species <- which(species_matrix[best_cell, ]==1)
    species_matrix[, counted_species] <- 0
    
    pd <-  phyloregion::PD(species_matrix, phylo)
    message(paste(names(best_cell), " : ", max(pd))) # returns max pd to show progress
  }
  return(list(cell_set, pd_number, area))
}

# Load data

load("PD_nullmodel/comm_and_phy.RData")

greed <- pd_greedy(submat, phylo=subphy[[rep]])
saveRDS(greed, paste0("data/pd_complementarity_", rep,".rds"))
