
# get rep number (=tree number)
args <- commandArgs()
print(args)
rep <- as.numeric(args[6])
model <- args[7]

library(phyloregion) # PD calculations etc
library(ecospat)
load('comm_and_phy.RData')

PE_ses <-function (x, phy, model = "tipshuffle", 
          reps = 1000, ...) 
{
  colnames(x) <- gsub(" ", "_", colnames(x))
  p <- keep.tip(phy, intersect(phy$tip.label, colnames(x)))
  x <- x[, intersect(p$tip.label, colnames(x))]
  PE_obs <- phylo_endemism(x, p)
  pe.rand <- switch(model, tipshuffle = lapply(seq_len(reps), function(i) phylo_endemism(x, rt(p))), 
                    rowwise = lapply(seq_len(reps), function(i) phylo_endemism(x[sample(nrow(x)), ], p)), 
                    colwise = lapply(seq_len(reps), function(i) phylo_endemism(x[, sample(ncol(x))], p)))
  y <- do.call(rbind, pe.rand)
  pe_rand_mean <- apply(X = y, MARGIN = 2, FUN = mean, na.rm = TRUE)
  pe_rand_sd <- apply(X = y, MARGIN = 2, FUN = var, na.rm = TRUE)
  zscore <- (PE_obs - pe_rand_mean)/sqrt(pe_rand_sd)
  pe_obs_rank <- apply(X = rbind(PE_obs, y), MARGIN = 2, FUN = rank)[1, ]
  pe_obs_rank <- ifelse(is.na(pe_rand_mean), NA, pe_obs_rank)
  m <- data.frame(grids = rownames(x), PE_obs, pe_rand_mean, 
                  pe_rand_sd, pe_obs_rank, zscore, pe_obs_p = pe_obs_rank/(reps + 1), 
                  reps = reps, row.names = row.names(x))
  z <- data.frame(table(sparse2long(x)$grids))
  names(z) <- c("grids", "richness")
  res <- Reduce(function(x, y) merge(x, y, by = "grids", all = TRUE), 
                list(z, m))
  res
}
environment(PE_ses) <- asNamespace('phyloregion')

# PE_ses
ind <- PE_ses(submat[[rep]], subphy[[rep]], reps=1000)
saveRDS(ind, file=paste0('PE_', rep, '.rds'))

print("PE_ses done")
