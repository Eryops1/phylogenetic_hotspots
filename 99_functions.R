
# functions collection

lunique <- function(x){length(unique(x))}

range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}



library(ggcorrplot)


#### plot setup ----
# my_pal1 <- c("#d3d3d3", "#a8aec3", "#7f88b2", "#5665a3", 
#              "#d3bc9a", "#a89b8e", "#7f7a82", "#565a77", 
#              "#d3a45e", "#a88756", "#7f6a4f", "#564e48", 
#              "#d38819", "#a87017", "#7f5816", "#564114")
# my_pal2 <- c("#e8e8e8", "#bddede", "#8ed4d4", "#5ac8c8",
#              "#dabdd4", "#bdbdd4", "#8ebdd4", "#5abdc8",
#              "#cc92c1", "#bd92c1", "#8e92c1", "#5a92c1",
#              "#be64ac", "#bd64ac", "#8e64ac", "#5a64ac")
# np_pal <- c("#e8e8e8", "#b0d5dd", "#76c1d1", "#35abc4", 
#             "#e8e8cc", "#b0d5cc", "#76c1cc", "#35abc4", 
#             "#e8e8a7", "#b0d5a7", "#76c1a7", "#35aba7", 
#             "#e8e840", "#b0d540", "#76c140", "#35ab40")
# np_pal2 <- c("#d3d3d3", "#a0c2c9", "#6bb0be", "#309cb2", "#d4d3ba", "#a1c2b1", "#6cb0a7", "#309c9d", "#d5d498", "#a2c391", "#6cb089", "#319c80", "#d9d53a", "#a4c337", "#6eb134", "#319d31")
# plot(rep(1:4, each=4), rep(1:4, 4), pch=20, cex=5, col=np_pal) # color test
# my_pal <- np_pal2
np_pal <- c("#d3d3d3", "#a0c2c9", "#6bb0be", "#309cb2", "#d4d3ba", "#a1c2b1", "#6cb0a7", "#309c9d", "#d5d498", "#a2c391", "#6cb089", "#319c80", "#d9d53a", "#a4c337", "#6eb134", "#319d31")
my_pal1 <- c("#d3d3d3", "#a0c2c9", "#6bb0be", "#309cb2", "#d4d3ba", "#a1c2b1", "#6cb0a7", "#309c9d", "#d5d498", "#a2c391", "#6cb089", "#319c80", "#d9d53a", "#a4c337", "#6eb134", "#319d31")
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





# function that modifies corrplot from ggcorrplot package
# to allow for visual highlighting of values above a certain thereshold
# adjusting corrplot function to highlight correlations higher lower than threshold
my.corrplot <- function (corr, method = c("square", "circle"), type = c("full", "lower", "upper"), 
                         ggtheme = ggplot2::theme_minimal, title = "", 
                         show.legend = TRUE, legend.title = "Pearson's r", show.diag = FALSE, 
                         colors = c("blue", "white", "red"), outline.color = "gray", 
                         hc.order = FALSE, hc.method = "complete", lab = FALSE, lab_col = "black", 
                         lab_size = 4, p.mat = NULL, sig.level = 0.05, insig = c("pch", "blank"), 
                         pch = 4, pch.col = "black", pch.cex = 5, tl.cex = 12, 
                         tl.col = "black", tl.srt = 45, digits = 2, highlight=TRUE, highthreshold=0.7) 
{
  type <- match.arg(type)
  method <- match.arg(method)
  insig <- match.arg(insig)
  if (inherits(corr, "cor_mat")) {
    cor.mat <- corr
    corr <- .tibble_to_matrix(cor.mat)
    p.mat <- .tibble_to_matrix(attr(cor.mat, "pvalue"))
  }
  if (!is.matrix(corr) & !is.data.frame(corr)) {
    stop("Need a matrix or data frame!")
  }
  corr <- as.matrix(corr)
  corr <- base::round(x = corr, digits = digits)
  if (hc.order) {
    ord <- .hc_cormat_order(corr)
    corr <- corr[ord, ord]
    if (!is.null(p.mat)) {
      p.mat <- p.mat[ord, ord]
      p.mat <- base::round(x = p.mat, digits = digits)
    }
  }
  if (type == "lower") {
    corr <- .get_lower_tri(corr, show.diag)
    p.mat <- .get_lower_tri(p.mat, show.diag)
  }
  else if (type == "upper") {
    corr <- .get_upper_tri(corr, show.diag)
    p.mat <- .get_upper_tri(p.mat, show.diag)
  }
  corr <- reshape2::melt(corr, na.rm = TRUE)
  colnames(corr) <- c("Var1", "Var2", "value")
  corr$pvalue <- rep(NA, nrow(corr))
  corr$signif <- rep(NA, nrow(corr))
  if (!is.null(p.mat)) {
    p.mat <- reshape2::melt(p.mat, na.rm = TRUE)
    corr$coef <- corr$value
    corr$pvalue <- p.mat$value
    corr$signif <- as.numeric(p.mat$value <= sig.level)
    p.mat <- subset(p.mat, p.mat$value > sig.level)
    if (insig == "blank") {
      corr$value <- corr$value * corr$signif
    }
  }
  corr$abs_corr <- abs(corr$value) * 10
  p <- ggplot2::ggplot(data = corr, mapping = ggplot2::aes_string(x = "Var1", 
                                                                  y = "Var2", fill = "value"))
  if (method == "square") {
    p <- p + ggplot2::geom_tile(color = outline.color)
  }
  else if (method == "circle") {
    p <- p + ggplot2::geom_point(color = outline.color, shape = 21, 
                                 ggplot2::aes_string(size = "abs_corr")) + ggplot2::scale_size(range = c(4, 
                                                                                                         10)) + ggplot2::guides(size = FALSE)
  }
  p <- p + ggplot2::scale_fill_gradient2(low = colors[1], high = colors[3], 
                                         mid = colors[2], midpoint = 0, limit = c(-1, 1), space = "Lab", 
                                         name = legend.title)
  if (class(ggtheme)[[1]] == "function") {
    p <- p + ggtheme()
  }
  else if (class(ggtheme)[[1]] == "theme") {
    p <- p + ggtheme
  }
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = tl.srt, 
                                                              vjust = 1, size = tl.cex, hjust = 1), axis.text.y = ggplot2::element_text(size = tl.cex)) + 
    ggplot2::coord_fixed()
  label <- round(x = corr[, "value"], digits = digits)
  if (!is.null(p.mat) & insig == "blank") {
    ns <- corr$pvalue > sig.level
    if (sum(ns) > 0) 
      label[ns] <- " "
  }
  if (lab) {
    if(highlight==TRUE) {
      nh <- abs(corr$value) >= highthreshold
      typo <- c(1, 2)[as.numeric(nh)+1]
      #   p <- p + ggplot2::aes_string(face = typo)
      label <- gsub("0\\.", "\\.", label)
      p <- p + ggplot2::geom_text(
        mapping = ggplot2::aes_string(x = "Var1", y = "Var2", fontface=typo), 
        label = label, color = lab_col, size = lab_size)
    }else
      p <- p + ggplot2::geom_text(
        mapping = ggplot2::aes_string(x = "Var1", y = "Var2"), label = label, color = lab_col, size = lab_size)
  }
  
  if (!is.null(p.mat) & insig == "pch") {
    p <- p + ggplot2::geom_point(data = p.mat, mapping = ggplot2::aes_string(x = "Var1", 
                                                                             y = "Var2"), shape = pch, size = pch.cex, color = pch.col)
  }
  if (title != "") {
    p <- p + ggplot2::ggtitle(title)
  }
  if (!show.legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p <- p + .no_panel()
  p
}
environment(my.corrplot) <- asNamespace('ggcorrplot')


# endemism resampled
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


normalized = function(x){(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))}


# PD_ses with normal distriution check? -------------------------------------------------


PD_ses_normal <- function (x, phy, model = c("tipshuffle", "rowwise", "colwise"),
          reps = 1000, ...)
{
  colnames(x) <- gsub(" ", "_", colnames(x))
  p <- keep.tip(phy, intersect(phy$tip.label, colnames(x)))
  x <- x[, intersect(p$tip.label, colnames(x))]

  PD_obs <- PD(x, p)
  PD_obs <- sqrt(PD_obs)
  
  pd.rand <- switch(model, 
                    tipshuffle = lapply(seq_len(reps), function(i) PD(x, rt(p))), 
                    rowwise = lapply(seq_len(reps), function(i) PD(x[sample(nrow(x)), ], p)), 
                    colwise = lapply(seq_len(reps), function(i) PD(x[, sample(ncol(x))], p)))
  pd_rand <- lapply(pd.rand, sqrt)
  y <- do.call(rbind, pd_rand)
  pd_rand_mean <- apply(X = y, MARGIN = 2, FUN = mean, na.rm = TRUE)
  # transform!
  #pd_rand_mean <- sqrt(pd_rand_mean)
  pd_rand_sd <- apply(X = y, MARGIN = 2, FUN = var, na.rm = TRUE)
  #pd_rand_sd <- sqrt(pd_rand_sd)
  zscore <- (PD_obs - pd_rand_mean)/sqrt(pd_rand_sd)
  pd_obs_rank <- apply(X = rbind(PD_obs, y), MARGIN = 2, FUN = rank)[1,]
  pd_obs_rank <- ifelse(is.na(pd_rand_mean), NA, pd_obs_rank)
  m <- data.frame(grids = rownames(x), 
                  PD_obs, 
                  pd_rand_mean,
                  pd_rand_sd, 
                  pd_obs_rank, 
                  zscore, 
                  pd_obs_p = pd_obs_rank/(reps +1), 
                  reps = reps, row.names = row.names(x))
  z <- data.frame(table(sparse2long(x)$grids))
  names(z) <- c("grids", "richness")
  res <- Reduce(function(x, y) merge(x, y, by = "grids", all = TRUE),
                list(z, m))
  res
}
environment(PD_ses_normal) <- asNamespace('phyloregion')



# 
# bi_wrap <- function(dat, x, y, style, dim){
#   # a wrapper for bi_class and bi_class_breaks
#   varx <- eval(substitute(x))
# #  vary <- substitute(y)
#   
#   classes <- biscale::bi_class(.data=dat, x=varx, y=substitute(y), style=style, dim=dim)$bi_class
#   # bi_breaks <- biscale::bi_class_breaks(.data=dat, x=x, y=y, style=style, dim=dim)
#   # # create numeric vector with case counts
#   # ct <- c(length(which(classes=="1-1")), length(which(classes=="1-2")), length(which(classes=="1-3")), length(which(classes=="1-4")), 
#   #         length(which(classes=="2-1")), length(which(classes=="2-2")), length(which(classes=="2-3")), length(which(classes=="2-4")), 
#   #         length(which(classes=="3-1")), length(which(classes=="3-2")), length(which(classes=="3-3")), length(which(classes=="3-4")), 
#   #         length(which(classes=="4-1")), length(which(classes=="4-2")), length(which(classes=="4-3")), length(which(classes=="4-4")))
#   # ret <- list(classes, bi_breaks, ct)
#   return(ret)
# }


class_col <- function(classes){
  # classes: bi_class vector
  ct <- c(length(which(classes=="1-1")), length(which(classes=="1-2")), length(which(classes=="1-3")), length(which(classes=="1-4")),
          length(which(classes=="2-1")), length(which(classes=="2-2")), length(which(classes=="2-3")), length(which(classes=="2-4")),
          length(which(classes=="3-1")), length(which(classes=="3-2")), length(which(classes=="3-3")), length(which(classes=="3-4")),
          length(which(classes=="4-1")), length(which(classes=="4-2")), length(which(classes=="4-3")), length(which(classes=="4-4")))
  return(ct)
}


# used R packages:
p <- c("biscale","terra","ggcorrplot","ggplot2","ncdf4","raster","exactextractr","data.table","sf",
  "castor","ecospat","phyloregion","treemapify","scico","stringr","gbm","caret","spatialRF")
sort(p)
paste0(p, ", ", collapse=" ")





sr_greedy <- function(species_matrix, m=0, n = nrow(species_matrix)) {
  # rows = area, columns = species. binary.
  # m = missing species tolerated? provide number of species
  # n = max number of botanical countries allowed
  
  # set to store the cells
  cell_set <- logical(nrow(species_matrix))
  names(cell_set) <- row.names(species_matrix)
  # store the number of species in each cell
  species_count <- rowSums(species_matrix)
  # store species count in each iteration
  species <- c()
  area <- c()
  
  # Loop until all species are represented (species count=0)
  while (sum(species_count) > m & sum(cell_set) < n) {
    # Find the cell with the maximum number of species and add it to set
    best_cell <- which.max(species_count)
    cell_set[best_cell] <- TRUE
    # store species count with country
    species <- c(species, max(species_count))
    area <- c(area, names(best_cell))
    # Update the species count for the remaining cells
    ## set species (=columns) that are represented in best cell(=row) to 0
    counted_species <- which(species_matrix[best_cell, ]==1)
    species_matrix[, counted_species] <- 0
    
    species_count <-  rowSums(species_matrix)
  }
  return(list(cell_set, species, area))
}

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


## ------
# #pde_greedy <- function(species_matrix, phylo, m=0, n = nrow(species_matrix)) {
#   # rows = area, columns = species. binary.
#   # m = missing species tolerated? provide number of species
#   # n = max number of botanical countries allowed
#   # phylo = phylo object
#   
#   # set to store the cells
#   cell_set <- logical(nrow(species_matrix))
#   names(cell_set) <- row.names(species_matrix)
#   # store the PD in each cell
#   pe <- phyloregion::phylo_endemism(species_matrix, phylo, weighted=F) # STRICT!
#   # store species count in each iteration
#   pe_number <- c()
#   area <- c()
#   
#   # Loop until all PE is represented (sum(pd)=0)
#   while (sum(pe) > m & sum(cell_set) < n) {
#     # Find the cell with the maximum total PE and add it to set
#     best_cell <- which.max(pe)
#     cell_set[best_cell] <- TRUE
#     # store pd with country
#     pe_number <- c(pe_number, max(pe))
#     area <- c(area, names(best_cell))
#     # Update the pd for the remaining cells
#     ## set species (=columns) that are represented in best cell(=row) to absent
#     counted_species <- which(species_matrix[best_cell, ]==1)
#     species_matrix[, counted_species] <- 0
#     
#     pd <-  phyloregion::phylo_endemism(species_matrix, phylo, weighted=F)
#     message(paste(names(best_cell), " : ", max(pd))) # returns max pd to show progress
#   }
#   return(list(cell_set, pe_number, area))
# }
# 
# #greed <- pde_greedy(submat, phylo=subphy[[1]])
# 
# # too slow, outsourced for parallel on cluster (01c_get_PDEcomplementarity_cluster.R). 
# # read results
# fnames <- dir("PDcomp")[grep("pde_complementarity", dir("PDcomp"))]
# lpd <- list()
# for(i in 1:100){
#   lpd[[i]] <- readRDS(paste0("PDcomp/", fnames[i]))  
# }
# 
# pde_complementarity = sapply(lpd, "[[", 1)
# rowSums(pde_complementarity, na.rm=T) # this is the same everywhere
# res <- data.table(pde_complementarity = pde_complementarity[,1],
#                   LEVEL3_COD = row.names(pde_complementarity))
# 
# pde_added = sapply(lpd, "[[", 2) 
# pde_area = sapply(lpd, "[[", 3)
# apply(pde_area, 1, function(x)length(unique(x)))
# # more than one option for close calls, using most common one though causes some
# # countries to be duplicates which messes up the order
# 
# most_common <- apply(pde_area, 1, function(x)names(sort(table(x),decreasing=TRUE))[1])
# # get average value for each country instead, ignore order
# tmp = data.table(pde_added = as.numeric(pde_added),
#                  LEVEL3_COD = as.character(pde_area))
# average_PDE_added <- tapply(tmp$pde_added, tmp$LEVEL3_COD, mean)
# res2 <- data.table(LEVEL3_COD = names(average_PDE_added),
#                    pde_added = average_PDE_added)
# res <- merge(res, res2, all.x=T, by="LEVEL3_COD")

## -----
