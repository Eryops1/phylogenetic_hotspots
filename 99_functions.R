# functions collection
library(ggcorrplot)


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

