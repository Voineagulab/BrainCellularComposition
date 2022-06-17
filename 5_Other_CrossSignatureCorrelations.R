## Setup
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Other/")

## Load signatures
  load("../../../Data/Preprocessed/Signatures - Brain.rda")
  load("../../../Data/Preprocessed/heart.rda"); heart <- heart[-1]
  load("../../../Data/Preprocessed/CrossTissue_Pancreas.rda"); pancreas <- pancreas[-1]
  
## Celltype filtering of signatures
  heart <- lapply(heart, function(x) {
    colnames(x) <- c("CM", "EC", "FB", "SMC")
    return(x)
  })
  
  sigsBrain <- lapply(sigsBrain, function(x) {
    x <- x[,which(colnames(x) %in% colnames(sigsBrain$IP))]
    return(x)
  })
    

## Functions
  # process signature lists
  process <- function(sigList) {
    genes <- lapply(sigList, rownames)
    genes <- do.call("c", genes)
    genes <- table(genes)
    genes <- names(genes)[which(genes == length(sigList))]
    
    # collect gene expression for core celltypes
    exp <- lapply(sigList, function(x) { x[genes,] })
    exp <- do.call("cbind", exp)
    exp <- log2(exp + 0.5)
    return(exp)
  }

  # heatmap
  my.heatmap <- function(exp, ct, key = FALSE, mar = 10) {
    colours <- viridis_pal()(10)
    
    e <- exp[,grep(ct, colnames(exp))]
    colnames(e) <- gsub(paste0("\\.", ct), "", colnames(e))
    
    heatmap.2(cor(e, method = "s"), dendrogram = "row", trace = "none",  density.info = "none", 
              key.title = "NA", key = key, keysize = 2.5, key.xlab = "Spearman", 
              main = NULL, col = colours, breaks = seq(0, 1, 0.1),
              cexRow = 1.2, cexCol = 1.2, lwid = c(2,mar), lhei = c(2,mar))
  }
  

## Brain signatures
  b <- process(sigsBrain)
  colnames(b) <- substr(colnames(b), 1, 6)
  
  pdf(file = "Signature Correlations - Brain.pdf", height = 3, width = 3)
  my.heatmap(b, "Neu")
  my.heatmap(b, "Ast")
  my.heatmap(b, "Oli")
  my.heatmap(b, "Mic")
  my.heatmap(b, "End")
  dev.off()
  
## Pancreas
  p <- process(pancreas)
  
  pdf(file = "Signature Correlations - Pancreas.pdf", height = 2, width = 2)
  my.heatmap(p, "alpha")
  my.heatmap(p, "beta")
  dev.off()
  
## Heart
  h <- process(heart)
  colnames(h) <- gsub("ENCODE", "EN", colnames(h))
  
  pdf(file = "Signature Correlations - Heart.pdf", height = 2, width = 2)
  my.heatmap(h, "CM")
  my.heatmap(h, "SMC")
  my.heatmap(h, "FB")
  my.heatmap(h, "EC")
  dev.off()
  
## Legend
  pdf(file = "Signature Correlations - Legend.pdf", height = 3, width = 3)
  my.heatmap(h, "EC", key = TRUE, mar = 1)
  dev.off()
  
  
## Okay, look at brain markers...
  # define union of genes
  get.markers <- function(sigList) {
    genes <- lapply(sigList, rownames)
    genes <- do.call("c", genes)
    genes <- table(genes)
    genes <- names(genes)[which(genes == length(sigList))]
    
    # collect markers
    marks <- lapply(sigList, function(x) {
      x <- x[genes,]
      nTopFeatures(x, 100)
    })
    
    # rearrange markers
    dat <- lapply(marks, function(x) {
      y <- list()
      for (j in names(table(x$forCelltype))) y[[j]] <- x$EnsID[which(x$forCelltype == j)]
      return(y)  
    })
    return(dat)
  }
  
  marks <- list()
  
  marks$brain <- get.markers(sigsBrain)
  marks$pancreas <- get.markers(pancreas)
  marks$heart <- get.markers(heart)
  
  overlap.plot <- function(tissue, ct, leg = "none") {
    # collect markers for that ct across signatures
    p <- sapply(marks[[tissue]], function(x) x[ct])
    p <- do.call("cbind", p) # auto remove NAs
    
    # run pairwise intersections
    res <- as.data.frame(matrix(nrow = ncol(p), ncol = ncol(p)))
    colnames(res) <- rownames(res) <- colnames(p)
    
    for(j in colnames(res)) {
      for(k in rownames(res)) {
        res[k,j] <- length(intersect(p[,j], p[,k])) / nrow(p)
      }
    }
    
    res$Ct <- colnames(res) <- substr(rownames(res), 1, 2)
    res <- melt(res)
    
    res$variable <- factor(res$variable, levels = levels(as.factor(res$variable))[order(levels(as.factor(res$variable)))])
    res$Ct <- factor(res$Ct, levels = levels(as.factor(res$Ct))[order(levels(as.factor(res$Ct)))])
    
    ggplot(res, aes(x = Ct, y = variable, fill = value*100, label = value*100)) +
      geom_tile(colour = "black") +
      geom_text(size = 2.5) +
      scale_fill_carto_c(palette = "Temps", limits = c(0,100)) +
      theme_bw() +
      labs(title = ct) +
      theme(axis.title = invis, panel.border = element_blank(), panel.grid = element_blank(), legend.position = leg)
  }
  
  pdf(file = "Signature Correlations - Brain Overlap in Top 100 Markers.pdf", height = 2.7, width = 2.7)
  overlap.plot("brain", "Neurons")
  overlap.plot("brain", "Astrocytes")
  overlap.plot("brain", "Oligodendrocytes")
  overlap.plot("brain", "Microglia")
  overlap.plot("brain", "Endothelia")
  dev.off()
  
  pdf(file = "Signature Correlations - Pancreas Overlap in Top 100 Markers.pdf", height = 2, width = 2)
  overlap.plot("pancreas", "alpha")
  overlap.plot("pancreas", "beta")
  dev.off()
  
  pdf(file = "Signature Correlations - Heart Overlap in Top 100 Markers.pdf", height = 2, width = 2)
  overlap.plot("heart", "CM")
  overlap.plot("heart", "SMC")
  overlap.plot("heart", "FB")
  overlap.plot("heart", "EC")
  dev.off()
  
  pdf(file = "Signature Correlations - Legend 2.pdf", height = 4, width = 4)
  overlap.plot("heart", "CM", leg = "right") 
  dev.off()
