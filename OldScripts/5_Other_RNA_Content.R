## In this script, we explore the role of RNA content


################################################################################################################################ #
## Setup ---- 

## Load
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Other/")
load("../../../Data/Preprocessed/Signatures - Brain (MuSiC).rda")
source("../../../Scripts/Fun_Composition.R")
source("../../../Scripts/Fun_Preprocessing.R")


################################################################################################################################ #
## First: is there an effect of cell-type on library size (if UMIs, it's a proxy for poly-A+ RNA content) ---- 

## Libsize
plot.data <- data.frame(libsize = colSums(sigsMuSiC$VL$counts),
                        celltype = substr(sigsMuSiC$VL$meta$brain.ct, 1, 3))

## Plot

pdf(file = "RNA Content - Lib Sizes Across Ct, VL.pdf", height = 2.5, width = 4)
ggplot(plot.data, aes(x = celltype, colour = celltype, y = log10(libsize))) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "log10(UMIs per cell)", x = "Celltype") +
  theme(legend.position = "none")
dev.off()


################################################################################################################################ #
## Simulate for deconvolution: is RNA or cell proportion estimated? ---- 

## Generate pseudobulk mixtures for each individual
  dat <- sigsMuSiC$VL$counts
  meta <- sigsMuSiC$VL$meta
  rm(sigsMuSiC)
  gc()

  pb <- list()
  props <- list()
  for (j in unique(meta$Individual)) {
    # index of all cells for a given individual
    print(j)
    w <- which(meta$Individual == j)
    
    # get pseudobulk
    exp <- rowSums(dat[,w])
    exp <- exp / (sum(exp) / 10^6)
    pb[[as.character(j)]] <- exp
    
    # rna and cell proportions
    rna <- aggregate(colSums(dat[,w]) ~ meta$brain.ct[w], FUN = sum) 
    rna <- rna[,2] / sum(dat[,w]) # fraction of UMIs coming from each cell-type
    
    cell <- table(meta$brain.ct[w])
    cell <- cell / sum(cell) # fraction of cells coming from each cell-type
    
    props[[as.character(j)]] <- data.frame(Cell = as.numeric(cell),
                                           RNA = rna,
                                           Celltype = names(cell))
  }

  pb <- do.call("cbind", pb)
  pb <- as.data.frame(pb)
  write.table(pb, file = "RNA Content - VL Individual Mixtures.txt", sep = "\t", quote = FALSE)
  
  
  props <- do.call("rbind", props)
  props$Individual <- as.character(substr(rownames(props), 1, 4))
  
  cors <- list()
  for (j in unique(props$Celltype)) {
    x <- props[which(props$Celltype == j),]
    cors[[j]] <- round(cor(x$Cell, x$RNA), 2)
  }
  do.call("c", cors)
  
  # scatterplot cell and rna proportions
  pdf(file = "RNA Content - VL Individual Proportions (2).pdf", height = 2.5, width = 4)
  ggplot(props, aes(x = Cell, y = RNA, colour = Celltype)) + 
    geom_point() +
    scale_x_continuous(limits = c(0,0.8), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,0.8), expand = c(0,0)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    theme_bw() +
    labs(x = "Cell Proportion In Mixture", y = "RNA Proportion In Mixture") +
    theme(panel.border = invis, axis.line = element_line())
  dev.off()
  
  # estimate proportion using CIB
  sig <- list()
  for (j in unique(meta$brain.ct)) {
    # index of all cells for a given individual
    print(j)
    w <- which(meta$brain.ct == j)
    
    # get pseudobulk
    exp <- rowSums(dat[,w])
    exp <- exp / (sum(exp) / 10^6)
    sig[[j]] <- exp
  }
  
  sig <- do.call("cbind", sig)
  sig <- as.data.frame(sig)
  
  write.CIB(data = pb, dir = "RNA Content - VL CIB Input.txt")
  cib <- run.CIB(directoryCIB = "", mixString = "RNA Content - VL CIB Input.txt", sigObject = sig, from.file = FALSE)
  write.table(cib, file = "RNA Content - VL CIB Output.txt", sep = "\t")

  ## Plot estimates versus true
    pdf(file = "RNA Content - VL Individual Proportions Scatterplot (2).pdf", height = 2.5, width = 4)
    
    p <- cib  
    p$Sample <- rownames(p)
    p <- melt(p)
    
    m <- match(paste0(props$Individual, props$Celltype), paste0(p$Sample, p$variable))
    props$Est <- p$value[m]    
  
    # cell vs estimate
    r <- cor(props$Cell, props$Est)
    r <- round(r, 2)
    ggplot(props, aes(x = Cell, y = Est, colour = Celltype)) + 
      geom_point() +
      scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      theme_bw() +
      labs(x = paste0("Cell Proportion In Mixture\nr=", r), y = "CIBERSORT Estimate") +
      theme(panel.border = invis, axis.line = element_line())
    
    # rna vs estimate
    r <- cor(props$RNA, props$Est)
    r <- round(r, 2)
    
    ggplot(props, aes(x = RNA, y = Est, colour = Celltype)) + 
      geom_point() +
      scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      theme_bw() +
      labs(x = paste0("RNA Proportion In Mixture\nr=", r), y = "CIBERSORT Estimate") +
      theme(panel.border = invis, axis.line = element_line())
  
    dev.off()
  
  ## Goodness of fit
    # in cells
    prop.cell <- dcast(props, Individual~Celltype, value.var = "Cell")
    prop.cell[(is.na(prop.cell))] <- 0
    gof.cell <- write.gof.v2(measuredExp = pb, estimatedComp = prop.cell, signatureUsed = sig, returnPred = FALSE)
    
    # in rna
    prop.rna <- dcast(props, Individual~Celltype, value.var = "RNA")
    prop.rna[(is.na(prop.rna))] <- 0
    gof.rna <- write.gof.v2(measuredExp = pb, estimatedComp = prop.rna, signatureUsed = sig, returnPred = FALSE)
    
    # plot
    pdf(file = "RNA Content - VL Individual GoF (2).pdf", height = 2, width = 2)
    qplot(gof.cell$r, gof.rna$r) + 
      scale_x_continuous(limits = c(0.92,1), expand = c(0,0)) +
      scale_y_continuous(limits = c(0.92,1), expand = c(0,0)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      theme_bw() +
      labs(x = "Goodness of Fit (r)\n(Cell Proportions)", y = "Goodness of Fit (r)\n(RNA Proportions)") +
      theme(panel.border = invis, axis.line = element_line())
    dev.off()
