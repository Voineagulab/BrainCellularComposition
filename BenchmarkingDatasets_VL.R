################################################################################################################################ #
## Setup ----

## Start!
  rm(list = ls()); gc()

## Directory
  # root.dir <- "INSERT.YOUR.ROOT.DIRECTORY.PLEASE/"
  root.dir <- "/Volumes/Data1/PROJECTS/BrainCellularComposition/BrainCellularComposition_GitHub/"
  setwd(paste0(root.dir, "Results/BenchmarkingDatasets/VL/"))

## Load
  load("../../../Data/Preprocessed/Signatures - SNME.rda")
  load("../../../Data/Preprocessed/Signatures - Brain.rda")
  load("../../../Data/Preprocessed/Signatures - Brain (MuSiC).rda")
  load("../../../Data/Preprocessed/bench_CA_VL.rda")
  load("../../../Data/Preprocessed/snme_grad_VL.rda")
  source("../../../Scripts/Fun_Composition.R")
  source("../../../Scripts/Fun_Preprocessing.R")
  library(rcartocolor)

## Unload the CA-based SNME
  rm(bench_CA)

## Parameters
  algs <- c("DRS", "DTA", "CIB", "MUS")
  sigs <- names(sigsSNME)
  m1.order <- sigs
  m2.order <- c("VL", "NG", "CA", "LK", "TS")

## Modify signatures
  ## Remove $VL from sigsMuSiC, as it's in sigsSNME
    sigsMuSiC <- sigsMuSiC[-grep("VL", names(sigsMuSiC))]
    
  ## A quick change to Vel$merge2 (removing NRGN)
    sigsSNME$VL$merge2 <- sigsSNME$VL$merge2[,-grep("NRGN", colnames(sigsSNME$VL$merge2))]
    sigsSNME$VL$Meta$merge2[grep("NRGN", sigsSNME$VL$Meta$orig.celltype)] <- "Drop"
    
    
  ## Remake the CA signature (because you need CPM)
    # extract
    x <- sigsSNME$CA$merge2
    x$Neurons <- sigsSNME$CA$merge1$Neurons
    colnames(x) <- gsub("OPCs", "OPC", colnames(x)) 
    
    # add Microglia
    y <- sigsMuSiC$CA$counts[,which(sigsMuSiC$CA$meta$brain.ct == "Microglia")]
    y <- addSymbol(y)
    y <- apply(y, 2, function(z) {
      libSize <- 10 ^ 6 / sum(z) 
      z <- z * libSize
      return(z)
    })
    
    common <- intersect(rownames(y), rownames(x))
    y <- y[common,]
    x <- x[common,]
    x$Microglia <- rowMeans(y)
    
    # renormalise to counts per million
    x <- apply(x, 2, function(z) {
      libSize <- 10 ^ 6 / sum(z) 
      z <- z * libSize
      return(z)
    })
    
    sigsSNME$CA <- as.data.frame(x)
  
  ## Rename celltypes in non-VL signatures
    sigsSNME[-which(names(sigsSNME) == "VL")] <- lapply(sigsSNME[-which(names(sigsSNME) == "VL")], function(x) {
      # change OPCs to OPC
      colnames(x) <- gsub("OPCs", "OPC", colnames(x)) 
      
      # rename neuronal subtypes
      colnames(x) <- gsub("Excitatory", "Neurons_Exc", colnames(x)) 
      colnames(x) <- gsub("Inhibitory", "Neurons_Inh", colnames(x)) 
      colnames(x) <- gsub("Endothelia", "Endothelial", colnames(x)) 
      
      return(x)
    })
    
    sigsMuSiC <- lapply(sigsMuSiC, function(x) {
      
      x$meta$brain.ct <- gsub("OPCs", "OPC", x$meta$brain.ct) 
      x$meta$brain.ct2 <- gsub("OPCs", "OPC", x$meta$brain.ct2) 
      
      x$meta$brain.ct <- gsub("Endothelia", "Endothelial", x$meta$brain.ct) 
      x$meta$brain.ct2 <- gsub("Endothelia", "Endothelial", x$meta$brain.ct2) 
      
      x$meta$brain.ct <- gsub("Excitatory", "Neurons_Exc", x$meta$brain.ct) 
      x$meta$brain.ct <- gsub("Inhibitory", "Neurons_Inh", x$meta$brain.ct) 

      return(x)
    })
    
 

  ## Remove End from CA and LK
    # music
    # remove <- grep("End", sigsMuSiC$CA$meta$brain.ct2)
    # sigsMuSiC$CA$counts <- sigsMuSiC$CA$counts[,-remove]
    # sigsMuSiC$CA$meta <- sigsMuSiC$CA$meta[-remove,]    
    # 
    # remove <- grep("End", sigsMuSiC$LK$meta$brain.ct2)
    # sigsMuSiC$LK$counts <- sigsMuSiC$LK$counts[,-remove]
    # sigsMuSiC$LK$meta <- sigsMuSiC$LK$meta[-remove,]    
    # 
    
    # sigsSNME
    sigsSNME$LK <- sigsSNME$LK[,-grep("End", colnames(sigsSNME$LK))]
    
  ## Finally, addSymbol to SigsMuSiC
    sigsMuSiC <- lapply(sigsMuSiC, function(x) {
      x$counts <- addSymbol(x$counts)
      return(x)
    })
    
## A quick function
  invis <- function() element_blank() # essentially: less to type. this removes an element from a ggplot!

################################################################################################################################ #
## Create new truth matrices ----
    
## Reflecting the merging of cell-sub-types
true <- list()
true$orig.celltype <- as.data.frame(bench_VL$true$orig.celltype)    
 

## Merge 1
x <- as.data.frame(true$orig.celltype)
true$merge1 <- data.frame(Astrocytes = x$`AST-FB` + x$`AST-PP`,
                          Endothelial = x$Endothelial,
                          Microglia = x$Microglia,
                          Neurons = rowSums(x[,c(4:11, 13:15)]),
                          Oligodendrocytes = x$Oligodendrocytes,
                          OPC = x$OPC)

## Merge 2
true$merge2 <- data.frame(Astrocytes = x$`AST-FB` + x$`AST-PP`,
                          Endothelial = x$Endothelial,
                          Microglia = x$Microglia,
                          Neurons_Exc = rowSums(x[,c(8:11, 13)]),
                          Neurons_Inh = rowSums(x[,4:7]),
                          # Neurons_NRGN = rowSums(x[,14:15]),
                          Oligodendrocytes = x$Oligodendrocytes,
                          OPC = x$OPC)

    
################################################################################################################################ #
## Deconvolve ----
    
## Setup
  est <- list()
  enrich <- list()
  complete <- list()
  
## Each level 
  
## Deconvolve using the VL signature and partial deconvolution
  est$VL <- list()
  
  for(k in names(sigsSNME$VL[4:6])) { # orig.celltype, merge1, and merge2
    
    # setup  
    est$VL[[k]] <- list()

    print(paste0(k, "_DRS"))
    est$VL[[k]]$DRS <- run.DRS(bench_VL$mixture$cpm, sigsSNME$VL[[k]])

    print(paste0(k, "_DTA"))
    est$VL[[k]]$DTA <- run.DTA(bench_VL$mixture$cpm, sigsSNME$VL[[k]], alg = "diff", q = 0.01)

    print(paste0(k, "_CIB"))
    est$VL[[k]]$CIB <- run.CIB(from.file = FALSE,
                                 sigObject = sigsSNME$VL[[k]],
                                 mixString = "bench_VL.txt")
    
    est$VL[[k]]$MUS <- run.music(mixture = bench_VL$mixture$counts, use.meta.column = k, drop.ct = FALSE,
                                signature = sigsSNME$VL$Full.counts, 
                                signature.meta = sigsSNME$VL$Meta)
    
    save(est, file = "Est Temp.rda")
  }


## For all other signatures
  ## Function
  run.deconv <- function(sig, subneurons = TRUE) {
    # setup
    s <- sigsSNME[[sig]]
    s1 <- s[,which(colnames(s) %in% colnames(true$merge1))]
    if (subneurons) s2 <- s[,which(colnames(s) %in% colnames(true$merge2))]
    
    e <- list()
    
    # deconvolve to the level of merge1
    e$merge1 <- list()
    e$merge1$DRS <- run.DRS(bench_VL$mixture$cpm, s1)
    e$merge1$DTA <- run.DTA(bench_VL$mixture$cpm, s1, alg = "diff", q = 0.01)
    print("Starting CIB merge1")
    e$merge1$CIB <- run.CIB(from.file = FALSE,
                            sigObject = s1,
                            mixString = "bench_VL.txt")
    
    if (sig %in% names(sigsMuSiC)) {
      e$merge1$MuSiC <- run.music(mixture = bench_VL$mixture$counts, 
                                use.meta.column = "brain.ct2", 
                                signature = sigsMuSiC[[sig]]$counts, 
                                signature.meta = sigsMuSiC[[sig]]$meta) 
    }
    
    # deconvolve to the level of merge2 (i.e. subneurons)
    if (subneurons) {
      e$merge2 <- list()
      e$merge2$DRS <- run.DRS(bench_VL$mixture$cpm, s2)
      e$merge2$DTA <- run.DTA(bench_VL$mixture$cpm, s2, alg = "diff", q = 0.01)
      print("Starting CIB merge2")
      e$merge2$CIB <- run.CIB(from.file = FALSE,
                              sigObject = s2,
                              mixString = "bench_VL.txt")
      if (sig %in% names(sigsMuSiC)) {
        e$merge2$MUS <- run.music(mixture = bench_VL$mixture$counts, 
                                  use.meta.column = "brain.ct", 
                                  signature = sigsMuSiC[[sig]]$counts, 
                                  signature.meta = sigsMuSiC[[sig]]$meta) 
      }
    }
    
    # fin
    return(e)
  }
    

  ## Apply
    for (j in names(sigsSNME)) {
      # skip VL
      if (j %in% c("VL")) next
      
      # set subneurons to FALSE for certain signatures
      if (j %in% c("F5", "IP", "MM", "DM")) {
        logi <- FALSE  
      } else {
        logi <- TRUE
      }
      
      print(paste0(j, ": ", Sys.time()))
      est[[j]] <- run.deconv(j, subneurons = logi)
      save(est, file = "Temp (NewEst).rda")
    }
    
  
    
## Run enrichment algorithms  
  enrich$xCell <- run.xCell(bench_VL$mixture$cpm)
  enrich$Blender <- run.Blender(bench_VL$mixture$cpm, OPCs = TRUE)
  
## Save
  save(est, enrich, file = "Final SNME VL Composition.rda")
  
################################################################################################################################ #
## Complete deconvolution ----



## Linseed function
  apply.linseed <- function(mixture, merge, title) {
    n <- ncol(true[[merge]])
    
    pdf(file = paste0("Linseed Plots (", title, ").pdf"), height = 2.5, width = 5)
    output <- run.linseed(mixture = mixture,
                          nCelltypes = n,
                          write.plots = TRUE,
                          write.data = TRUE) 
    output$Data$svdPlot() + 
      labs(y = "Cumulative Variance Explained", x = "Number of Dimensions")
    dev.off()
    
    output <- output$Transformed  
    return(output)
  }
  
## A coex function
  analyse.coex <- function(coex, merge) {
    # enrichment tests for ct assignment
    coex <- define.celltypes(x = coex, algorithm.type = "Coex", celltype.signature = sigsSNME$VL[[merge]], 
                             return.confidence = TRUE, nMarks.x = 10000, nMarks.sig = 100, bg = rownames(bench_VL$mixture$cpm))
    
    # reannotate     
    m <- match(colnames(coex$comp), names(coex$assignments))
    colnames(coex$comp) <- coex$assignments[m]
    coex <- coex$comp
    
    # return
    return(coex)
  }
    
## Linseed: on the standard mixtures
  complete$linseed <- apply.linseed(bench_VL$mixture$cpm, "merge2", "merge2")

## Linseed: on the gradient mixtures
  complete$linseed.grad <- apply.linseed(grad_VL$confound.p, "merge2", "merge2.grad")
    
  

## Coex: on the standard mixtures
  ## Run
    coex.standardOutput <- run.coex(mixture = bench_VL$mixture$cpm, 
                                    signature = sigsSNME$VL$merge2, 
                                    only.threshold = FALSE, 
                                    sft = "auto", 
                                    output.all = TRUE)
  
  ## Extract information about merge1 cell-type abundance
    complete$coex <- analyse.coex(coex.standardOutput, "merge2")
    
## Coex: on the gradient mixtures
  ## Run
    coex.gradientOutput <- run.coex(mixture = grad_VL$confound.p,
                                      signature = sigsSNME$VL$merge2,
                                      only.threshold = FALSE,
                                      sft = "auto", 
                                      output.all = TRUE)
    
  ## Analyse
    complete$coex.grad <- analyse.coex(coex.gradientOutput, "merge2")
  
## Save
  save(complete, file = "Complete Deconvolution Estimates.rda")

 

################################################################################################################################ #
## Formal plots of deconvolution using Vel signatures ----

## Some plotting parameters
  ## A list of plotting parameters
  plot.information <- lapply(est$VL, function(w) {
    x <- colnames(w[[1]])
    
    names(x)[x == "Neurons"] <- "8Neu"
    names(x)[grep("AST|Ast", x)] <- "1Ast"
    names(x)[grep("^L|mat|Exc", x)] <- "5Neu-Exc"
    names(x)[grep("IN|Inh", x)] <- "6Neu-Inh"
    names(x)[grep("NRGN", x)] <- "7Neu-NRGN"
    names(x)[grep("^Oli", x)] <- "4Oli"
    names(x)[grep("^OPC", x)] <- "4OPC"
    names(x)[grep("End", x)] <- "2End"
    names(x)[grep("Mic", x)] <- "3Mic"
    
    x <- x[order(names(x))]
    return(x)
  })
  
  plot.information$universal.order <- c("1Ast", "2End", "3Mic", "4Oli", "4OPC", "5Neu-Exc", "6Neu-Inh", "7Neu-NRGN", "8Neu")
  plot.information$universal.palette <- c("#009392", "#39b185", "#9ccb86", "#f6edbd", "#e9e29c", "#eeb479", "#e88471", "#cf597e", "#e88471")  # this is carto_pal(7, "Temps"), with carto_pal(7, "Geyser")[4] at index 4, and carto_pal(7, "Temps")[6] also stuck at the end
  # plot.information$universal.palette <- c(plot.information$universal.palette, plot.information$universal.palette[6])
  names(plot.information$universal.palette) <- plot.information$universal.order
  
  
## Calculate stats
  
  errors <- list()
  for(k in sigs) {
    errors[[k]] <- list()
    for (j in names(est[[k]])) {
      errors[[k]][[j]] <- lapply(est[[k]][[j]], function(x) {
        return(t(write.stats(t = true[[j]], e = x, "Alg", TRUE)))
      })
    }  
  }
  
## Plot
  plot.formal.nmae <- function(merge) {
    
    
    # parameters: based on merge
    if (merge == "merge1") {
      l <- m1.order
      nrow <- 2
      leg <- "right"
    } else {
      l <- m2.order
      nrow <- 1
      leg <- "right"
    }
    
    # parameters
    b <- c(0,0.2,0.5,0.995,2,100)
    lab <- c("< 0.2", "0.2 - 0.5", "0.5 - 1", "1 - 2", "> 2")
    fills <- rev(c("#008080","#70a494","#b4c8a8","#de8a5a","#ca562c"))

    
    
    plot.data <- melt(errors)
    plot.data <- plot.data[which(plot.data$L2 == merge),]
    plot.data <- plot.data[which(plot.data$L1 == "VL"),]
    plot.data <- plot.data[which(plot.data$Var2 == paste0("nmae_Alg")),]
    colnames(plot.data) <- c("Celltype", "Metric", "Value", "Algorithm", "merge", "Signature")
    
    plot.data$Value <- round(plot.data$Value,3)
    plot.data$nmae <- plot.data$Value
    
    
    plot.data$nmae <- cut(plot.data$nmae, 
                              breaks = b,
                              labels = lab)
    
    # plot
    ggplot(plot.data, aes(x = Celltype, y = Algorithm, fill = nmae, label = Value)) +
      geom_tile(colour = "black") +
      geom_text(size = 2.8) +
      scale_fill_manual(values = fills) +
      theme_bw() +
      theme(legend.position = leg, panel.border = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
            axis.title.x = element_blank()) 
  }
 
  pdf(file = "NMAE - merge1.pdf", height = 3.5, width = 3.5)
  plot.formal.nmae("merge1") + theme(legend.position = "none")
  dev.off()
  
  pdf(file = "NMAE - merge2.pdf", height = 3.5, width = 4)
  plot.formal.nmae("merge2") + theme(legend.position = "none", axis.title.y = element_blank())
  dev.off()

  pdf(file = "NMAE - All Clusters.pdf", height = 3.5, width = 8)
  plot.formal.nmae("orig.celltype") 
  dev.off()

  
## Barplot
  plot.formal.barplot <- function(sig, merge, add.enrich.algs = TRUE, ncol = 1, lower.ylim = NA) {
    dat <- sapply(est[[sig]][[merge]], function(x) { diag(cor(x, true[[merge]][,colnames(x)])) })
    dat <- as.data.frame(dat)
    
    # add xCell and Blender
    if (add.enrich.algs) {
      dat$xCell <- NA
      dat["Astrocytes", "xCell"] <- cor(enrich$xCell$Astrocytes, true$merge1$Astrocytes)
      dat["Neurons", "xCell"] <- cor(enrich$xCell$Neurons, true$merge1$Neurons)
      
      dat$Blender <- NA
      dat["Astrocytes", "Blender"] <- cor(enrich$Blender$Astrocytes, true$merge1$Astrocytes)
      dat["Neurons", "Blender"] <- cor(enrich$Blender$Neurons, true$merge1$Neurons)
      dat["Oligodendrocytes", "Blender"] <- cor(enrich$Blender$Oligodendrocytes, true$merge1$Oligodendrocytes)
      dat["OPC", "Blender"] <- cor(enrich$Blender$OPC, true$merge1$OPC)
      dat["Endothelial", "Blender"] <- cor(enrich$Blender$Endothelia, true$merge1$Endothelial)
      dat["Microglia", "Blender"] <- cor(enrich$Blender$Microglia, true$merge1$Microglia)
    }
    
    dat$Celltype <- rownames(dat)
    dat <- melt(dat)
    
    colnames(dat)[2] <- "Algorithm"
    dat$Algorithm <- factor(dat$Algorithm, levels = rev(c("CIB", "DRS", "DTA", "MUS", "Blender", "xCell")))
    levels(dat$Algorithm) <- c("xCell", "Blender", "MuSiC", "dtangle", "DeconRNASeq",  "CIBERSORT")
    
    # recolour...
    m <- match(dat$Celltype, plot.information[[merge]])
    dat$Class <- names(plot.information[[merge]][m])
    
    # ...and reorder the x-axis
    dat$Celltype <- factor(dat$Celltype, levels = plot.information[[merge]])
    
    # now plot!
    ggplot(dat, aes(x = Celltype, y = value, fill = Class)) +
      geom_col(position = "dodge", width = 0.7) +
      theme_bw() +
      geom_hline(yintercept = 0) +
      facet_wrap(~Algorithm, ncol = ncol, strip.position = "top") +
      scale_y_continuous(limits = c(lower.ylim, 1), expand = c(0,0)) +
      labs(y = "Estimated vs. True (r)") +
      scale_fill_manual(values = plot.information$universal.palette) +
      geom_hline(yintercept = 0.8, linetype = 2) +
      # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.border = element_blank(), axis.line.y = element_line(),
      #       axis.title.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank())
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), axis.line.y = element_line(),
            axis.title.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank())
    
  }
  
  y <- -0.05
  
  pdf(file = "Correlation Barplot - merge1.pdf", height = 6.2, width = 2.3)
  plot.formal.barplot("VL", "merge1", add.enrich.algs = FALSE, lower.ylim = y) 
  dev.off()
  
  pdf(file = "Correlation Barplot - merge2.pdf", height = 6.16, width = 2.3)
  plot.formal.barplot("VL", "merge2", add.enrich.algs = FALSE, lower.ylim = y) + labs(y = "") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  dev.off()
  
  pdf(file = "Correlation Barplot - All Clusters.pdf", height = 6.18, width = 3.4)
  plot.formal.barplot("VL", "orig.celltype", add.enrich.algs = FALSE, lower.ylim = y) + labs(y = "") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  dev.off()
  

## On the enrichment algorithms...
  plot.formal.enrich.barplot <- function(alg) {
    if (alg == "xCell") {
      dat <- data.frame(Celltype = c("Astrocytes", "Neurons"),
                        value = c(cor(enrich$xCell$Astrocytes, true$merge1$Astrocytes),
                                  cor(enrich$xCell$Neurons, true$merge1$Neurons)))
    }
    
    if (alg == "Blender") {
      dat <- enrich$Blender
      colnames(dat)[4] <- "Endothelial"
      dat <- diag(cor(dat, true$merge1[,colnames(dat)]))
      dat <- as.data.frame(dat)
      dat$Celltype <- rownames(dat)
      colnames(dat)[1] <- "value"
    }
    
    # recolour...
    m <- match(dat$Celltype, plot.information$merge1)
    dat$Class <- names(plot.information$merge1[m])
    
    # ...and reorder the x-axis
    dat$Celltype <- factor(dat$Celltype, levels = plot.information$merge1)

    # now plot!
    ggplot(dat, aes(x = Celltype, y = value, fill = Class)) +
      geom_col(position = "dodge", width = 0.7) +
      theme_bw() +
      geom_hline(yintercept = 0) +
      # facet_wrap(~Algorithm, ncol = ncol, strip.position = "top") +
      scale_y_continuous(limits = c(-0.07, 1), expand = c(0,0)) +
      labs(y = "Estimated vs. True (r)") +
      scale_fill_manual(values = plot.information$universal.palette) +
      geom_hline(yintercept = 0.8, linetype = 2) +
      labs(title = alg) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.border = element_blank(), axis.line.y = element_line(),
            axis.title.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5))
  
  }
  
  pdf(file = "Final Barplot - xCell.pdf", height = 1.78, width = 1.2)
  plot.formal.enrich.barplot("xCell") + labs(y = "") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  dev.off()
  
  pdf(file = "Final Barplot - Blender.pdf", height = 2, width = 2.3)
  plot.formal.enrich.barplot("Blender")
  dev.off()
  
  

## Formal exploration of factors affecting performance
  ## First, quantify accuracy
    plot.data <- sapply(est$VL$orig.celltype, function(x) {
      x <- x[,colnames(true$orig.celltype)]
      diag(cor(x, true$orig.celltype)) 
      })
    
    
    plot.data <- plot.data > 0.8 # binary colour for accuracy
    plot.data <- as.data.frame(plot.data)  
    
  ## Add celltype information
    # plot.data$Celltype <- rownames(plot.data)
    m <- match(rownames(plot.data), plot.information$orig.celltype)
    plot.data$Class <- names(plot.information$orig.celltype)[m]
    plot.data$Class <- substr(plot.data$Class, 2, 100)
    plot.data$Class <- gsub("Neu-", "", plot.data$Class)
  

  ## Add abundance
    plot.data$Abundance <- colMeans(true$orig.celltype)

  ## Add each celltype maximum non-self correlation
    cor <- cor(sigsSNME$VL$orig.celltype, method = "s")
    plot.data$Sim <- apply(cor, 1, function(x) max(x[which(x != max(x))]))
    
    # jitter
    plot.data$SimJitter <- plot.data$Sim + 0.01
    plot.data$SimJitter[which(plot.data$Class == "NRGN")] <- plot.data$SimJitter[which(plot.data$Class == "NRGN")] - 0.02
    plot.data$SimJitter[9] <- plot.data$SimJitter[9] - 0.02
    plot.data$Label <- letters[1:nrow(plot.data)]
    
    leg <- data.frame(Label = plot.data$Label, Celltype = rownames(plot.data))
    
  ## Plot
    
    plot.data <- melt(plot.data, id.vars = c("Class", "Sim", "SimJitter", "Abundance", "Label"))
    plot.data$variable <- as.character(plot.data$variable)
  
  pdf(file = "Celltype Similarity vs. Abundance.pdf", height = 6, width = 8)
  ggplot(plot.data, aes(x = (Abundance*100), y = Sim, colour = value, label = Label)) +
    geom_text(size = 4 ) +

    scale_colour_manual(values = c("firebrick1", "dodgerblue1")) +
    facet_wrap(~variable) +
    scale_y_continuous(limits = c(NA, 1)) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_continuous(trans = pseudo_log_trans(sigma = 1), limits = c(0, 50), 
                       expand = c(0,0), breaks = c(0, 1, 2, 5, 10, 25, 40), labels = function(x) paste0(x, "%")) +
    labs(y = "Highest Non-self Correlation to Other Cell-types (Rho)", x = "Abundance in Mixture") 
    dev.off()
    
   

  
################################################################################################################################ #
## Formal plots of complete deconvolution ----  
  
  
## First, fix "true"
  # by "true", the truth matrices need a column for pan-neuronal abundance
  complete.true <- list(A = grad_VL$meta.prop,
                        B = grad_VL$meta.prop)
    
    complete.true$A$Neurons <- complete.true$A$Neurons_Exc + complete.true$A$Neurons_Inh
    complete.true$A <- complete.true$A[,-grep("_", colnames(complete.true$A))]
    
    # and rename columns
    complete.true <- lapply(complete.true, function(x) {
      colnames(x) <- gsub("OPCs", "OPC", colnames(x))
      return(x)
    })
    
## On Linseed
  ## Heatmap function
    linseed.heatmap <- function(e, t, filename) {
      plot.data <- cor(e, t, method = "p")
      
      plot.data <- as.data.frame(plot.data)
      plot.data$Linseed <- LETTERS[1:nrow(plot.data)]
      plot.data <- melt(plot.data)
      colnames(plot.data)[3] <- "Correlation"
      plot.data$variable <- gsub("Neurons_", "", plot.data$variable)
      
      # pdf(file = filename, height = 2.7, width = 2.5)
      print(ggplot(plot.data, aes(x = variable, y = Linseed, fill = Correlation, label = round(Correlation, 2))) +
              geom_tile(colour = "black") +
              geom_text(size = 2.5) +
              theme_bw() +
              scale_fill_carto_c(palette = "ArmyRose", limits = c(-1,1)) +
              scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
              labs(x = "True Cell-type", y = "Linseed Cell-type") +
              NoLegend() +
              theme(panel.border = invis(), axis.ticks = invis(), panel.grid = invis(), axis.text.x = element_text(angle = 90, hjust = 1)))
      # dev.off()
    }  
    
  ## Apply
    linseed.heatmap(complete$linseed, true$merge2, "Complete Linseed Heatmap, merge2.pdf")
    linseed.heatmap(complete$linseed.grad, complete.true$B, "Complete Linseed Heatmap, gradientB, merge2.pdf")
    
  ## Scatterplot function
    complete.scatter <- function(t, l, ct, colour, ylab, axis.min = 0, abline = TRUE) {
      r <- round(cor(t, l), 2)
      label <- paste0("r=", r)
      p <- qplot(t, l, colour = I(colour)) +
        theme_bw() +
        # geom_smooth(method = "lm", se = FALSE) +
        scale_y_continuous(limits = c(axis.min, NA), expand = c(0,0)) +
        labs(y = ylab, x = paste0("True ", ct, "\n", label)) +
        scale_x_continuous(limits = c(axis.min, NA), expand = c(0,0)) +
        theme(panel.border = invis(), axis.line = element_line()) +
        NoLegend()
      
      if (abline) {
        return(p + geom_abline(slope = 1, intercept = 0, linetype = 2))
      } else {
        return(p)
      }
    }
    
    
    
## Coex
  ## Function
    coex.barplot <- function(e, u, t = "merge1", filename) {
      dat <- diag(cor(e, u[,colnames(e)]))
      
      absent <- colnames(u)[-which(colnames(u) %in% colnames(e))]
      zeros <- rep(0, length(absent))
      names(zeros) <- absent
      dat <- c(dat, zeros)
      
      m <- match(names(dat), plot.information[[t]])
      
      
      
      dat <- data.frame(Cor = dat, Celltype = names(dat), Class = names(plot.information[[t]][m]))
      dat$Celltype <- gsub("Neurons_", "", dat$Celltype)
      
      # pdf(file = filename, height = 2.5, width = 2.5)
      print(ggplot(dat, aes(x = Celltype, y = Cor, fill = Class, colour = Class)) +
              geom_col() +
              theme_bw() +
              scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
              scale_fill_manual(values = plot.information$universal.palette) +
              scale_colour_manual(values = plot.information$universal.palette) +
              scale_y_continuous(limits = c(NA, 1)) +
              labs(y = "Pearson Correlation") +
              geom_hline(yintercept = 0.8, linetype = 2) +
              theme(panel.border = invis(), axis.line = element_line(), legend.position = "none"))
      # dev.off()
    }
    
    
    coex.barplot(e = complete$coex, u = true$merge2, t = "merge2", filename = "Complete Coex Barplot, merge2.pdf")
    coex.barplot(e = complete$coex.grad, u = complete.true$B, t = "merge2", filename = "Complete Coex Barplot, gradient, merge2.pdf")
 
   
  

################################################################################################################################ #
## Formal plots of deconvolution using mismatched signatures ----
  
  
## Scatterplots
  ## Function
  plot.formal.mismatch.scatter <- function(sig, alg, merge, no.xlab = TRUE, no.ylab = TRUE) {

    e <- est[[sig]][[merge]][[alg]]
    t <- true[[merge]]
    t <- t[,colnames(e)]
    
    # reformat for plotting
    e <- melt(e)
    t <- melt(t)
    
    plot.data <- data.frame(Celltype = e$variable,
                            Estimate = e$value,
                            True = t$value)
    
    # add colouring information
    m <- match(plot.data$Celltype, plot.information[[merge]])
    plot.data$Class <- names(plot.information[[merge]][m])
    
    # shuffle order
    plot.data <- plot.data[sample(rownames(plot.data), nrow(plot.data), replace = FALSE),]
    
    # add extra annotation to the title
    lab <- sig
    if (sig == "VL") lab <- paste0(lab, " (Matched)")
    
    # merge-specific parameters
    if (merge == "merge1") {
      sigma <- 0.5
      annot.x <- 5
      annot.y <- 80
    } else {
      sigma <- 5 
      annot.x <- 20
      annot.y <- 95
    }
    
    # plot
    plot <- ggplot(plot.data, aes(x = True*100, y = Estimate*100, colour = Class)) +
      geom_point(alpha = 1, size = 2) +
      theme_bw() +
      annotate("text", x = annot.x, y = annot.y, label = lab) +
      scale_colour_manual(values = plot.information$universal.palette) +
      scale_x_continuous(trans = pseudo_log_trans(sigma = sigma), limits = c(0, 120), expand = c(0,0), breaks = c(1, 2, 5, 10, 25, 50, 100)) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = sigma), limits = c(0, 120), expand = c(0,0), breaks = c(1, 2, 5, 10, 25, 50, 100)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
      theme(legend.position = "none", panel.grid = invis, panel.border = invis,
            plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill = "grey80")) +
      labs(y = paste0("Estimated Proportion"), x = "True Proportion")
    
    if (no.xlab) plot <- plot + theme(axis.title.x = invis, axis.text.x = invis, axis.ticks.x = invis)
    if (no.ylab) plot <- plot + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
    
    return(plot)
  }
  
  ## Apply function
  for (j in algs) {
    
    pdf(file = paste0("Mismatch Signature - ", j, " merge2.pdf"), height = 2.5, width = 8)
     pA <- plot.formal.mismatch.scatter(sig = "VL", alg = j, merge = "merge2", no.ylab = FALSE, no.xlab = FALSE)
    pB <- plot.formal.mismatch.scatter(sig = "NG", alg = j, merge = "merge2", no.xlab = FALSE)
    pC <- plot.formal.mismatch.scatter(sig = "CA", alg = j, merge = "merge2", no.xlab = FALSE)
    pD <- plot.formal.mismatch.scatter(sig = "LK", alg = j, merge = "merge2", no.xlab = FALSE)
    pE <- plot.formal.mismatch.scatter(sig = "TS", alg = j, merge = "merge2", no.xlab = FALSE)
    pF <- get_legend(pA + theme(legend.position = "right"))
    # pF <- plot.formal.mismatch.scatter(sig = "Dar", alg = "CIB", merge = "merge1", no.xlab = FALSE) # because it's merge1
    # print(plot_grid(pA, pB, pC, pD, pE, pF, nrow = 2, rel_widths = c(1, 0.8, 0.8), rel_heights = c(1, 1.2)))
    print(plot_grid(pA, pB, pC, pD, pE, nrow = 1, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8)))
    dev.off()
    
    if (j == "MUS") next
    
    pdf(file = paste0("Mismatch Signature - ", j, " merge1.pdf"), height = 2.8, width = 8)
    pA <- plot.formal.mismatch.scatter(sig = "VL", alg = j, merge = "merge1", no.ylab = FALSE)
    pB <- plot.formal.mismatch.scatter(sig = "NG", alg = j, merge = "merge1")
    pC <- plot.formal.mismatch.scatter(sig = "CA", alg = j, merge = "merge1")
    pD <- plot.formal.mismatch.scatter(sig = "LK", alg = j, merge = "merge1")
    pE <- plot.formal.mismatch.scatter(sig = "TS", alg = j, merge = "merge1")
    pF <- plot.formal.mismatch.scatter(sig = "F5", alg = j, merge = "merge1", no.xlab = FALSE, no.ylab = FALSE) + labs(y = "")
    pG <- plot.formal.mismatch.scatter(sig = "IP", alg = j, merge = "merge1", no.xlab = FALSE)
    pH <- plot.formal.mismatch.scatter(sig = "MM", alg = j, merge = "merge1", no.xlab = FALSE)
    pI <- plot.formal.mismatch.scatter(sig = "DM", alg = j, merge = "merge1", no.xlab = FALSE)
    pJ <- get_legend(pA + theme(legend.position = "right") + guides(colour = guide_legend(ncol = 2)))  
    print(plot_grid(pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ, nrow = 2, ncol = 5, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8), rel_heights = c(1, 1.2)))  
    dev.off()
  }
  
  
  ## Separately plot MUS
   pdf(file = paste0("Mismatch Signature - MUS merge1.pdf"), height = 2.8, width = 5)
    pA <- plot.formal.mismatch.scatter(sig = "VL", alg = "MUS", merge = "merge1", no.ylab = FALSE)
    pB <- plot.formal.mismatch.scatter(sig = "NG", alg = "MUS", merge = "merge1")
    pC <- plot.formal.mismatch.scatter(sig = "CA", alg = "MUS", merge = "merge1")
    pD <- plot.formal.mismatch.scatter(sig = "LK", alg = "MUS", merge = "merge1", no.xlab = FALSE, no.ylab = FALSE) + labs(y = "")
    pE <- plot.formal.mismatch.scatter(sig = "TS", alg = "MUS", merge = "merge1", no.xlab = FALSE)
    # pI <- plot.formal.mismatch.scatter(sig = "DM", alg = "MUS", merge = "merge1", no.xlab = FALSE)
    pJ <- get_legend(pA + theme(legend.position = "right") + guides(colour = guide_legend(ncol = 2)))  
    print(plot_grid(pA, pB, pC, pD, pE, pI, pJ, nrow = 2, ncol = 3, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8), rel_heights = c(1, 1.2)))  
  dev.off()
  
  
  
  
## Alternative view: heatmap!!
  ## Function
    stat.heatmap <- function(merge, metric) {
      # parameters: based on merge
      if (merge == "merge1") {
        l <- m1.order
        nrow <- 2
        leg <- c(0.9, 0.2)
      } else {
        l <- m2.order
        nrow <- 1
        leg <- "bottom"
      }
      
      # parameters: based on metric
      if (metric == "nmae") {
        b <- c(0,0.2,0.5,0.995,2,100)
        lab <- c("< 0.2", "0.2 - 0.5", "0.5 - 1", "1 - 2", "> 2")
        fills <- rev(c("#008080","#70a494","#b4c8a8","#de8a5a","#ca562c"))
      } else {
        b <- c(-10,-1,0.5,0.8,0.9,0.95,1)
        lab <- c("NA", "< 0.5", "0.5 - 0.8", "0.8 - 0.9", "0.9 - 0.95", "> 0.95")
        fills <- c("grey50", "#008080","#70a494","#b4c8a8","#de8a5a","#ca562c")
      }
      
      
      plot.data <- melt(errors)
      plot.data <- plot.data[which(plot.data$L2 == merge),]
      plot.data <- plot.data[which(plot.data$Var2 == paste0(metric, "_Alg")),]
      colnames(plot.data) <- c("Celltype", "Metric", "value", "Algorithm", "merge", "Signature")
      
      plot.data[,metric] <- plot.data$value
      plot.data$value <- round(plot.data$value,2)
      if (metric == "r") {
        plot.data$value <- as.character(plot.data$value)
        plot.data$value[which(is.na(plot.data$value))] <- "NA"
        # plot.data$bin[which(is.na(plot.data$bin))] <- -2  
        
        plot.data[which(is.na(plot.data[,metric])),metric] <- -2  
      }
      
      plot.data[,metric] <- cut(plot.data[,metric], 
                                breaks = b,
                                labels = lab)
      
  
      plot.data$Signature <- factor(plot.data$Signature, levels = l)
      
      # plot
      ggplot(plot.data, aes_string(x = "Algorithm", y = "Celltype", fill = metric, label = "value")) +
        geom_tile(colour = "black") +
        geom_text(size = 2.5) +
        facet_wrap(~Signature, nrow = nrow) +
        scale_fill_manual(values = fills) +
        theme_bw() +
        theme(legend.position = leg, axis.title.y = invis) 
    }
  
  ## Apply function
    pdf(file = "Mismatch Signature - Heatmap, merge1.pdf", height = 4, 8)
    stat.heatmap("merge1", "r")
    stat.heatmap("merge1", "nmae")
    dev.off()
    
    pdf(file = "Mismatch Signature - Heatmap, merge2.pdf", height = 4, 8)
    stat.heatmap("merge2", "r")
    stat.heatmap("merge2", "nmae")
    dev.off()
    

