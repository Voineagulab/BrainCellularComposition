################################################################################################################################ #
## Setup ----

## Start!
  rm(list = ls())
  setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/snme/CA/")

## Load
  load("../../../../Data/Preprocessed/Signatures - SNME.rda")
  load("../../../../Data/Preprocessed/Signatures - Brain (MuSiC).rda")
  load("../../../../Data/Preprocessed/snme.rda")
  load("../../../../Data/Preprocessed/snme_grad_CA.rda")
  source("../../../../Scripts/Fun_Composition.R")
  source("../../../../Scripts/Fun_Preprocessing.R")
  library(rcartocolor)

## Unload the CA-based SNME
  snme <- snme$HCA

## Parameters
  algs <- c("DRS", "DTA", "CIB", "MUS")
  sigs <- names(sigsSNME)
  m1.order <- sigs
  m2.order <- c("CA", "NG", "VL", "LK", "TS")
  
## Modify signatures
   ## Remove $CA from sigsMuSiC, as it's in sigsSNME
    sigsMuSiC <- sigsMuSiC[-3]
    
  ## Add dummy columns to sigsMuSiC DM
    sigsMuSiC$DM$meta$brain.ct2 <- sigsMuSiC$DM$meta$brain.ct <- sigsMuSiC$DM$meta$orig.celltype
    
  ## Remove End and Mic
    sigsMuSiC <- lapply(sigsMuSiC, function(x) {
      remove <- grep("^End|^Mic", x$meta$brain.ct)
      
      x$counts <- x$counts[,-remove]
      x$meta <- x$meta[-remove,]
      
      return(x)
    })
    
  ## Finally, addSymbol to SigsMuSiC
    sigsMuSiC <- lapply(sigsMuSiC, function(x) {
      x$counts <- addSymbol(x$counts)
      return(x)
    })
  
    
## Rename celltypes
  sigsSNME[-c(1:2)] <- lapply(sigsSNME[-c(1:2)], function(x) {
    # rename neuronal subtypes
    colnames(x) <- gsub("Excitatory", "Neurons_Exc", colnames(x)) 
    colnames(x) <- gsub("Inhibitory", "Neurons_Inh", colnames(x)) 
    return(x)
  })
  
  sigsMuSiC <- lapply(sigsMuSiC, function(x) {
      x$meta$brain.ct <- gsub("Excitatory", "Neurons_Exc", x$meta$brain.ct) 
      x$meta$brain.ct <- gsub("Inhibitory", "Neurons_Inh", x$meta$brain.ct) 

      return(x)
    })
  
## Modify the formatting of the VL signature
  x <- sigsSNME$VL$merge2
  x$Neurons <- sigsSNME$VL$merge1$Neurons
  x <- x[,-grep("NRGN", colnames(x))]
  colnames(x)[6:7] <- c("Endothelia", "OPCs")
  
  sigsSNME$VL <- x

## A quick function
  invis <- element_blank() # essentially: less to type. this removes an element from a ggplot!

################################################################################################################################ #
## Create new truth matrices ----
    
## Reflecting the merging of cell-sub-types
true <- list()
true$orig.celltype <- as.data.frame(snme$true$orig.celltype)    

## Merge 1
x <- as.data.frame(true$orig.celltype)
true$merge1 <- data.frame(Neurons = rowSums(x[,grep("Exc|Inh", colnames(x))]),
                          Astrocytes = x$`Astro L1-6 FGFR3 SLC14A1`,
                          Oligodendrocytes = x$`Oligo L1-6 OPALIN`,
                          OPCs = x$`OPC L1-6 PDGFRA`)

## Merge 2
true$merge2 <- data.frame(Neurons_Exc = rowSums(x[,grep("Exc", colnames(x))]),
                          Neurons_Inh = rowSums(x[,grep("Inh", colnames(x))]),
                          Astrocytes = x$`Astro L1-6 FGFR3 SLC14A1`,
                          Oligodendrocytes = x$`Oligo L1-6 OPALIN`,
                          OPCs = x$`OPC L1-6 PDGFRA`)


    
################################################################################################################################ #
## Deconvolve ----
    

## Setup
  est <- list()
  enrich <- list()
  complete <- list()
  
## Deconvolve using the CA (matched) signatures
  est$CA <- list()
  
  for(k in names(sigsSNME$CA[-c(1:2)])) {
    
    # setup  
    est$CA[[k]] <- list()
    
    print(paste0(k, "_DRS"))
    est$CA[[k]]$DRS <- run.DRS(snme$mixture$cpm, sigsSNME$CA[[k]])
    
    print(paste0(k, "_DTA"))
    est$CA[[k]]$DTA <- run.DTA(snme$mixture$cpm, sigsSNME$CA[[k]], alg = "diff", q = 0.01)
    
    print(paste0(k, "_CIB"))
    est$CA[[k]]$CIB <- run.CIB(from.file = FALSE,
                               sigObject = sigsSNME$CA[[k]],
                               mixString = "snme_Vel.txt")
    
    print(paste0(k, "_MUS"))
    if (substr(k, 1, 4) == "drop" | (j == "CA" & k == "merge2")) {
      drop <- TRUE 
      print("Dropping ct in MUS")
    } else {
      drop <- FALSE
    }
    
    est$CA[[k]]$MUS <- run.music(mixture = snme$mixture$counts, use.meta.column = k, drop.ct = drop,
                                 signature = sigsSNME$CA$Full.counts, 
                                 signature.meta = sigsSNME$CA$Meta)
    
    save(est, file = "Est Temp.rda")
  }
  

## For non-matched (i.e. non-CA signatures)
  ## Function 
  run.deconv <- function(sig, subneurons = TRUE) {
    # setup
    s <- sigsSNME[[sig]]
    s1 <- s[,which(colnames(s) %in% colnames(true$merge1))]
    if (subneurons) s2 <- s[,which(colnames(s) %in% colnames(true$merge2))]
    
    e <- list()
    
    # deconvolve to the level of merge1
    e$merge1 <- list()
    e$merge1$DRS <- run.DRS(snme$mixture$cpm, s1)
    e$merge1$DTA <- run.DTA(snme$mixture$cpm, s1, alg = "diff", q = 0.01)
    print("Starting CIB merge1")
    e$merge1$CIB <- run.CIB(from.file = FALSE,
                            sigObject = s1,
                            mixString = "snme_HCA.txt")
    
    if (sig %in% names(sigsMuSiC)) {
      e$merge1$MUS <- run.music(mixture = snme$mixture$counts, 
                                use.meta.column = "brain.ct2", 
                                signature = sigsMuSiC[[j]]$counts, 
                                signature.meta = sigsMuSiC[[j]]$meta) 
    }
    
    # deconvolve to the level of merge2 (i.e. subneurons)
    if (subneurons) {
      e$merge2 <- list()
      e$merge2$DRS <- run.DRS(snme$mixture$cpm, s2)
      e$merge2$DTA <- run.DTA(snme$mixture$cpm, s2, alg = "diff", q = 0.01)
      print("Starting CIB merge2")
      e$merge2$CIB <- run.CIB(from.file = FALSE,
                              sigObject = s2,
                              mixString = "snme_HCA.txt")
      if (sig %in% names(sigsMuSiC)) {
        e$merge2$MUS <- run.music(mixture = snme$mixture$counts, 
                                  use.meta.column = "brain.ct", 
                                  signature = sigsMuSiC[[j]]$counts, 
                                  signature.meta = sigsMuSiC[[j]]$meta) 
      }
    }
    
    # fin
    return(e)
  }
  
  
  ## Apply
  for (j in names(sigsSNME)) {
    # skip CA
    if (j == "CA") next
    
    # set subneurons to FALSE for certain signatures
    if (j %in% c("F5", "IP", "MM", "DM")) {
      logi <- FALSE  
    } else {
      logi <- TRUE
    }
    
    print(paste0(j, ": ", Sys.time()))
    est[[j]] <- run.deconv(j, subneurons = logi)
    save(est, file = "Temp (Est).rda")
  }
  
  ## Note that MuSiC automatically dropped DM's OPCs, as it's only present in a single individual
  est$DM$merge1$MUS$OPCs <- 0
  est$DM <- est$DM[1]
  

## Enrichment
  enrich$xCell <- run.xCell(snme$mixture$cpm)
  enrich$Blender <- run.Blender(snme$mixture$cpm)
    
  
  save(est, file = "Prop (Nine Signatures, New).rda")
  save(enrich, file = "Enrichment.rda")
  
  
################################################################################################################################ #
## Complete Deconvolution ----
  
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
    coex <- define.celltypes(x = coex, algorithm.type = "Coex", celltype.signature = sigsSNME$CA[[merge]], 
                             return.confidence = TRUE, nMarks.x = 10000, nMarks.sig = 100, bg = rownames(snme$mixture$cpm))
    
    # reannotate     
    m <- match(colnames(coex$comp), names(coex$assignments))
    colnames(coex$comp) <- coex$assignments[m]
    coex <- coex$comp
    
    # return
    return(coex)
  }
    
## Linseed: on the standard mixtures
  complete$linseed.merge1 <- apply.linseed(snme$mixture$cpm, "merge1", "merge1")
  complete$linseed.merge2 <- apply.linseed(snme$mixture$cpm, "merge2", "merge2")

## Linseed: on the gradient mixtures
  complete$linseed.merge1.gA <- apply.linseed(snme.gradients$A$confound.p, "merge1", "merge1.gA")
  complete$linseed.merge2.gA <- apply.linseed(snme.gradients$A$confound.p, "merge2", "merge2.gA")
  
  complete$linseed.merge1.gB <- apply.linseed(snme.gradients$B$confound.p, "merge1", "merge1.gB")
  complete$linseed.merge2.gB <- apply.linseed(snme.gradients$B$confound.p, "merge2", "merge2.gB")
    
    

## Save incomplete
  save(complete, file = "Complete Deconvolution Estimates (In Progress).rda")
    
## Coex: on the standard mixtures
  ## Run
    coex.standardOutput <- run.coex(mixture = snme$mixture$cpm, 
                                    signature = sigsSNME$CA$merge1, 
                                    only.threshold = FALSE, 
                                    sft = "auto", 
                                    output.all = TRUE)
  
  ## Extract information about merge1 cell-type abundance
    complete$coex.merge1 <- analyse.coex(coex.standardOutput, "merge1")
    complete$coex.merge2 <- analyse.coex(coex.standardOutput, "merge2")
    
## Coex: on the gradient mixtures
  ## Run
    coex.gradientOutput.A <- run.coex(mixture = snme.gradients$A$confound.p,
                                      signature = sigsSNME$CA$merge1,
                                      only.threshold = FALSE,
                                      sft = "auto", 
                                      output.all = TRUE)
    
    coex.gradientOutput.B <- run.coex(mixture = snme.gradients$B$confound.p,
                                      signature = sigsSNME$CA$merge1,
                                      only.threshold = FALSE,
                                      sft = "auto", 
                                      output.all = TRUE)
    
    
  ## Analyse
    complete$coex.merge1.gA <- analyse.coex(coex.gradientOutput.A, "merge1")
    complete$coex.merge2.gA <- analyse.coex(coex.gradientOutput.A, "merge2")
    
    complete$coex.merge1.gB <- analyse.coex(coex.gradientOutput.B, "merge1")
    complete$coex.merge2.gB <- analyse.coex(coex.gradientOutput.B, "merge2")
    
  ## Save raw coex
  save(coex.gradientOutput.A, coex.gradientOutput.B, coex.standardOutput, file = "Coexpression Data.rda")
    
## Save
  save(complete, file = "Complete Deconvolution Estimates.rda")

## Note that all files were manually moved to ./CompleteDeconv for decluttering purposes

################################################################################################################################ #
## Formal plots ----

## List to hold all plotting information!
 plot.information <- lapply(est$CA, function(w) {
    x <- colnames(w[[1]])
    
    names(x)[x == "Neurons"] <- "5Neu"
    names(x)[grep("AST|Ast", x)] <- "1Ast"
    names(x)[grep("^Ol", x)] <- "2Oli"
    names(x)[grep("^OP", x)] <- "2OPC"
    names(x)[grep("Exc", x)] <- "3Neu-Exc"
    names(x)[grep("Inh", x)] <- "4Neu-Inh"
    
    
    x <- x[order(names(x))]
    return(x)
  })

# plot.information$universal.order <- c("1Ast", "2End", "3Mic", "4Oli", "5Neu-Exc", "6Neu-Inh", "7Neu-NRGN", "8Neu")
  
  plot.information$universal.order <- c("1Ast", "2Oli", "2OPC", "3Neu-Exc", "4Neu-Inh", "5Neu")
  plot.information$universal.palette <- c("#009392","#f6edbd", "#e9e29c", "#eeb479", "#e88471", "#e88471")  
  names(plot.information$universal.palette) <- plot.information$universal.order

  # note: this convoluted set-up allows the same colour scheme as in the Vel analyses

## Plot barplot
  plot.formal.barplot <- function(merge, add.enrich.algs = TRUE, ncol = 1, substring = TRUE) {
    dat <- sapply(est[[merge]], function(x) { diag(cor(x, true[[merge]][,colnames(x)])) })
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
      dat["OPCs", "Blender"] <- cor(enrich$Blender$OPCs, true$merge1$OPCs)
      
    }
    
    dat$Celltype <- rownames(dat)
    dat <- melt(dat)
    
    colnames(dat)[2] <- "Algorithm"
    # levels(dat$Algorithm)[grep("music", levels(dat$Algorithm))] <- "Music"
    dat$Algorithm <- factor(dat$Algorithm, levels = rev(c("CIB", "DRS", "dtangle", "Music", "Blender", "xCell")))
    levels(dat$Algorithm) <- c("xCell", "Blender", "MuSiC", "dtangle", "DeconRNASeq", "CIBERSORT")
    
    # recolour...
    m <- match(dat$Celltype, plot.information[[merge]])
    dat$Class <- names(plot.information[[merge]][m])
    
    # ...and reorder the x-axis
    dat$Celltype <- factor(dat$Celltype, levels = plot.information[[merge]])
    
    # now plot!
    ggplot(dat, aes(x = Celltype, y = value, fill = Class)) +
            geom_col(position = "dodge", width = 0.7) +

            theme_bw() +
            facet_wrap(~Algorithm, ncol = ncol, strip.position = "top") +
            scale_y_continuous(limits = c(-0.06, 1), expand = c(0,0)) +
            labs(y = "Estimated vs. True (r)") +
            scale_fill_manual(values = plot.information$universal.palette) +
            geom_hline(yintercept = 0.8, linetype = 2) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), panel.border = element_blank(), axis.line = element_line(),
                  axis.title.x = element_blank(), legend.position = "none")
  
  }
  
  pdf(file = "Final Barplot - HCA, merge1.pdf", height = 9, width = 2)
  plot.formal.barplot("merge1", add.enrich.algs = TRUE) #+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + NoLegend()
  dev.off()
  
  pdf(file = "Final Barplot - HCA, merge2.pdf", height = 6.38, width = 2)
  plot.formal.barplot("merge2", add.enrich.algs = FALSE) + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) 
  dev.off()
  
  # pdf(file = "Final Barplot (Tall) - Vel, merge3.pdf", height = 3, width = 8) # note, different structure as it's for the supplements!
  # plot.formal.barplot2("Vel", "merge3", add.enrich.algs = FALSE, substring = FALSE, ncol = 4) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + NoLegend()
  # dev.off()
  
  pdf(file = "Final Barplot - HCA, All Clusters.pdf", height = 7.07, width = 3.4)
  plot.formal.barplot("orig.celltype", add.enrich.algs = FALSE, substring = FALSE)  + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) 
  dev.off()
  
  
  
## Plot nmae  
  # calculate
  errors <- list()
  for(k in sigs) {
    errors[[k]] <- list()
    for (j in names(est[[k]])) {
      errors[[k]][[j]] <- lapply(est[[k]][[j]], function(x) {
        return(t(write.stats(t = true[[j]], e = x, "Alg", TRUE)))
      })
    }  
  }

  # old function, making barplots
  # plot.formal.nmae <- function(sig = "HCA", merge, nrow = 1) {
  #   dat <- melt(stats[[merge]])
  #   dat <- dat[grep("nmae", dat$Var2),]
  #   colnames(dat) <- c("Celltype", "Error", "value", "Algorithm")
  #   dat$Algorithm <- factor(dat$Algorithm)
  #   # levels(dat$Algorithm)[4] <- "Music"
  #   
  #       # recolour...
  #   m <- match(dat$Celltype, plot.information[[merge]])
  #   dat$Class <- names(plot.information[[merge]][m])
  #   
  #   # ...and reorder the x-axis
  #   dat$Celltype <- factor(dat$Celltype, levels = plot.information[[merge]])
  #   
  #   ggplot(dat, aes(x = Celltype, y = value, fill = Class)) +
  #     geom_col(position = "dodge", colour = "black") +
  #     theme_bw() +
  #     facet_wrap(~Algorithm, nrow = nrow) +
  #     # scale_fill_carto_d(palette = palette) +
  #     scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
  #     labs(y = "Normalised Mean\nAbsolute Error") +
  #     scale_fill_manual(values = plot.information$universal.palette) +
  #     geom_hline(yintercept = 1, linetype = 2) +
  #     theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.border = element_blank(), axis.line = element_line(),
  #           axis.title.x = element_blank(), legend.position = "none")
  # }
  
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
      plot.data <- plot.data[which(plot.data$L1 == "CA"),]
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
 
  pdf(file = "Final NMAE - merge1.pdf", height = 3.5, width = 3.4)
  plot.formal.nmae(merge = "merge1") + theme(legend.position = "none")
  dev.off()
  
  pdf(file = "Final NMAE - merge2.pdf", height = 3.5, width = 3.9)
  plot.formal.nmae(merge = "merge2") + theme(legend.position = "none")
  dev.off()

  
  pdf(file = "Final NMAE - All Clusters.pdf", height = 4, width = 8)
  plot.formal.nmae(merge = "orig.celltype") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  dev.off()

## A formal scatterplot...
   plot.formal.scatter <- function(e = NA, sig = "CA", t = NA, alg, merge, autopath = TRUE, nrow = NULL, calc.error = TRUE) {
    if (autopath) {
      e <- est[[sig]][[merge]][[alg]]
      t <- true[[merge]]
    }
    
    ct <- colnames(e[which(colnames(e) %in% colnames(t))])
    t <- t[,ct]
    e <- e[,ct]
    
    # get correlations!
    cor <- diag(cor(e, t))
    cor <- round(cor, 2)
    
     # get nmae
    if (calc.error) {
      err <- list()
      for (j in colnames(e)) { err[[j]] <- nmae(t[,j], e[,j]) }
      err <- do.call("c", err)
      err <- signif(err, 2)  
      err <- paste0("; nmae=", err)
    } else {
      err <- ""
    }
    
    # reformat for plotting
    e <- melt(e)
    t <- melt(t)
    
    plot.data <- data.frame(Celltype = e$variable,
                            Estimate = e$value,
                            True = t$value)
  
    # add colouring information
        # recolour...
    m <- match(plot.data$Celltype, plot.information[[merge]])
    plot.data$Class <- names(plot.information[[merge]][m])
    
    # add correlation information
    levels(plot.data$Celltype) <- paste0(levels(plot.data$Celltype), "\nr=", cor, err)
    
    # plot
    ggplot(plot.data, aes(x = True, y = Estimate, fill = Class)) +
      geom_point(colour = "black", shape = 21) +
      facet_wrap(~Celltype, scales = "free", nrow = nrow) +
      theme_bw() +
      geom_smooth(method = "lm", se = FALSE, colour = "black") +
      scale_fill_manual(values = plot.information$universal.palette) +
      scale_x_continuous(n.breaks = 3) +
      scale_y_continuous(n.breaks = 3) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, size = 6), strip.background = invis,
            strip.text = element_text(size = 6), axis.text.y = element_text(size = 6), panel.grid = invis) +
      labs(y = paste0("Estimated Proportion (", alg, "+", sig, ")"), y = "True Proportion")
  }
  
  # plot all algorithms using orig.clusters.average 
  pdf(file = "Final Scatterplot - All Celltypes.pdf", height = 5.5, width = 8)
  for (j in c("CIB", "DRS", "DTA", "MUS")) {
    print(plot.formal.scatter(alg = j, merge = "orig.celltype", nrow = 3) + labs(y = paste0(j, " Proportion")))          
  }
  dev.off()
  
  # plot all algorithms (including enrichment) using merge1
  pdf(file = "Final Scatterplot - merge1.pdf", height = 1.8, width = 6)
  for (j in c("CIB", "DRS", "DTA", "MUS")) {
    print(plot.formal.scatter(alg = j, merge = "merge1", nrow = 1) + labs(y = paste0(j, " Proportion")))     
  }
  dev.off()
  
  pdf(file = "Final Scatterplot - merge1 (xCell).pdf", height = 1.8, width = 3)
  plot.formal.scatter(e = enrich$xCell, t = true$merge1, sig = "", alg = "xCell", merge = "merge1", nrow = 1, autopath = FALSE, calc.error = FALSE) + labs (y = "xCell Enrichment")
  dev.off()
  
  pdf(file = "Final Scatterplot - merge1 (Blender).pdf", height = 1.8, width = 6)
  temp <- enrich$Blender; colnames(temp)[4] <- "Endothelial"
  plot.formal.scatter(e = temp, t = true$merge1, sig = "", alg = "Blender", merge = "merge1", nrow = 1, autopath = FALSE, calc.error = FALSE) + labs (y = "Blender Enrichment")
  dev.off()
  
  
  # plot all algorithms using merge2
  pdf(file = "Final Scatterplot - merge2.pdf", height = 2, width = 8)
  for (j in c("CIB", "DRS", "DTA", "MUS")) {
    print(plot.formal.scatter(alg = j, merge = "merge2", nrow = 1) + labs(y = paste0(j, " Proportion")))       
  }
  dev.off()
  
  
## Formal exploration of factors affecting performance
  ## First, quantify accuracy
    plot.data <- sapply(est$CA$orig.celltype, function(x) {
      x <- x[,colnames(true$orig.celltype)]
      diag(cor(x, true$orig.celltype)) 
    })
    
    
    plot.data <- plot.data > 0.8 # binary colour for accuracy
    plot.data <- as.data.frame(plot.data)  
    
    
  ## Add celltype information
    m <- match(rownames(plot.data), plot.information$orig.celltype)
    plot.data$Class <- names(plot.information$orig.celltype)[m]
    plot.data$Class <- substr(plot.data$Class, 2, 100)
    plot.data$Class <- gsub("Neu-", "", plot.data$Class)
    plot.data$Class[18] <- "OPC"
    
  ## Add abundance
    plot.data$Abundance <- colMeans(true$orig.celltype)
  
  ## Add each celltype maximum non-self correlation
    cor <- cor(sigsSNME$CA$orig.celltype, method = "s")
    plot.data$Sim <- apply(cor, 1, function(x) max(x[which(x != max(x))]))
    
    plot.data$SimJitter <- plot.data$Sim + 0.008
    plot.data$Label <- letters[1:nrow(plot.data)]
    
    leg <- data.frame(Label = plot.data$Label, Celltype = rownames(plot.data))
    
  ## Plot
    plot.data <- melt(plot.data, id.vars = c("Class", "Sim", "SimJitter", "Abundance", "Label"))
    plot.data$variable <- as.character(plot.data$variable)
    
    pdf(file = "Celltype Similarity vs. Abundance.pdf", height = 6, width = 8)
    ggplot(plot.data, aes(x = (Abundance*100), y = Sim, colour = value, label = Label)) +
      geom_text(size = 4 ) +
      # geom_point(size = 1) +
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
  complete.true <- list(Am1 = snme.gradients$A$meta.prop,
                        Am2 = snme.gradients$A$meta.prop,
                        Bm1 = snme.gradients$B$meta.prop,
                        Bm2 = snme.gradients$B$meta.prop)
    
    complete.true$Am1$Neurons <- complete.true$Am1$Neurons_Exc + complete.true$Am1$Neurons_Inh
    complete.true$Am1 <- complete.true$Am1[,-grep("_", colnames(complete.true$Am1))]
    
    complete.true$Bm1$Neurons <- complete.true$Bm1$Neurons_Exc + complete.true$Bm1$Neurons_Inh
    complete.true$Bm1 <- complete.true$Bm1[,-grep("_", colnames(complete.true$Bm1))]
    
   
  
## On Linseed
  ## Heatmap function
    linseed.heatmap <- function(e, t, filename) {
      plot.data <- cor(e, t, method = "p")
      
      plot.data <- as.data.frame(plot.data)
      plot.data$Linseed <- LETTERS[1:nrow(plot.data)]
      plot.data <- melt(plot.data)
      colnames(plot.data)[3] <- "Correlation"
      plot.data$variable <- gsub("Neurons_", "", plot.data$variable)
      
      pdf(file = filename, height = 2.7, width = 2.3)
      print(ggplot(plot.data, aes(x = variable, y = Linseed, fill = Correlation, label = round(Correlation, 2))) +
              geom_tile(colour = "black") +
              geom_text(size = 2.5) +
              theme_bw() +
              scale_fill_carto_c(palette = "ArmyRose", limits = c(-1,1)) +
              scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
              labs(x = "True Cell-type", y = "Linseed Cell-type") +
              NoLegend() +
              theme(panel.border = invis, axis.title.y = invis, axis.ticks = invis, panel.grid = invis, axis.text.x = element_text(angle = 90, hjust = 1)))
      dev.off()
    }  
    
  ## Apply
    linseed.heatmap(complete$linseed.merge2, true$merge2, "Complete Linseed Heatmap, merge2.pdf")
    linseed.heatmap(complete$linseed.merge2.gB, complete.true$Bm2, "Complete Linseed Heatmap, gradientB, merge2.pdf")
    
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
        theme(panel.border = invis, axis.line = element_line()) +
        NoLegend()
      
      if (abline) {
        return(p + geom_abline(slope = 1, intercept = 0, linetype = 2))
      } else {
        return(p)
      }
    }
    

     ## Apply
    pdf(file = "Complete Linseed Scatter.pdf", height = 2, width = 2)
    complete.scatter(true$merge2$Astrocytes, complete$linseed.merge2$`Cell type 2`, ct = "Astrocytes", ylab = "Linseed Celltype B" , colour = plot.information$universal.palette["1Ast"])
    complete.scatter(true$merge2$Oligodendrocytes, complete$linseed.merge2$`Cell type 3`, ct = "Oligodendrocytes", ylab = "Linseed Celltype C", colour = plot.information$universal.palette["2Oli"])
    dev.off()
    
    pdf(file = "Complete Linseed Scatter, gB.pdf", height = 2, width = 2)
    complete.scatter(complete.true$Bm2$Astrocytes, complete$linseed.merge2.gB$`Cell type 1`, ct = "Astrocytes", ylab = "Linseed Celltype A" , colour = plot.information$universal.palette["1Ast"])
    complete.scatter(complete.true$Bm2$OPC, complete$linseed.merge2.gB$`Cell type 4`, ct = "OPC", ylab = "Linseed Celltype D" , colour = plot.information$universal.palette["2OPC"])
    complete.scatter(complete.true$Bm2$Neurons_Exc, complete$linseed.merge2.gB$`Cell type 2`, ct = "Excitatory Neurons", ylab = "Linseed Celltype B" , colour = plot.information$universal.palette["3Neu-Exc"])
    complete.scatter(complete.true$Bm2$Oligodendrocytes, complete$linseed.merge2.gB$`Cell type 3`, ct = "Oligodendrocytes", ylab = "Linseed Celltype C" , colour = plot.information$universal.palette["2Oli"])
    complete.scatter(complete.true$Bm2$Neurons_Inh, complete$linseed.merge2.gB$`Cell type 5`, ct = "Inhibitory Neurons", ylab = "Linseed Celltype E" , colour = plot.information$universal.palette["4Neu-Inh"])
    
    dev.off()
    
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
      
      pdf(file = filename, height = 2.5, width = 2.3)
      print(ggplot(dat, aes(x = Celltype, y = Cor, fill = Class, colour = Class)) +
              geom_col() +
              theme_bw() +
              scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
              scale_fill_manual(values = plot.information$universal.palette) +
              scale_colour_manual(values = plot.information$universal.palette) +
              scale_y_continuous(limits = c(NA, 1)) +
              labs(y = "Pearson Correlation") +
              geom_hline(yintercept = 0.8, linetype = 2) +
              theme(panel.border = invis, axis.line = element_line(), legend.position = "none"))
      dev.off()
    }
    
    
    coex.barplot(e = complete$coex.merge2, u = true$merge2, t = "merge2", filename = "Complete Coex Barplot, merge2.pdf")
    coex.barplot(e = complete$coex.merge2.gB, u = complete.true$Bm2, t = "merge2", filename = "Complete Coex Barplot, gradientB, merge2.pdf")
 
  ## Heatmap
    coex.heatmap <- function(e, u, filename) {
      dat <- diag(cor(e, u[,colnames(e)]))
      
      absent <- colnames(u)[-which(colnames(u) %in% colnames(e))]
      zeros <- rep(0, length(absent))
      names(zeros) <- absent
      dat <- c(dat, zeros)
      
      dat <- data.frame(Cor = dat, Celltype = names(dat), y = "")
      
      dat$Celltype <- gsub("Neurons_", "", dat$Celltype)
      
      pdf(file = filename, height = 1, width = 2.2)
      print(ggplot(dat, aes(x = Celltype, y = y, fill = Cor, label = round(Cor, 2))) +
              geom_tile(colour = "black") +
              geom_text(size = 2.5) +
              theme_bw() +
              scale_fill_carto_c(palette = "ArmyRose", limits = c(-1,1)) +
              scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
              labs(x = "True Cell-type") +
              NoLegend() +
              theme(panel.border = invis, axis.ticks = invis, panel.grid = invis, axis.title.y = invis, axis.text.x = element_text(angle = 90, hjust = 1)))
      dev.off()
    }
    
    coex.heatmap(e = complete$coex.merge2, u = true$merge2, filename = "Complete Coex Heatmap.pdf")
    coex.heatmap(e = complete$coex.merge2.gB, u = complete.true$Bm2, filename = "Complete Coex Heatmap, gradientB.pdf")
    
    
    
  ## Scatterplot
      pdf(file = "Complete Coex Scatterplot.pdf", height = 2, width = 2)
      complete.scatter(true$merge2$Astrocytes, complete$coex.merge2$Astrocytes, ct = "Astrocytes", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["1Ast"])
      complete.scatter(true$merge2$Neurons_Inh, complete$coex.merge2$Neurons_Inh, ct = "Inhibitory Neurons", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["4Neu-Inh"])
      complete.scatter(true$merge2$Oligodendrocytes, complete$coex.merge2$Oligodendrocytes, ct = "Oligodendrocytes", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["2Oli"])
      dev.off()
      
      
      pdf(file = "Complete Coex Scatterplot, gradientB.pdf", height = 2, width = 2)
      complete.scatter(complete.true$Bm2$Astrocytes, complete$coex.merge2.gB$Astrocytes, ct = "Astrocytes", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["1Ast"])
      complete.scatter(complete.true$Bm2$Neurons_Inh, complete$coex.merge2.gB$Neurons_Inh, ct = "Inhibitory Neurons", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["4Neu-Inh"])
      complete.scatter(complete.true$Bm2$Oligodendrocytes, complete$coex.merge2.gB$Oligodendrocytes, ct = "Oligodendrocytes", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["4Oli"])
      complete.scatter(complete.true$Bm2$OPC, complete$coex.merge2.gB$OPC, ct = "OPCs", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["4OPC"])
      dev.off()
      
      
## Aside: properties of mixtures
  ## Abundance
      pdf(file = "Complete QC Abundance.pdf", height = 1.8, width = 2.5)
      dat <- melt(true$merge2)
      dat$variable <- gsub("Neurons_", "", dat$variable)
      ggplot(dat, aes(x = variable, y = value)) +
        geom_violin(scale = "width", width = 0.8, draw_quantiles = 0.5, colour = "black", fill = "darkorange1") +
        theme_bw() +
        scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
        scale_y_continuous(limits = c(0, 0.85), expand = c(0,0)) +
        theme(panel.grid = invis, axis.title.x = invis, panel.border = invis, axis.line = element_line(), legend.position = "none") +
        labs(y = "Simulated Abundance")
      
      dat <- melt(complete.true$Bm2)
      dat$variable <- gsub("Neurons_", "", dat$variable)
      ggplot(dat, aes(x = variable, y = value)) +
        geom_violin(scale = "width", width = 0.8, draw_quantiles = 0.5, colour = "black", fill = "darkorange1") +
        theme_bw() +
        scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
        scale_y_continuous(limits = c(0, 0.85), expand = c(0,0)) +
        theme(panel.grid = invis, axis.title.x = invis, panel.border = invis, axis.line = element_line(), legend.position = "none") +
        labs(y = "Simulated Abundance")
      
      dev.off()
    
    
    
    
  ## CorMx
    # random
      qc.corMx <- function(t, filename, flip = FALSE) {
      plot.data <- cor(t, method = "p")
      plot.data <- as.data.frame(plot.data)
      rownames(plot.data) <- gsub("Neurons_", "", rownames(plot.data))

      plot.data$Ct <- colnames(plot.data) <- substr(rownames(plot.data), 1, 3)
      
      plot.data <- melt(plot.data, na.rm = TRUE)
      plot.data$Ct <- factor(plot.data$Ct)
      plot.data$variable <- factor(plot.data$variable, levels = levels(plot.data$Ct))
      
      colnames(plot.data)[3] <- "Correlation"
      
  
      a <- ggplot(plot.data, aes(x = variable, y = Ct, fill = Correlation, label = round(Correlation, 2))) +
              geom_tile(colour = "black") +
              geom_text(size = 2.5) +
              theme_bw() +
              # coord_flip() +
              scale_fill_carto_c(palette = "Earth", limits = c(-1,1)) +
              scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
              scale_y_discrete(position = "right") +
              NoLegend() +
              theme(panel.border = invis, axis.ticks = invis, panel.grid = invis, axis.title = invis, axis.text.x = element_text(angle = 90, hjust = 1))
      
      if (flip) {
        pdf(file = filename, height = 3, width = 2.7)
        print(a + theme(axis.text.y = invis()))
        dev.off()
      } else {
        pdf(file = filename, height = 3, width = 2.7)
        print(a)
        dev.off()
      }
    }
    
    
   qc.corMx(true$merge2[,colnames(complete.true$Bm2)], "Complete QC Random Cors.pdf")
   qc.corMx(complete.true$Bm2, "Complete QC Non-Random Cors.pdf")
   
  

################################################################################################################################ #
## Formal plots of deconvolution using mismatched signatures ----
  
  
# ## Scatterplots
#   ## Function
#   plot.formal.mismatch.scatter <- function(sig, alg, merge, autopath = TRUE, no.xlab = TRUE, no.ylab = TRUE) {
#     if (autopath) {
#       e <- est[[sig]][[merge]][[alg]]
#       t <- true[[merge]]
#     }
#     
#     # colnames(e) <- gsub("OPCs", "OPC", colnames(e)) 
#     t <- t[,colnames(e)]
#     
#     # reformat for plotting
#     e <- melt(e)
#     t <- melt(t)
#     
#     plot.data <- data.frame(Celltype = e$variable,
#                             Estimate = e$value,
#                             True = t$value)
#     
#     # add colouring information
#     m <- match(plot.data$Celltype, plot.information[[merge]])
#     plot.data$Class <- names(plot.information[[merge]][m])
#     
#     # shuffle order
#     plot.data <- plot.data[sample(rownames(plot.data), nrow(plot.data), replace = FALSE),]
#     
#     # plot
#     plot <- ggplot(plot.data, aes(x = True*100, y = Estimate*100, colour = Class)) +
#       geom_point(alpha = 1, size = 1) +
#       theme_bw() +
#       scale_colour_manual(values = plot.information$universal.palette) +
#       scale_x_continuous(trans = pseudo_log_trans(sigma = 5), limits = c(0, 120), expand = c(0,0), breaks = c(0, 2, 5, 10, 25, 50, 100)) +
#       scale_y_continuous(trans = pseudo_log_trans(sigma = 5), limits = c(0, 120), expand = c(0,0), breaks = c(0, 2, 5, 10, 25, 50, 100)) +
#       geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
#       theme(legend.position = "none", panel.grid = invis, panel.border = invis,
#             axis.line = element_line(), plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill = "grey80")) +
#       labs(y = paste0("Estimated Proportion"), x = "True Proportion", title = sig)
#     
#     if (no.xlab) plot <- plot + theme(axis.title.x = invis, axis.text.x = invis, axis.ticks.x = invis)
#     if (no.ylab) plot <- plot + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#     
#     return(plot)
#   }
#   
#   ## Apply function
#   for (j in algs) {
#     if (j == "MUS") next
#     pdf(file = paste0("Mismatch Signature - ", j, " merge2.pdf"), height = 5, width = 8)
#     pA <- plot.formal.mismatch.scatter(sig = "CA", alg = j, merge = "merge2", autopath = TRUE, no.ylab = FALSE) + labs(title = "CA (Matched)")
#     pB <- plot.formal.mismatch.scatter(sig = "NG", alg = j, merge = "merge2", autopath = TRUE)
#     pC <- plot.formal.mismatch.scatter(sig = "VL", alg = j, merge = "merge2", autopath = TRUE)
#     pD <- plot.formal.mismatch.scatter(sig = "LK", alg = j, merge = "merge2", autopath = TRUE, no.ylab = FALSE, no.xlab = FALSE)
#     pE <- plot.formal.mismatch.scatter(sig = "TS", alg = j, merge = "merge2", autopath = TRUE, no.xlab = FALSE)
#     pF <- get_legend(pA + theme(legend.position = "right"))
#     # pF <- plot.formal.mismatch.scatter(sig = "Dar", alg = "CIB", merge = "merge1", autopath = TRUE, no.xlab = FALSE) # because it's merge1
#     print(plot_grid(pA, pB, pC, pD, pE, pF, nrow = 2, rel_widths = c(1, 0.8, 0.8), rel_heights = c(1, 1.2)))
#     dev.off()
#     
#     pdf(file = paste0("Mismatch Signature - ", j, " merge1.pdf"), height = 5.5, width = 8)
#     pA <- plot.formal.mismatch.scatter(sig = "CA", alg = j, merge = "merge1", autopath = TRUE, no.ylab = FALSE) + labs(title = "CA (Matched)")
#     pB <- plot.formal.mismatch.scatter(sig = "NG", alg = j, merge = "merge1", autopath = TRUE)
#     pC <- plot.formal.mismatch.scatter(sig = "VL", alg = j, merge = "merge1", autopath = TRUE)
#     pD <- plot.formal.mismatch.scatter(sig = "LK", alg = j, merge = "merge1", autopath = TRUE, no.ylab = FALSE)
#     pE <- plot.formal.mismatch.scatter(sig = "TS", alg = j, merge = "merge1", autopath = TRUE)
#     pF <- plot.formal.mismatch.scatter(sig = "F5", alg = j, merge = "merge1", autopath = TRUE) 
#     pG <- plot.formal.mismatch.scatter(sig = "IP", alg = j, merge = "merge1", autopath = TRUE, no.ylab = FALSE, no.xlab = FALSE)
#     pH <- plot.formal.mismatch.scatter(sig = "MM", alg = j, merge = "merge1", autopath = TRUE, no.xlab = FALSE)
#     pI <- plot.formal.mismatch.scatter(sig = "DM", alg = j, merge = "merge1", autopath = TRUE, no.xlab = FALSE)
#     pJ <- get_legend(pA + theme(legend.position = "bottom"))
#     print(plot_grid(pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ, nrow = 4, ncol = 3, rel_widths = c(1, 0.8, 0.8), rel_heights = c(1, 1, 1.2, 0.25)))
#     dev.off()
#   }
#   
#   
# ## Barplots of error
#   plot.formal.mismatch.nmae <- function(sig, merge, alg, ylim = 2) {
#     dat <- melt(errors[[sig]][[merge]])
#     dat <- dat[grep("nmae", dat$Var2),]
#     colnames(dat) <- c("Celltype", "Error", "value", "Algorithm")
#     dat$Algorithm <- factor(dat$Algorithm)
#     levels(dat$Algorithm)[4] <- "MUS"
#     
#     dat <- dat[which(dat$Algorithm == alg),]
#     
#     
#     # recolour...
#     m <- match(dat$Celltype, plot.information[[merge]])
#     dat$Class <- names(plot.information[[merge]][m])
#     
#     # ...and reorder the x-axis
#     dat$Celltype <- factor(dat$Celltype, levels = plot.information[[merge]])
#     levels(dat$Celltype) <- gsub(pattern = "Neurons_", replacement = "", levels(dat$Celltype))
#     levels(dat$Celltype) <- substr(levels(dat$Celltype), 1, 3)
#     
#     # add a plot title
#     # label <- paste0(sig, " (mean=", round(mean(dat$value),2), ")")
#     label <- sig
#     
#     # limits
#     if (max(dat$value) > ylim) {
#       ylim <- max(dat$value) + 1
#     }
#     
#     # plot
#     ggplot(dat, aes(x = Celltype, y = value, fill = Class)) +
#       geom_col(position = "dodge", colour = "black") +
#       theme_bw() +
#       scale_y_continuous(limits = c(0, ylim), expand = c(0,0)) +
#       labs(y = "Normalised Mean\nAbsolute Error", title = label) +
#       scale_fill_manual(values = plot.information$universal.palette) +
#       geom_hline(yintercept = 1, linetype = 2) +
#       theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), panel.border = invis, axis.line = element_line(),
#             axis.title.x = invis, legend.position = "none", plot.title = element_text(hjust = 0.5), title = element_text(size = 8))
#   }
#   
#   for (j in algs) {
#     # setup
#     if (j == "MUS") next
#     if (j == "CIB") {
#       ylim <- 2 
#     } else {
#       ylim <- 5
#     }
#     
#     # run
#     pdf(file = paste0("Mismatch Signature Error - ", j, " merge2.pdf"), height = 2, width = 8)
#     pA <- plot.formal.mismatch.nmae("CA", "merge2", j, ylim = ylim) + labs(title = "CA (Matched)")
#     pB <- plot.formal.mismatch.nmae("NG", "merge2", j, ylim = ylim) + theme(axis.title.y = invis)
#     pC <- plot.formal.mismatch.nmae("VL", "merge2", j, ylim = ylim) + theme(axis.title.y = invis)
#     pD <- plot.formal.mismatch.nmae("LK", "merge2", j, ylim = ylim) + theme(axis.title.y = invis)
#     pE <- plot.formal.mismatch.nmae("TS", "merge2", j, ylim = ylim) + theme(axis.title.y = invis)
#     print(plot_grid(pA, pB, pC, pD, pE, nrow = 1, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)))
#     dev.off()
#     
#     pdf(file = paste0("Mismatch Signature Error - ", j, " merge1.pdf"), height = 3, width = 8)
#     pA <- plot.formal.mismatch.nmae("CA", "merge1", j, ylim = ylim) + labs(title = "CA (Matched)")
#     pB <- plot.formal.mismatch.nmae("NG", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
#     pC <- plot.formal.mismatch.nmae("VL", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
#     pD <- plot.formal.mismatch.nmae("LK", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
#     pE <- plot.formal.mismatch.nmae("TS", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
#     pF <- plot.formal.mismatch.nmae("F5", "merge1", j, ylim = ylim) 
#     pG <- plot.formal.mismatch.nmae("IP", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
#     pH <- plot.formal.mismatch.nmae("DM", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
#     pI <- plot.formal.mismatch.nmae("MM", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
#     pJ <- get_legend(pA + theme(legend.position = "right", legend.title = invis))
#     print(plot_grid(pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ, nrow = 2, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8)))
#     dev.off()
#   }
#   
# 
# 
#   
# ## Barplots of correlation
#   ## Function
#  plot.formal.mismatch.cor <- function(e = NA, t = NA, sig, alg, merge, autopath = TRUE, no.xlab = TRUE, no.ylab = TRUE) {
#    if (autopath) {
#       e <- est[[sig]][[merge]][[alg]]
#       t <- true[[merge]]
#     }
#     
#     # colnames(e) <- gsub("OPCs", "OPC", colnames(e)) 
#     t <- t[,colnames(e)]
#     
#     cor <- diag(cor(e, t, method = "p"))
#     
#     if(anyNA(cor)) cor[which(is.na(cor))] <- 0
#     
#     dat <- data.frame(Celltype = names(cor), value = cor)
#     
#     # recolour...
#     m <- match(dat$Celltype, plot.information[[merge]])
#     dat$Class <- names(plot.information[[merge]][m])
#     
#   # ...and reorder the x-axis
#     dat$Celltype <- factor(dat$Celltype, levels = plot.information[[merge]])
#     levels(dat$Celltype) <- gsub(pattern = "Neurons_", replacement = "", levels(dat$Celltype))
#     levels(dat$Celltype) <- substr(levels(dat$Celltype), 1, 3)
#     
#     # add a plot title
#     # label <- paste0(sig, " (mean=", round(mean(dat$value),2), ")")
#     label <- sig
#     
#     # plot
#     ggplot(dat, aes(x = Celltype, y = value, fill = Class)) +
#       geom_col(position = "dodge", colour = "black") +
#       theme_bw() +
#       scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
#       labs(y = "Pearson Correlation", title = label) +
#       scale_fill_manual(values = plot.information$universal.palette) +
#       geom_hline(yintercept = 0.8, linetype = 2) +
#       theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), panel.border = invis, axis.line = element_line(),
#             axis.title.x = invis, legend.position = "none", plot.title = element_text(hjust = 0.5), title = element_text(size = 8))
#  }
#  
#   ## Apply
#     for (j in algs) {
#       if (j == "MUS") next
#       pdf(file = paste0("Mismatch Signature Correlation - ", j, " merge2.pdf"), height = 2, width = 8)
#       pA <- plot.formal.mismatch.cor(sig = "CA", alg = "CIB", merge = "merge2", autopath = TRUE, no.ylab = FALSE) + labs(title = "CA (Matched)")
#       pB <- plot.formal.mismatch.cor(sig = "NG", alg = "CIB", merge = "merge2", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#       pC <- plot.formal.mismatch.cor(sig = "VL", alg = "CIB", merge = "merge2", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#       pD <- plot.formal.mismatch.cor(sig = "LK", alg = "CIB", merge = "merge2", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#       pE <- plot.formal.mismatch.cor(sig = "TS", alg = "CIB", merge = "merge2", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#       print(plot_grid(pA, pB, pC, pD, pE, nrow = 1, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)))
#       dev.off()
#       
#       pdf(file = paste0("Mismatch Signature Correlation - ", j, " merge1.pdf"), height = 3, width = 8)
#       pA <- plot.formal.mismatch.cor(sig = "CA", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + labs(title = "CA (Matched)")
#       pB <- plot.formal.mismatch.cor(sig = "NG", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#       pC <- plot.formal.mismatch.cor(sig = "VL", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#       pD <- plot.formal.mismatch.cor(sig = "LK", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#       pE <- plot.formal.mismatch.cor(sig = "TS", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#       pF <- plot.formal.mismatch.cor(sig = "F5", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) 
#       pG <- plot.formal.mismatch.cor(sig = "IP", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#       pH <- plot.formal.mismatch.cor(sig = "DM", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#       pI <- plot.formal.mismatch.cor(sig = "MM", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
#       pJ <- get_legend(pA + theme(legend.position = "right", legend.title = invis))
#       print(plot_grid(pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ, nrow = 2, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8)))
#       dev.off()
#     }
  
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
    if (sig == "CA") lab <- paste0(lab, " (Matched)")
    
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
      geom_point(alpha = 1, size = 1) +
      theme_bw() +
      annotate("text", x = annot.x, y = annot.y, label = lab) +
      scale_colour_manual(values = plot.information$universal.palette) +
      scale_x_continuous(trans = pseudo_log_trans(sigma = sigma), limits = c(0, 120), expand = c(0,0), breaks = c(2, 5, 10, 25, 50, 100)) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = sigma), limits = c(0, 120), expand = c(0,0), breaks = c(2, 5, 10, 25, 50, 100)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
      theme(legend.position = "none", panel.grid = invis, panel.border = invis,
            axis.line = element_line(), plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill = "grey80")) +
      labs(y = paste0("Estimated Proportion"), x = "True Proportion")
    
    if (no.xlab) plot <- plot + theme(axis.title.x = invis, axis.text.x = invis, axis.ticks.x = invis)
    if (no.ylab) plot <- plot + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
    
    return(plot)
  }
  
  ## Apply function
  for (j in algs) {
    
    pdf(file = paste0("Mismatch Signature - ", j, " merge2.pdf"), height = 2.5, width = 8)
    pA <- plot.formal.mismatch.scatter(sig = "CA", alg = j, merge = "merge2", no.ylab = FALSE, no.xlab = FALSE)
    pB <- plot.formal.mismatch.scatter(sig = "NG", alg = j, merge = "merge2", no.xlab = FALSE)
    pC <- plot.formal.mismatch.scatter(sig = "VL", alg = j, merge = "merge2", no.xlab = FALSE)
    pD <- plot.formal.mismatch.scatter(sig = "LK", alg = j, merge = "merge2", no.xlab = FALSE)
    pE <- plot.formal.mismatch.scatter(sig = "TS", alg = j, merge = "merge2", no.xlab = FALSE)
    pF <- get_legend(pA + theme(legend.position = "right"))
    print(plot_grid(pA, pB, pC, pD, pE, nrow = 1, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8)))
    dev.off()
    
    if (j == "MUS") next
    
    pdf(file = paste0("Mismatch Signature - ", j, " merge1.pdf"), height = 2.8, width = 8)
    pA <- plot.formal.mismatch.scatter(sig = "CA", alg = j, merge = "merge1", no.ylab = FALSE)
    pB <- plot.formal.mismatch.scatter(sig = "NG", alg = j, merge = "merge1")
    pC <- plot.formal.mismatch.scatter(sig = "VL", alg = j, merge = "merge1")
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
  
  
  ## Separately plot MUS merge1
   pdf(file = paste0("Mismatch Signature - MUS merge1.pdf"), height = 2.8, width = 5)
    pA <- plot.formal.mismatch.scatter(sig = "CA", alg = "MUS", merge = "merge1", no.ylab = FALSE)
    pB <- plot.formal.mismatch.scatter(sig = "NG", alg = "MUS", merge = "merge1")
    pC <- plot.formal.mismatch.scatter(sig = "VL", alg = "MUS", merge = "merge1")
    pD <- plot.formal.mismatch.scatter(sig = "LK", alg = "MUS", merge = "merge1", no.xlab = FALSE, no.ylab = FALSE) + labs(y = "")
    pE <- plot.formal.mismatch.scatter(sig = "TS", alg = "MUS", merge = "merge1", no.xlab = FALSE)
    pI <- plot.formal.mismatch.scatter(sig = "DM", alg = "MUS", merge = "merge1", no.xlab = FALSE)
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
    
## Goodness of fit across signatures!
  # here, I am calculating gof for exc/inh subtypes (where possible), for CIB-based deconvolutions
  gof <- list()
    
  # function
  get.gof <- function(e, sig) {
    sig <- sig[,colnames(e)]
    output <- write.gof.v2(estimatedComp = e, measuredExp = snme$mixture$cpm, signatureUsed = sig)
    output <- output$r
  }
  
  # calculate
  gof$VL <- get.gof(est$VL$merge2$CIB, sigsSNME$VL)
  gof$CA <- get.gof(est$CA$merge2$CIB, sigsSNME$CA$merge2)
  gof$TS <- get.gof(est$TS$merge2$CIB, sigsSNME$TS)
  gof$LK <- get.gof(est$LK$merge2$CIB, sigsSNME$LK)
  gof$NG <- get.gof(est$NG$merge2$CIB, sigsSNME$NG)
  gof$F5 <- get.gof(est$F5$merge1$CIB, sigsSNME$F5)
  gof$IP <- get.gof(est$IP$merge1$CIB, sigsSNME$IP)
  gof$MM <- get.gof(est$MM$merge1$CIB, sigsSNME$MM)
  gof$DM <- get.gof(est$DM$merge1$CIB, sigsSNME$DM)
  
  # generate plot data
  plot.data <- do.call("cbind", gof)
  plot.data <- melt(plot.data)
  order <- names(gof)[rev(order(sapply(gof, mean)))]
  plot.data$Var2 <- factor(plot.data$Var2, levels = order)
  
  
  
  # and plot! 
  pdf(file = "GoF.pdf", height = 3, width = 4)
  ggplot(plot.data, aes(x = Var2, y = value)) +
    geom_point(position = position_jitter(), alpha = 0.2, size = 0.5) +
    theme_bw() +
    theme(panel.grid = invis, panel.border = invis, axis.line = element_line()) +
    geom_hline(yintercept = c(0.5, 0.7), linetype = 2) +
    labs(y = "Goodness of Fit (r)", x = "Signature")
  dev.off()