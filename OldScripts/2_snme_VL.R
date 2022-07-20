################################################################################################################################ #
## Setup ----

## Start!
rm(list = ls())
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/snme/VL/")

## Load
load("../../../../Data/Preprocessed/Signatures - SNME.rda")
load("../../../../Data/Preprocessed/Signatures - Brain (MuSiC).rda")
load("../../../../Data/Preprocessed/snme.rda")
source("../../../../Scripts/Fun_Composition.R")
source("../../../../Scripts/Fun_Preprocessing.R")
library(rcartocolor)

## Unload the CA-based SNME
snme <- snme$Vel

## Parameters
algs <- c("DRS", "DTA", "CIB", "MUS")
sigs <- names(sigsSNME)
m1.order <- sigs
m2.order <- c("VL", "NG", "CA", "LK", "TS")

## Modify signatures
  ## Remove $VL from sigsMuSiC, as it's in sigsSNME
    sigsMuSiC <- sigsMuSiC[-1]

  ## Add dummy columns to sigsMuSiC DM
    sigsMuSiC$DM$meta$brain.ct2 <- sigsMuSiC$DM$meta$brain.ct <- sigsMuSiC$DM$meta$orig.celltype

  ## A quick change to Vel$merge2 (removing NRGN)
    sigsSNME$VL$merge2 <- sigsSNME$VL$merge2[,-3]
    sigsSNME$VL$Meta$merge2[grep("NRGN", sigsSNME$VL$Meta$orig.celltype)] <- "Drop"
  
  ## Rename celltypes
    sigsSNME[-c(1:2)] <- lapply(sigsSNME[-c(1:2)], function(x) {
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
    
  ## And modify the formatting of CA, as well as add Microglia
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

  ## Remove End from CA and LK
    # music
    remove <- grep("End", sigsMuSiC$CA$meta$brain.ct2)
    sigsMuSiC$CA$counts <- sigsMuSiC$CA$counts[,-remove]
    sigsMuSiC$CA$meta <- sigsMuSiC$CA$meta[-remove,]    
    
    remove <- grep("End", sigsMuSiC$LK$meta$brain.ct2)
    sigsMuSiC$LK$counts <- sigsMuSiC$LK$counts[,-remove]
    sigsMuSiC$LK$meta <- sigsMuSiC$LK$meta[-remove,]    
    
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
true$orig.celltype <- as.data.frame(snme$true$orig.celltype)    
 

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

## Merge 3
# true$merge3 <- data.frame(Astrocytes = x$`AST-FB` + x$`AST-PP`,
#                           Endothelial = x$Endothelial,
#                           Microglia = x$Microglia,
#                           Neurons_Exc1 = rowSums(x[,c("L2/3", "L4", "L5/6-CC", "Neu-mat")]),
#                           Neurons_Exc2 = x$`L5/6`,
#                           Neurons_Inh1 = rowSums(x[,c("IN-SST", "IN-VIP")]),
#                           Neurons_Inh2 = rowSums(x[,c("IN-PV","IN-SV2C")]),
#                           Neurons_NRGN = rowSums(x[14:15]),
#                           Oligodendrocytes = x$Oligodendrocytes,
#                           OPC = x$OPC)

## Drop 1: removes all Inh Neurons
true$drop1 <- true$merge2[,colnames(sigsSNME$VL$drop1)[-2]]

## Drop 2: removes L4 neurons from Exc Neurons
true$drop2 <- true$merge2

## Drop 3: removes L4 neurons without remerging
true$drop3 <- true$orig.celltype[,colnames(sigsSNME$VL$drop3)] 

## Drop 4: removes Oli
true$drop4 <- true$merge1[,colnames(sigsSNME$VL$drop4)]

################################################################################################################################ #
## Quick QC ----


## This plot is for Fig 1A; it is a summary of the dataset!
  dat <- data.frame(Celltype = plot.information$merge2, Label = names(plot.information$merge2), n = c(2229, 523, 450, 4721, 2677, 8601, 4328), Dummy = "X")
  dat <- rbind(dat, c("Neurons_NRGN", "7Neu-NRGN", 1117, "X"))
  dat$n <- as.numeric(dat$n)
  dat$Celltype <- paste0(dat$Celltype, " (",  dat$n, ")")
  # dat$Celltype <- paste0(dat$n, " ", dat$Celltype)
  
  dat$Celltype[1] <- paste0(dat$Celltype[1], " (2 subtypes)")
  dat$Celltype[6] <- paste0(dat$Celltype[6], " (5 subtypes)")
  dat$Celltype[7] <- paste0(dat$Celltype[7], " (4 subtypes)")
  dat$Celltype[8] <- paste0(dat$Celltype[8], " (2 subtypes)")
  
  pdf(file = "QC - Number of Nuclei.pdf", height = 4, width = 2.6)
  ggplot(dat, aes(x = Dummy, y = n, fill = Label)) +
    geom_col(colour = "black") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = plot.information$universal.palette, labels = dat$Celltype) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    labs(y = "Number of Nuclei") +
    guides(fill = guide_legend(ncol = 1)) +
    theme(panel.grid = element_blank(), panel.border = element_blank(), axis.ticks  = element_blank(), axis.title = element_blank(), 
          axis.text = element_blank(), legend.title = element_blank(), legend.position = "bottom")
  dev.off()  
    
################################################################################################################################ #
## Deconvolve ----
    
## Setup
  est <- list()
  enrich <- list()
  complete <- list()
  
## Loop for partial deconvolution algorithms
  ## For VL, run separately
    est$VL <- list()
    
    for(k in names(sigsSNME$VL[-c(1:2)])) {
      
      # setup  
      est$VL[[k]] <- list()

      print(paste0(k, "_DRS"))
      est$VL[[k]]$DRS <- run.DRS(snme$mixture$cpm, sigsSNME$VL[[k]])

      print(paste0(k, "_DTA"))
      est$VL[[k]]$DTA <- run.DTA(snme$mixture$cpm, sigsSNME$VL[[k]], alg = "diff", q = 0.01)

      print(paste0(k, "_CIB"))
      est$VL[[k]]$CIB <- run.CIB(from.file = FALSE,
                                   sigObject = sigsSNME$VL[[k]],
                                   mixString = "snme_Vel.txt")
      
      print(paste0(k, "_MUS"))
      if (substr(k, 1, 4) == "drop" | (j == "Vel" & k == "merge2")) {
       drop <- TRUE 
       print("Dropping ct in MUS")
      } else {
        drop <- FALSE
      }
      
      est$VL[[k]]$MUS <- run.music(mixture = snme$mixture$counts, use.meta.column = k, drop.ct = drop,
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
      e$merge1$DRS <- run.DRS(snme$mixture$cpm, s1)
      e$merge1$DTA <- run.DTA(snme$mixture$cpm, s1, alg = "diff", q = 0.01)
      print("Starting CIB merge1")
      e$merge1$CIB <- run.CIB(from.file = FALSE,
                              sigObject = s1,
                              mixString = "snme_Vel.txt")
      
      if (sig %in% names(sigsMuSiC)) {
        e <- run.music(mixture = snme$mixture$counts, 
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
                                mixString = "snme_Vel.txt")
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
    
  ## Note that MuSiC automatically dropped DM's OPCs, as it's only present in a single individual
    est$DM$merge1$MUS$OPC <- 0
    
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
  enrich$xCell <- run.xCell(snme$mixture$cpm)
  enrich$Blender <- run.Blender(snme$mixture$cpm, OPCs = TRUE)
  
## Save
  save(est, enrich, file = "Final SNME VL Composition.rda")
  
################################################################################################################################ #
## Complete deconvolution ----

  
## First, load the gradient mixtures
  # load("../../../../Data/Preprocessed/snme_gradient.rda")
  load("../../../../Data/Preprocessed/VL_snme_gradient.rda")

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
                                    signature = sigsSNME$VL$merge1, 
                                    only.threshold = FALSE, 
                                    sft = "auto", 
                                    output.all = TRUE)
  
  ## Extract information about merge1 cell-type abundance
    complete$coex.merge1 <- analyse.coex(coex.standardOutput, "merge1")
    complete$coex.merge2 <- analyse.coex(coex.standardOutput, "merge2")
    
## Coex: on the gradient mixtures
  ## Run
    coex.gradientOutput.A <- run.coex(mixture = snme.gradients$A$confound.p,
                                      signature = sigsSNME$VL$merge1,
                                      only.threshold = FALSE,
                                      sft = "auto", 
                                      output.all = TRUE)
    
    coex.gradientOutput.B <- run.coex(mixture = snme.gradients$B$confound.p,
                                      signature = sigsSNME$VL$merge1,
                                      only.threshold = FALSE,
                                      sft = "auto", 
                                      output.all = TRUE)

  ## Analyse
    complete$coex.merge1.gA <- analyse.coex(coex.gradientOutput.A, "merge1")
    complete$coex.merge2.gA <- analyse.coex(coex.gradientOutput.A, "merge2")
    
    complete$coex.merge1.gB <- analyse.coex(coex.gradientOutput.B, "merge1")
    complete$coex.merge2.gB <- analyse.coex(coex.gradientOutput.B, "merge2")
  
## Save
  save(complete, file = "Complete Deconvolution Estimates.rda")

 

################################################################################################################################ #
## Formal plots of deconvolution using Vel signatures ----

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
  
  # function: old version of the plot, which was a barplot
  # plot.formal.nmae <- function(sig, merge, nrow = 1) {
  #   dat <- melt(errors[[sig]][[merge]])
  #   dat <- dat[grep("nmae", dat$Var2),]
  #   colnames(dat) <- c("Celltype", "Error", "value", "Algorithm")
  #   dat$Algorithm <- factor(dat$Algorithm)
  #   
  #   # recolour...
  #   m <- match(dat$Celltype, plot.information[[merge]])
  #   dat$Class <- names(plot.information[[merge]][m])
  #   
  #   # ...and reorder the x-axis
  #   dat$Celltype <- factor(dat$Celltype, levels = plot.information[[merge]])
  #   
  #   # plot
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
 
  pdf(file = "Final NMAE - Vel, merge1.pdf", height = 3.5, width = 3.5)
  plot.formal.nmae("merge1") + theme(legend.position = "none")
  dev.off()
  
  pdf(file = "Final NMAE - Vel, merge2.pdf", height = 3.5, width = 4)
  plot.formal.nmae("merge2") + theme(legend.position = "none", axis.title.y = element_blank())
  dev.off()

  pdf(file = "Final NMAE - Vel, All Clusters.pdf", height = 3.5, width = 8)
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
  
  pdf(file = "Final Barplot - Vel, merge1.pdf", height = 6.2, width = 2.3)
  plot.formal.barplot("VL", "merge1", add.enrich.algs = FALSE, lower.ylim = y) 
  dev.off()
  
  pdf(file = "Final Barplot - Vel, merge2.pdf", height = 6.16, width = 2.3)
  plot.formal.barplot("VL", "merge2", add.enrich.algs = FALSE, lower.ylim = y) + labs(y = "") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  dev.off()
  
  pdf(file = "Final Barplot - Vel, All Clusters.pdf", height = 6.18, width = 3.4)
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
  
  
## A formal scatterplot...
   plot.formal.scatter <- function(e = NA, t = NA, sig, alg, merge, autopath = TRUE, nrow = NULL, calc.error = TRUE) {
    if (autopath) {
      e <- est[[sig]][[merge]][[alg]]
      t <- true[[merge]]
    }
  
    t <- t[,colnames(e)]
    
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
      # scale_x_continuous(n.breaks = 3) +
      # scale_y_continuous(n.breaks = 3) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, size = 6), strip.background = invis,
            strip.text = element_text(size = 6), axis.text.y = element_text(size = 6), panel.grid = invis) +
      labs(y = paste0("Estimated Proportion (", alg, "+", sig, ")"), y = "True Proportion")
  }
  
  # plot all algorithms using orig.clusters.average 
  pdf(file = "Final Scatterplot - All Celltypes.pdf", height = 5.5, width = 8)
  for (j in algs) {
    print(plot.formal.scatter(sig = "Vel", alg = j, merge = "orig.celltype", nrow = 3) + labs(y = paste0(j, " Proportion")))          
  }
  dev.off()
  
  # plot all algorithms (including enrichment) using merge1
  pdf(file = "Final Scatterplot - merge1.pdf", height = 1.8, width = 8)
  for (j in algs) {
    print(plot.formal.scatter(sig = "Vel", alg = j, merge = "merge1", nrow = 1) + theme(axis.title.x = invis) + labs(y = paste0(j, " Proportion")))     
  }
  dev.off()
  
  pdf(file = "Final Scatterplot - merge1 (xCell).pdf", height = 2, width = 3.5)
  plot.formal.scatter(e = enrich$xCell, t = true$merge1, sig = "", alg = "xCell", merge = "merge1", nrow = 1, autopath = FALSE, calc.err = FALSE) +
    labs(y = "xCell Score")
  dev.off()
  
  pdf(file = "Final Scatterplot - merge1 (Blender).pdf", height = 1.8, width = 8)
  temp <- enrich$Blender; colnames(temp)[4] <- "Endothelial"
  plot.formal.scatter(e = temp, t = true$merge1, sig = "", alg = "Blender", merge = "merge1", nrow = 1, autopath = FALSE, calc.err = FALSE) +
    theme(axis.title.x = invis) + labs(y = "Blender Estimate")
  dev.off()
  
  
  # plot all algorithms using merge2
  pdf(file = "Final Scatterplot - merge2.pdf", height = 5.5, width = 3.5)
  for (j in algs) {
    print(plot.formal.scatter(sig = "Vel", alg = j, merge = "merge2", nrow = 4) + labs(y = paste0(j, " Estimate")))     
  }
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
  
    # plot.data$ct <- plot.data$celltype <- substr(rownames(plot.data), start = 1, stop = 3)
    # plot.data$celltype[which(plot.data$celltype == "Exc")] <- paste0(LETTERS, plot.data$celltype[which(plot.data$celltype == "Exc")])
    # plot.data$celltype[which(plot.data$celltype == "Inh")] <- paste0(LETTERS, plot.data$celltype[which(plot.data$celltype == "Inh")])
    # # 
  
  ## Add abundance
    plot.data$Abundance <- colMeans(true$orig.celltype)

  ## Add each celltype maximum non-self correlation
    cor <- cor(sigsSNME$VL$orig.celltype, method = "s")
    plot.data$Sim <- apply(cor, 1, function(x) max(x[which(x != max(x))]))
    
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
    # geom_vline(xintercept = 2, linetype = 2, alpha = 0.6) +
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
      
      pdf(file = filename, height = 2.7, width = 2.5)
      print(ggplot(plot.data, aes(x = variable, y = Linseed, fill = Correlation, label = round(Correlation, 2))) +
              geom_tile(colour = "black") +
              geom_text(size = 2.5) +
              theme_bw() +
              scale_fill_carto_c(palette = "ArmyRose", limits = c(-1,1)) +
              scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
              labs(x = "True Cell-type", y = "Linseed Cell-type") +
              NoLegend() +
              theme(panel.border = invis(), axis.ticks = invis(), panel.grid = invis(), axis.text.x = element_text(angle = 90, hjust = 1)))
      dev.off()
    }  
    
  ## Apply
    # linseed.heatmap(complete$linseed.merge1, true$merge1, "Complete Linseed Heatmap, merge1.pdf")
    linseed.heatmap(complete$linseed.merge2, true$merge2, "Complete Linseed Heatmap, merge2.pdf")
    # linseed.heatmap(complete$linseed.merge1.gA, complete.true$Am1, "Complete Linseed Heatmap, gradientA, merge1.pdf")
    linseed.heatmap(complete$linseed.merge2.gA, complete.true$Am2, "Complete Linseed Heatmap, gradientA, merge2.pdf")
    # linseed.heatmap(complete$linseed.merge1.gB, complete.true$Bm1, "Complete Linseed Heatmap, gradientB, merge1.pdf")
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
        theme(panel.border = invis(), axis.line = element_line()) +
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
    complete.scatter(true$merge2$Neurons_Exc, complete$linseed.merge2$`Cell type 4`, ct = "Excitatory Neurons", ylab = "Linseed Celltype D", colour = plot.information$universal.palette["5Neu-Exc"]) 
    complete.scatter(true$merge2$Oligodendrocytes, complete$linseed.merge2$`Cell type 5`, ct = "Oligodendrocytes", ylab = "Linseed Celltype E", colour = plot.information$universal.palette["4Oli"])
    complete.scatter(true$merge2$OPC, complete$linseed.merge2$`Cell type 6`, ct = "OPCs", ylab = "Linseed Celltype F", colour = plot.information$universal.palette["4OPC"]) 
    dev.off()
    
    pdf(file = "Complete Linseed Scatter, gB.pdf", height = 2, width = 2)
    complete.scatter(complete.true$Bm2$Astrocytes, complete$linseed.merge2.gB$`Cell type 1`, ct = "Astrocytes", ylab = "Linseed Celltype A" , colour = plot.information$universal.palette["1Ast"])
    complete.scatter(complete.true$Bm2$OPC, complete$linseed.merge2.gB$`Cell type 2`, ct = "OPC", ylab = "Linseed Celltype B" , colour = plot.information$universal.palette["4OPC"])
    complete.scatter(complete.true$Bm2$Neurons_Exc, complete$linseed.merge2.gB$`Cell type 3`, ct = "Excitatory Neurons", ylab = "Linseed Celltype C" , colour = plot.information$universal.palette["5Neu-Exc"])
    complete.scatter(complete.true$Bm2$Endothelial, complete$linseed.merge2.gB$`Cell type 4`, ct = "Endothelia", ylab = "Linseed Celltype D" , colour = plot.information$universal.palette["2End"])
    complete.scatter(complete.true$Bm2$Oligodendrocytes, complete$linseed.merge2.gB$`Cell type 5`, ct = "Oligodendrocytes", ylab = "Linseed Celltype E" , colour = plot.information$universal.palette["4Oli"])
    complete.scatter(complete.true$Bm2$Neurons_Inh, complete$linseed.merge2.gB$`Cell type 6`, ct = "Inhibitory Neurons", ylab = "Linseed Celltype F" , colour = plot.information$universal.palette["6Neu-Inh"])
    complete.scatter(complete.true$Bm2$Microglia, complete$linseed.merge2.gB$`Cell type 7`, ct = "Microglia", ylab = "Linseed Celltype G" , colour = plot.information$universal.palette["3Mic"])
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
      
      pdf(file = filename, height = 2.5, width = 2.5)
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
              theme(panel.border = invis(), axis.ticks = invis(), panel.grid = invis(), axis.title.y = invis(), axis.text.x = element_text(angle = 90, hjust = 1)))
      dev.off()
    }
    
    coex.heatmap(e = complete$coex.merge2, u = true$merge2, filename = "Complete Coex Heatmap.pdf")
    coex.heatmap(e = complete$coex.merge2.gB, u = complete.true$Bm2, filename = "Complete Coex Heatmap, gradientB.pdf")
    
    
  
  ## Scatterplot
      pdf(file = "Complete Coex Scatterplot.pdf", height = 2, width = 2)
      complete.scatter(true$merge2$Astrocytes, complete$coex.merge2$Astrocytes, ct = "Astrocytes", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["1Ast"])
      complete.scatter(true$merge2$Neurons_Inh, complete$coex.merge2$Neurons_Inh, ct = "Inhibitory Neurons", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["6Neu-Inh"])
      complete.scatter(true$merge2$Oligodendrocytes, complete$coex.merge2$Oligodendrocytes, ct = "Oligodendrocytes", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["4Oli"])
      dev.off()
      
      
      pdf(file = "Complete Coex Scatterplot, gradientB.pdf", height = 2, width = 2)
      complete.scatter(complete.true$Bm2$Astrocytes, complete$coex.merge2.gB$Astrocytes, ct = "Astrocytes", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["1Ast"])
      complete.scatter(complete.true$Bm2$Neurons_Exc, complete$coex.merge2.gB$Neurons_Exc, ct = "Excitatory Neurons", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["5Neu-Exc"])
      complete.scatter(complete.true$Bm2$Oligodendrocytes, complete$coex.merge2.gB$Oligodendrocytes, ct = "Oligodendrocytes", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["4Oli"])
      complete.scatter(complete.true$Bm2$OPC, complete$coex.merge2.gB$OPC, ct = "OPCs", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["4OPC"])
      complete.scatter(complete.true$Bm2$Endothelial, complete$coex.merge2.gB$Endothelial, ct = "Endothelia", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["2End"])
      complete.scatter(complete.true$Bm2$Microglia, complete$coex.merge2.gB$Microglia, ct = "Microglia", ylab = "Coex Enrichment", axis.min = NA, abline = FALSE, colour = plot.information$universal.palette["3Mic"])
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
              theme(panel.border = invis(), axis.ticks = invis(), panel.grid = invis(), axis.title = invis(), axis.text.x = element_text(angle = 90, hjust = 1))
      
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
    
   qc.corMx(true$merge2, "Complete QC Random Cors.pdf")
   qc.corMx(complete.true$Bm2, "Complete QC Non-Random Cors.pdf")
   
        
  
################################################################################################################################ #
## Formal plots of deconvolution using dropped signatures ----
   
## First, rerun!
   drop <- drop.sub <- list()
   
   for (j in colnames(sigsSNME$VL$merge1)) {
     print(Sys.time())
     drop[[j]] <- run.CIB(from.file = FALSE,
                          sigObject = sigsSNME$VL$merge1[,-which(colnames(sigsSNME$VL$merge1) == j)],
                          mixString = "snme_Vel.txt")
   }
   
      drop.sub$Neurons_Exc <- run.CIB(from.file = FALSE,
                                   sigObject = sigsSNME$VL$merge2[,-which(colnames(sigsSNME$VL$merge2) == "Neurons_Exc")],
                                   mixString = "snme_Vel.txt")
   
   drop.sub$Neurons_Inh <- run.CIB(from.file = FALSE,
                                   sigObject = sigsSNME$VL$merge2[,-which(colnames(sigsSNME$VL$merge2) == "Neurons_Inh")],
                                   mixString = "snme_Vel.txt")
   
   save(drop, drop.sub, file = "Dropped Estimates.rda")
   
## Plots for drops in the main deconvolution list (saved in "est")  
# plot.formal.drop <- function(metric, drop.id = "drop1", nondrop.id = "merge2", match.celltypes = TRUE, nrow = NA, y.lab.seed) {
#   ## Collect stats for the dropped condition
#     drop <- lapply(est$VL[[drop.id]], function(x) {
#       x <- write.stats(e = x,
#                        t = true[[drop.id]],
#                        error = TRUE,
#                        alg = "")
# 
#       x <- x[grep("nmae_|r_", rownames(x)),]
#       x <- t(x)
#       x <- as.data.frame(x)
#       x$Condition <- drop.id
#       return(x)
#     })
# 
#     for(j in names(drop)) drop[[j]]$Alg <- j
#     drop <- do.call("rbind", drop)
#     drop$Celltype <- sapply(strsplit(rownames(drop), "\\."), "[", 2)
# 
#   ## Collect stats for the equivalent deconvolution using a non-dropped condition
#     non.drop <- lapply(est$VL[[nondrop.id]], function(x) {
#       x <- write.stats(e = x,
#                        t = true[[nondrop.id]],
#                        error = TRUE,
#                        alg = "")
# 
#       x <- x[grep("nmae_|r_", rownames(x)),]
#       x <- t(x)
#       x <- as.data.frame(x)
#       x$Condition <- nondrop.id
#       return(x)
#     })
# 
#     for(j in names(non.drop)) non.drop[[j]]$Alg <- j
#     non.drop <- do.call("rbind", non.drop)
#     non.drop$Celltype <- sapply(strsplit(rownames(non.drop), "\\."), "[", 2)
# 
#   ## Filter to CIB
#     drop <- drop[which(drop$Alg == "CIB"),]
#     non.drop <- non.drop[which(non.drop$Alg == "CIB"),]
# 
#   ## Remove any cell-types exclusive to merge
#     if (match.celltypes) non.drop <- non.drop[-which(is.na(match(non.drop$Celltype, drop$Celltype))),]
# 
#   ## Combine these data
#     dat <- data.frame(drop = drop$r_, non.drop = non.drop$r_, Celltype = drop$Celltype)
#     dat$Celltype <- gsub("Neurons_", "", dat$Celltype)
#     dat$Celltype <- substr(dat$Celltype, 1, 3)
# 
#     axis.size <- min(c(dat$drop, dat$non.drop))
#     lab.x <- paste0("Pearson Correlation\n(All Cell-types in Signature)\nMean r = ", round(mean(non.drop$r_), 2))
#     lab.y <- paste0("Pearson Correlation\n(", y.lab.seed, ")\nMean r = ", round(mean(drop$r_), 2))
# 
# 
#     ggplot(dat, aes(x = non.drop, y = drop, label = Celltype)) +
#       geom_text() +
#       theme_bw() +
#       theme(panel.border = invis(), axis.line = element_line()) +
#       geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red") +
#       scale_y_continuous(limits = c(axis.size, 1)) +
#       scale_x_continuous(limits = c(axis.size, 1)) +
#       labs(x = lab.x, y = lab.y)
# }
#   
#   pdf(file = "Dropped Scatterplot.pdf", height = 4, width = 4)
#   plot.formal.drop(metric = "Pearson", "drop1", "merge2", nrow = 1, y.lab.seed = "Inhibitory Neurons Absent From Signature") 
#   plot.formal.drop(metric = "Pearson", "drop2", "merge2", match.celltypes = FALSE, nrow = 1, y.lab.seed = "L4 Neurons Absent From Signature") 
#   plot.formal.drop(metric = "Pearson", "drop3", "orig.celltype", nrow = 1, y.lab.seed = "L4 Neurons Absent From Signature") 
#   plot.formal.drop(metric = "Pearson", "drop4", "merge1", nrow = 1, match.celltypes = TRUE, y.lab.seed = "Oligodendrocytes Absent From Signature") 
#   dev.off()
#   
 
## Plot
  cols <- ggsci::pal_locuszoom()(7)
  cols <- c(cols, "darkorange1")
  names(cols) <- c(gsub("Neurons_", "", colnames(true$merge2))[c(4,5,1,2,3,6,7)], "Neu")
  names(cols) <- substr(names(cols), 1, 3)
  
  
  plot.formal.drop <- function(e = drop, drop.id = "Oligodendrocytes", nondrop.id = "merge1") {
    ## Collect stats for the dropped condition
    d <- write.stats(e = e[[drop.id]], t = true[[nondrop.id]], error = TRUE, alg = "CIB")
    d <- d[grep("nmae_|r_", rownames(d)),]
    d <- t(d)
    d <- as.data.frame(d)
    d$Condition <- "Dropped"
    
    ## Collect stats for the equivalent deconvolution using a non-dropped condition
    non.drop <- write.stats(e = est$VL[[nondrop.id]]$CIB[,rownames(d)], t = true[[nondrop.id]], error = TRUE, alg = "CIB")
    non.drop <- non.drop[grep("nmae_|r_", rownames(non.drop)),]
    non.drop <- t(non.drop)
    non.drop <- as.data.frame(non.drop)
    non.drop$Condition <- "Full"
    
    
    ## Combine these data
    dat <- data.frame(dr = d$r_, nr = non.drop$r_, de = d$nmae_CIB, ne = non.drop$nmae_CIB, Celltype = rownames(d))
    dat$Celltype <- gsub("Neurons_", "", dat$Celltype)
    dat$Celltype <- substr(dat$Celltype, 1, 3)
    
    # axis.sizeA <- min(c(dat$dr, dat$nr))
    axis.sizeA <- 0
    lab.xA <- paste0("r (All Ct in Sig)\nMean = ", round(mean(non.drop$r_), 2))
    lab.yA <- paste0("r (", substr(drop.id, 1, 3), " Absent from Sig)\nMean = ", round(mean(d$r_), 2))
    
    pA <- ggplot(dat, aes(x = nr, y = dr, colour = Celltype)) +
      geom_point(size = 2) +
      scale_colour_manual(values = cols) +
      theme_bw() +
      theme(panel.border = invis(), axis.line = element_line(), axis.title = element_text(size = 8), 
            axis.text = element_text(size = 8), legend.position = "none") +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red") +
      scale_y_continuous(limits = c(axis.sizeA, 1)) +
      scale_x_continuous(limits = c(axis.sizeA, 1)) +
      labs(x = lab.xA, y = lab.yA)
    
    axis.sizeB <- max(c(dat$de, dat$ne)) + 0.1
    lab.xB <- paste0("NMAE (All Ct in Sig)\nMean = ", round(mean(non.drop$nmae_CIB), 2))
    lab.yB <- paste0("NMAE (", substr(drop.id, 1, 3), " Absent from Sig)\nMean = ", round(mean(d$nmae_CIB), 2))
    
    pB <- ggplot(dat, aes(x = ne, y = de, colour = Celltype)) +
      geom_point(size = 2) +
      scale_colour_manual(values = cols) +
      theme_bw() +
      theme(panel.border = invis(), axis.line = element_line(), axis.title = element_text(size = 8), 
            axis.text = element_text(size = 8), legend.position = "none") +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red") +
      scale_y_continuous(limits = c(0, axis.sizeB)) +
      scale_x_continuous(limits = c(0, axis.sizeB)) +
      labs(x = lab.xB, y = lab.yB)
    
    plot_grid(pA, pB, nrow = 1)
  }
  
  pdf(file = "Dropped Scatterplots V2.pdf", height = 2, width = 4)
  for (j in names(drop)) print(plot.formal.drop(drop.id = j))
  plot.formal.drop(e = drop.sub, drop.id = "Neurons_Exc", nondrop.id = "merge2")
  plot.formal.drop(e = drop.sub, drop.id = "Neurons_Inh", nondrop.id = "merge2")
  show_col(cols, labels = FALSE, ncol = 1)
  dev.off()
  
  
 
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
    # pA <- plot.formal.mismatch.scatter(sig = "VL", alg = j, merge = "merge2", no.ylab = FALSE)
    # pB <- plot.formal.mismatch.scatter(sig = "NG", alg = j, merge = "merge2")
    # pC <- plot.formal.mismatch.scatter(sig = "CA", alg = j, merge = "merge2")
    # pD <- plot.formal.mismatch.scatter(sig = "LK", alg = j, merge = "merge2", no.ylab = FALSE, no.xlab = FALSE)
    # pE <- plot.formal.mismatch.scatter(sig = "TS", alg = j, merge = "merge2", no.xlab = FALSE)
    # pF <- get_legend(pA + theme(legend.position = "right"))
    # # pF <- plot.formal.mismatch.scatter(sig = "Dar", alg = "CIB", merge = "merge1", no.xlab = FALSE) # because it's merge1
    # print(plot_grid(pA, pB, pC, pD, pE, pF, nrow = 2, rel_widths = c(1, 0.8, 0.8), rel_heights = c(1, 1.2)))
    
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
    pI <- plot.formal.mismatch.scatter(sig = "DM", alg = "MUS", merge = "merge1", no.xlab = FALSE)
    pJ <- get_legend(pA + theme(legend.position = "right") + guides(colour = guide_legend(ncol = 2)))  
    print(plot_grid(pA, pB, pC, pD, pE, pI, pJ, nrow = 2, ncol = 3, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8), rel_heights = c(1, 1.2)))  
  dev.off()
  
  ## CIB merge1 for the main figure
    pdf(file = paste0("Mismatch Signature - CIB merge1 (MAIN).pdf"), height = 5, width = 6)
    pA <- plot.formal.mismatch.scatter(sig = "VL", alg = "CIB", merge = "merge1", no.ylab = FALSE) + labs(y = "")
    pB <- plot.formal.mismatch.scatter(sig = "NG", alg = "CIB", merge = "merge1")
    pC <- plot.formal.mismatch.scatter(sig = "CA", alg = "CIB", merge = "merge1")
    pD <- plot.formal.mismatch.scatter(sig = "LK", alg = "CIB", merge = "merge1", no.ylab = FALSE)
    pE <- plot.formal.mismatch.scatter(sig = "TS", alg = "CIB", merge = "merge1")
    pF <- plot.formal.mismatch.scatter(sig = "F5", alg = "CIB", merge = "merge1") 
    pG <- plot.formal.mismatch.scatter(sig = "IP", alg = "CIB", merge = "merge1", no.xlab = FALSE, no.ylab = FALSE) + labs(y = "")
    pH <- plot.formal.mismatch.scatter(sig = "MM", alg = "CIB", merge = "merge1", no.xlab = FALSE)
    pI <- plot.formal.mismatch.scatter(sig = "DM", alg = "CIB", merge = "merge1", no.xlab = FALSE)
    # pJ <- get_legend(pA + theme(legend.position = "right") + guides(colour = guide_legend(ncol = 2)))  
    print(plot_grid(pA, pB, pC, pD, pE, pF, pG, pH, pI, nrow = 3, ncol = 3, rel_widths = c(1, 0.8, 0.8), rel_heights = c(1, 1, 1.2)))  
    dev.off()
    
    pdf(file = paste0("Mismatch Signature - CIB merge1 (MAIN LEGEND).pdf"), height = 5, width = 2)
    leg <- get_legend(pA + theme(legend.position = "right"))
    print(plot_grid(leg))
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
    

## Barplots of error
  plot.formal.mismatch.nmae <- function(sig, merge, alg, ylim = 2) {
    dat <- melt(errors[[sig]][[merge]])
    dat <- dat[grep("nmae", dat$Var2),]
    colnames(dat) <- c("Celltype", "Error", "value", "Algorithm")
    dat$Algorithm <- factor(dat$Algorithm)
    levels(dat$Algorithm)[4] <- "MUS"

    dat <- dat[which(dat$Algorithm == alg),]


    # recolour...
    m <- match(dat$Celltype, plot.information[[merge]])
    dat$Class <- names(plot.information[[merge]][m])

    # ...and reorder the x-axis
    dat$Celltype <- factor(dat$Celltype, levels = plot.information[[merge]])
    levels(dat$Celltype) <- gsub(pattern = "Neurons_", replacement = "", levels(dat$Celltype))
    levels(dat$Celltype) <- substr(levels(dat$Celltype), 1, 3)

    # add a plot title
    # label <- paste0(sig, " (mean=", round(mean(dat$value),2), ")")
    label <- sig

    # limits
    if (max(dat$value) > ylim) {
      ylim <- max(dat$value) + 1
    }

    
    dat$Dummy <- "x"
    
    # plot
    ggplot(dat, aes(x = Dummy, y = value, fill = Class, colour = Class)) +
      geom_hline(yintercept = 1, linetype = 2) +
      # geom_col(position = "dodge", colour = "black", width = 0.01) +
      # # geom_col(position = position_dodge(width = 0.8), width = 0.1) +
      # geom_point(position = "dodge", size = 2, shape = 21, colour = "black") +
       geom_col(position = position_dodge(width = 1), colour = "black") +
        # geom_col(position = position_dodge(width = 0.8), width = 0.1) +
        # geom_point(position = position_dodge(width = 0.7), size = 3, shape = 21, colour = "black") +
        stat_summary(geom = "point", fun = function(x) x, position = position_dodge(width = 1), size = 3, shape = 23, colour = "black") +
      
      theme_bw() +
      scale_y_continuous(limits = c(0, ylim), expand = c(0,0)) +
      labs(y = "Normalised Mean\nAbsolute Error", title = label) +
      scale_fill_manual(values = plot.information$universal.palette) +
      scale_colour_manual(values = plot.information$universal.palette) +
      
      theme(axis.text.x = invis, panel.border = invis, axis.line.y = element_line(),
              axis.title.x = invis, legend.position = "none", plot.title = element_text(hjust = 0.5), title = element_text(size = 8))
  }

    pdf(file = paste0("Mismatch Signature Error - CIB merge1 (MAIN).pdf"), height = 1.5, width = 8)
    pA <- plot.formal.mismatch.nmae("VL", "merge1", "CIB", ylim = 2) 
    pB <- plot.formal.mismatch.nmae("NG", "merge1", "CIB", ylim = 2) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pC <- plot.formal.mismatch.nmae("CA", "merge1", "CIB", ylim = 2)+ theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pD <- plot.formal.mismatch.nmae("LK", "merge1", "CIB", ylim = 2) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pE <- plot.formal.mismatch.nmae("TS", "merge1", "CIB", ylim = 2) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pF <- plot.formal.mismatch.nmae("F5", "merge1", "CIB", ylim = 2) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pG <- plot.formal.mismatch.nmae("IP", "merge1", "CIB", ylim = 2) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pH <- plot.formal.mismatch.nmae("DM", "merge1", "CIB", ylim = 2) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pI <- plot.formal.mismatch.nmae("MM", "merge1", "CIB", ylim = 2) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    # pJ <- get_legend(pA + theme(legend.position = "right", legend.title = invis))
    # print(plot_grid(pA, pB, pC, pD, pE, pF, pG, pH, pI, nrow = 3, ncol = 3, rel_widths = c(1.2, 0.8, 0.8)))
    print(plot_grid(pA, pB, pC, pD, pE, pF, pG, pH, pI, nrow = 1, ncol = 9, rel_widths = c(1.2, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)))
    dev.off()
  
  # for (j in algs) {
  #   # setup
  #   if (j == "MUS") next
  #   if (j == "CIB") {
  #     ylim <- 2
  #   } else {
  #     ylim <- 5
  #   }
  # 
  #   # run
  #   pdf(file = paste0("Mismatch Signature Error - ", j, " merge2.pdf"), height = 2, width = 8)
  #   pA <- plot.formal.mismatch.nmae("VL", "merge2", j, ylim = ylim)
  #   pB <- plot.formal.mismatch.nmae("NG", "merge2", j, ylim = ylim) + theme(axis.title.y = invis)
  #   pC <- plot.formal.mismatch.nmae("CA", "merge2", j, ylim = ylim) + theme(axis.title.y = invis)
  #   pD <- plot.formal.mismatch.nmae("LK", "merge2", j, ylim = ylim) + theme(axis.title.y = invis)
  #   pE <- plot.formal.mismatch.nmae("TS", "merge2", j, ylim = ylim) + theme(axis.title.y = invis)
  #   print(plot_grid(pA, pB, pC, pD, pE, nrow = 1, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)))
  #   dev.off()
  # 
  #   pdf(file = paste0("Mismatch Signature Error - ", j, " merge1.pdf"), height = 3, width = 8)
  #   pA <- plot.formal.mismatch.nmae("VL", "merge1", j, ylim = ylim)
  #   pB <- plot.formal.mismatch.nmae("NG", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
  #   pC <- plot.formal.mismatch.nmae("CA", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
  #   pD <- plot.formal.mismatch.nmae("LK", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
  #   pE <- plot.formal.mismatch.nmae("TS", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
  #   pF <- plot.formal.mismatch.nmae("F5", "merge1", j, ylim = ylim)
  #   pG <- plot.formal.mismatch.nmae("IP", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
  #   pH <- plot.formal.mismatch.nmae("DM", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
  #   pI <- plot.formal.mismatch.nmae("MM", "merge1", j, ylim = ylim) + theme(axis.title.y = invis)
  #   pJ <- get_legend(pA + theme(legend.position = "right", legend.title = invis))
  #   print(plot_grid(pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ, nrow = 2, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8)))
  #   dev.off()
  # }




## Barplots of correlation
  ## Function
    plot.formal.mismatch.cor <- function(e = NA, t = NA, sig, alg, merge, autopath = TRUE, no.xlab = TRUE, no.ylab = TRUE) {
      if (autopath) {
        e <- est[[sig]][[merge]][[alg]]
        t <- true[[merge]]
      }
      
      # colnames(e) <- gsub("OPCs", "OPC", colnames(e))
      t <- t[,colnames(e)]
      
      cor <- diag(cor(e, t, method = "p"))
      
      if(anyNA(cor)) cor[which(is.na(cor))] <- 0
      
      dat <- data.frame(Celltype = names(cor), value = cor)
      
      # recolour...
      m <- match(dat$Celltype, plot.information[[merge]])
      dat$Class <- names(plot.information[[merge]][m])
      
      # ...and reorder the x-axis
      dat$Celltype <- factor(dat$Celltype, levels = plot.information[[merge]])
      levels(dat$Celltype) <- gsub(pattern = "Neurons_", replacement = "", levels(dat$Celltype))
      levels(dat$Celltype) <- substr(levels(dat$Celltype), 1, 3)
      
      # add a plot title
      # label <- paste0(sig, " (mean=", round(mean(dat$value),2), ")")
      label <- sig
      
      # plot
      # ggplot(dat, aes(x = Celltype, y = value, fill = Class)) +
      #   geom_col(position = "dodge", colour = "black", width = 0.01) +
      #   # geom_col(position = position_dodge(width = 0.8), width = 0.1) +
      #   geom_point(position = "dodge", size = 2, shape = 21, colour = "black") +
      #   theme_bw() +
      #   scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
      #   labs(y = "Pearson Correlation", title = label) +
      #   scale_fill_manual(values = plot.information$universal.palette) +
      #   geom_hline(yintercept = 0.8, linetype = 2) +
      #   theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), panel.border = invis, axis.line = element_line(),
      #         axis.title.x = invis, legend.position = "none", plot.title = element_text(hjust = 0.5), title = element_text(size = 8))
      dat$Dummy <- "x"
      
      ggplot(dat, aes(x = Dummy, y = value, fill = Class, colour = Class)) +
        geom_col(position = position_dodge(width = 1), colour = "black") +
        # geom_col(position = position_dodge(width = 0.8), width = 0.1) +
        # geom_point(position = position_dodge(width = 0.7), size = 3, shape = 21, colour = "black") +
        stat_summary(geom = "point", fun = function(x) x, position = position_dodge(width = 1), size = 3, shape = 23, colour = "black") +
        theme_bw() +
        scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
        labs(y = "Pearson Correlation", title = label) +
        scale_fill_manual(values = plot.information$universal.palette) +
        scale_colour_manual(values = plot.information$universal.palette) +
        geom_hline(yintercept = 0.8, linetype = 2) +
        theme(axis.text.x = invis, panel.border = invis, axis.line.y = element_line(),
              axis.title.x = invis, legend.position = "none", plot.title = element_text(hjust = 0.5), title = element_text(size = 8))
    }
    
    ## Plot for the main figure
    pdf(file = paste0("Mismatch Signature Correlation - CIB merge1 (MAIN).pdf"), height = 1.5, width = 8)
    pA <- plot.formal.mismatch.cor(sig = "VL", merge = "merge1", alg = "CIB", autopath = TRUE, no.ylab = FALSE)
    pB <- plot.formal.mismatch.cor(sig = "NG", merge = "merge1", alg = "CIB", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pC <- plot.formal.mismatch.cor(sig = "CA", merge = "merge1", alg = "CIB", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pD <- plot.formal.mismatch.cor(sig = "LK", merge = "merge1", alg = "CIB", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pE <- plot.formal.mismatch.cor(sig = "TS", merge = "merge1", alg = "CIB", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pF <- plot.formal.mismatch.cor(sig = "F5", merge = "merge1", alg = "CIB", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pG <- plot.formal.mismatch.cor(sig = "IP", merge = "merge1", alg = "CIB", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pH <- plot.formal.mismatch.cor(sig = "DM", merge = "merge1", alg = "CIB", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    pI <- plot.formal.mismatch.cor(sig = "MM", merge = "merge1", alg = "CIB", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis, axis.line.y = invis)
    # pJ <- get_legend(pA + theme(legend.position = "right", legend.title = invis))
    # print(plot_grid(pA, pB, pC, pD, pE, pF, pG, pH, pI, nrow = 3, ncol = 3, rel_widths = c(1.2, 0.8, 0.8)))
    print(plot_grid(pA, pB, pC, pD, pE, pF, pG, pH, pI, nrow = 1, ncol = 9, rel_widths = c(1.2, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)))
    dev.off()
 
  # ## Apply
  #   for (j in algs) {
  #     if (j == "MUS") next
  #     pdf(file = paste0("Mismatch Signature Correlation - ", j, " merge2.pdf"), height = 2, width = 8)
  #     pA <- plot.formal.mismatch.cor(sig = "VL", alg = "CIB", merge = "merge2", autopath = TRUE, no.ylab = FALSE)
  #     pB <- plot.formal.mismatch.cor(sig = "NG", alg = "CIB", merge = "merge2", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
  #     pC <- plot.formal.mismatch.cor(sig = "CA", alg = "CIB", merge = "merge2", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
  #     pD <- plot.formal.mismatch.cor(sig = "LK", alg = "CIB", merge = "merge2", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
  #     pE <- plot.formal.mismatch.cor(sig = "TS", alg = "CIB", merge = "merge2", autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
  #     print(plot_grid(pA, pB, pC, pD, pE, nrow = 1, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)))
  #     dev.off()
  # 
  #     pdf(file = paste0("Mismatch Signature Correlation - ", j, " merge1.pdf"), height = 3, width = 8)
  #     pA <- plot.formal.mismatch.cor(sig = "VL", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE)
  #     pB <- plot.formal.mismatch.cor(sig = "NG", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
  #     pC <- plot.formal.mismatch.cor(sig = "CA", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
  #     pD <- plot.formal.mismatch.cor(sig = "LK", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
  #     pE <- plot.formal.mismatch.cor(sig = "TS", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
  #     pF <- plot.formal.mismatch.cor(sig = "F5", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE)
  #     pG <- plot.formal.mismatch.cor(sig = "IP", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
  #     pH <- plot.formal.mismatch.cor(sig = "DM", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
  #     pI <- plot.formal.mismatch.cor(sig = "MM", merge = "merge1", alg = j, autopath = TRUE, no.ylab = FALSE) + theme(axis.title.y = invis, axis.text.y = invis, axis.ticks.y = invis)
  #     pJ <- get_legend(pA + theme(legend.position = "right", legend.title = invis))
  #     print(plot_grid(pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ, nrow = 2, rel_widths = c(1, 0.8, 0.8, 0.8, 0.8)))
  #     dev.off()
  #   }
  # 
    
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
  gof$VL <- get.gof(est$VL$merge2$CIB, sigsSNME$VL$merge2)
  gof$CA <- get.gof(est$CA$merge2$CIB, sigsSNME$CA)
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
    # geom_violin(scale = "width", width = 0.9, fill = "white", colour = "black") +
    # geom_boxplot(fill = "grey50", colour = "black") +
    geom_point(position = position_jitter(), alpha = 0.2, size = 0.5) +
    # stat_summary(geom = "point", fun = function(x) mean(x), shape = "-", size = 8, show.legend = FALSE, colour = "firebrick1") +
    theme_bw() +
    theme(panel.grid = invis(), panel.border = invis(), axis.line = element_line()) +
    geom_hline(yintercept = c(0.5, 0.7), linetype = 2) +
    labs(y = "Goodness of Fit (r)", x = "Signature")
  dev.off()
