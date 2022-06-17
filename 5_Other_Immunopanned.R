## Using Zhang's purified cells as mixtures.

############################################################################################################################### #
## Setup ----

## Generic
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)


## Set directory
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/")

## Functions and libraries
source("/Volumes/Data1/PROJECTS/BrainCellularComposition/Scripts/Fun_Preprocessing.R")
source("/Volumes/Data1/PROJECTS/BrainCellularComposition/Scripts/Fun_Composition.R")
library(Seurat)
library(rcartocolor)
library(ggplot2)
invis <- element_blank()

## Files for preprocessing 
load("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/geneInfo.rda") 
load("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/exonicLength.rda") 

## New function
  # to correct for length!
  length.correct <- function(cpm) {
    cpm <- cpm[which(rownames(cpm) %in% rownames(exonicLength)),]  
    
    length <- transpose(exonicLength)[[1]] # note 2: the file "exonicLength" must be loaded into the global environment,
    names(length) <- rownames(exonicLength)
    
    m <- match(rownames(cpm), names(length))
    length <- length[m]/1000
    
    rpkm <- apply(cpm, 2, function(x) x / length)  
    rpkm <- as.data.frame(rpkm)
    return(rpkm)
  }
    

################################################################################################################################ #
## Preprocess Zhang mixtures ----

## These are immunopanned pure brain cell-types, serving as a mixture dataset to deconvolve. True proportion is 1 of the annotated cell-type
  
# downloaded as fpkm data in Table S4, Zhang et al. 2016, and modified to keep only healthy, human, pure cell-type transcriptomes
zhang <- read.csv("../../../Data/Raw/Zhang2016_RNASeq.csv")

dup <- which(duplicated(zhang$GeneSymbol))
zhang <- zhang[-dup,]

# reannotate gene names
rownames(zhang) <- zhang$GeneSymbol
zhang <- zhang[,-which(colnames(zhang) == "GeneSymbol")]
zhang <- addENSID(zhang)

# filter out cell-types
zhang <- zhang[,-grep("Foetal", colnames(zhang))]

# apply minimum threshold of 1 fpkm in any one cell-type
exp_thresh <- 1
zhang <- zhang[which(apply(zhang, 1, max) > exp_thresh),] # 12k genes is reasonable. only necessary in one sample as there's only one neuronal sample...

# write to disk
write.CIB(zhang, "../../../Data/Preprocessed/CIB/Zhang_PureSamples.txt")


################################################################################################################################ #
## Deconvolve ----

## Load signatures
  # load
  load("../../../Data/Preprocessed/Signatures - Brain.rda") 

  # filter to key signatures: IP, SC, and HCA in sigsBrain
  rm(sigsRME)
    
  sigsBrain <- sigsBrain[c(2,3,5,6,7,9)]
  # sigsBrain$LK <- length.correct(sigsBrain$LK)
  
## Load second batch of signatures
  load("../../../Data/Preprocessed/snme_Signatures.rda")
  ct <- c("Astrocytes", "Microglia", "Neurons", "Oligodendrocytes")
  sigsBrain$VL <- length.correct(addENSID(sigs.snme$Vel$merge1[,ct]))
  sigsBrain$NG <- length.correct(addENSID(sigs.snme$Nagy$merge1[,ct]))
  rm(sigs.snme); gc()
  
  
# ## Load third batch of signatures
#   load("SN/Lister Signature, Integrated, Nagy Annot.rda")
#   
#   # add ensemble gene ID
#   listerXnagy <- addENSID(listerXnagy)
#   listerXnagy <- listerXnagy[,-3]
#   colnames(listerXnagy) <- c("Oligodendrocytes", "Astrocytes", "Neurons")
#   
#   # adjust cpm for exonic length
#   sigsBrain$LT <- length.correct(listerXnagy) # also removes OPC
#   
# ## Filter to Ast, Mic, Neu, Oli
#   sigsBrain <- lapply(sigsBrain, function(x) {
#     x <- x[,which(colnames(x) %in% c("Oligodendrocytes", "Astrocytes", "Neurons", "Microglia"))]
#   })
  
## Deconvolve
  # setup
  est <- list()
  est$DRS <- est$CIB <- list()
  
  # deconvolve(~15min per interation)
  for (j in names(sigsBrain)) {
    print(Sys.time())
    
    est$DRS[[j]] <- run.DRS(zhang, sigsBrain[[j]])
    # rownames(est$DRS[[j]]) <- colnames(zhang)

    est$CIB[[j]] <- run.CIB(from.file = FALSE,
                                sigObject = sigsBrain[[j]],
                                mixString = "Zhang_PureSamples.txt")

    est$DTA[[j]] <- run.dtangle(mixture = zhang, 
                                    signature = sigsBrain[[j]],
                                    alg = "diff",
                                    q = 0.01)
  }
  
  save(est, file = "Zhang Estimates (Update 2).rda")
  load(file = "Zhang Estimates (Update 2).rda")
  
  # extract estimate for the relevant celltype
  get.purity <- function(search, dat = y) {
    if (length(grep(search, colnames(dat))) == 0) return()
    data.frame(Proportion = dat[grep(search, rownames(dat)), grep(search, colnames(dat))],
               Celltype = colnames(dat)[grep(search, colnames(dat))])
  }
  
  ct <- substr(colnames(sigsBrain$IP), 1, 3)

  est <- lapply(est, function(x) {
    x <- lapply(x, function(y) {
      y <- lapply(ct, function(z) {
        get.purity(search = z, dat = y)
      }) 
      y <- do.call("rbind", y)
    })
  })
  
## Plot 
  purity.plot <- function(e, alg) {
    plot.data <- melt(e)
    plot.data$Celltype <- substr(plot.data$Celltype, 1, 3)
    plot.data$Celltype <- factor(plot.data$Celltype)
    levels(plot.data$Celltype) <- paste0(levels(plot.data$Celltype), "\n(n=", c(12,3,1,5), ")")
    
    plot.data <- plot.data[which(plot.data$L2 != "LT"),]
    plot.data$L2 <- factor(plot.data$L2, levels = c("IP", "SC", "CA", "LK", "VL", "NG"))
    levels(plot.data$L2) <- c("Whole Cell (IP) (Matched)", "Whole Cell (SC)", "Nuclear (CA)", "Nuclear (LK)",
                              "Nuclear (VL)", "Nuclear (NG)")
    
    ggplot(plot.data[which(plot.data$L1 == alg),], aes(x = L2, y = value, fill = L2, colour = L2)) +
      stat_summary(geom = "point", fun = function(x) {x}, position = position_jitter(width = 0.2, height = 0), alpha = 0.4, size = 1.5) +
      stat_summary(geom = "point", fun = function(x) mean(x), shape = "-", size = 8, show.legend = FALSE) +
      facet_wrap(~Celltype, strip.position = "bottom", nrow = 1) +
      theme_bw() +
      labs(y = "Estimated Proportion\nin Pure Whole Cell Samples", title = alg) +
      # theme(panel.border = invis, axis.line.y = invis) +
      # geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 2, colour = "grey50") +
      scale_colour_manual(values = c("firebrick1", "darkorange1", "dodgerblue", "deepskyblue", "darkcyan", "darkolivegreen4", "limegreen")) +
      scale_y_continuous(expand = c(0,0), limits = c(0, 1.01)) +
      theme(axis.text.x = invis, axis.title.x = invis, axis.ticks.x = invis, panel.grid.major.x = invis, legend.title = invis,
            panel.grid.minor.y = invis) 
  }
  
  pdf(file = "Zhang Scatterplot (Update).pdf", height = 3, width = 8)
  purity.plot(est, "CIB") + labs(title = "")
  purity.plot(est, "DRS")
  purity.plot(est, "DTA")
  dev.off()  
  
## Plot just IP
  plot.data <- as.data.frame(sapply(est, function(x) x[["IP"]][,1]))
  plot.data$Celltype <- est$CIB$IP$Celltype
    plot.data$Celltype <- substr(plot.data$Celltype, 1, 3)
    plot.data$Celltype <- factor(plot.data$Celltype)
    levels(plot.data$Celltype) <- paste0(levels(plot.data$Celltype), "\n(n=", c(12,3,1,5), ")")
  
  plot.data <- melt(plot.data)
  
  pdf(file = "Zhang IP Scatterplot.pdf", height = 3, width = 8)
  ggplot(plot.data, aes(x = variable, y = value, colour = variable)) +
    geom_point(position = position_jitter(width = 0.2)) +
    stat_summary(geom = "point", fun = function(x) mean(x), shape = "-", size = 8, show.legend = FALSE, colour = "black") +
    theme_bw() +
    facet_wrap(~Celltype, nrow = 1) +
    labs(y = "Estimated Proportion of Pure Cell-type") +
    theme(legend.position = "none", axis.title.x = invis) +
    scale_y_continuous(limits = c(0, 1.05), expand = c(0,0))
  dev.off()
  
  
  
  
################################################################################################################################ #
## Strategies to improve performance #1: Stable genes! ----
  
## Determine stable genes
  # first, load our nuclear/wholecell data
  load("Normalised RNAseq.rda")

  ## Run a nuclear vs. whole cell comparison to determine compartment-specificity of genes!
    dat <- t(log2(norm$cpm + 0.5))
    prep <- c(rep("nuclei", 5), rep("cell", 5))
    
    mod <- lm(dat ~ prep)
    sum <- summary(mod)
    
    pvals <- sapply(sum, function(y)  { y$coefficients["prepnuclei", "Pr(>|t|)"] } )
    log2fc <- sapply(sum, function(y)  { y$coefficients["prepnuclei", "Estimate"] } )
    # t <- sapply(sum, function(y)  { y$coefficients["group1", "t value"] } )
    res <- data.frame(pvals = pvals, log2fc = log2fc)
    rownames(res) <- sapply(strsplit(rownames(res), " "), "[", 2)
    
    res$padj <- p.adjust(res$pvals, method = "fdr")
    
  ## Get external lists of stable genes
    stable.genes <- list()
    all.genes <- rownames(res) # same as rownames(res)
    
    ## First: load in results from Price et al., 2019.
      # load
      price <- read.csv("Price2019_S2.csv")
      price <- price$Ensembl.ID[grep("Adult|Both", price$Group)] # the EnsID of any gene with compartment specificity in adult brains
      
      # list 1 is the notPrice
      stable.genes$notPrice <- all.genes[-which(all.genes %in% price)]  
      
    ## Second: genes that are not significantly DE at padj <0.05 and fc > 1.3
      x <- which(res$padj < 0.05 & abs(res$log2fc) > log2(1.3))
      stable.genes$notDE <- all.genes[-x]
      
    ## Third: genes that are p > 0.05 & fc < 1.3
      x <- which(res$pvals > 0.05 & abs(res$log2fc) < log2(1.3))
      stable.genes$stable <- all.genes[x]
      
    ## Output
      save(stable.genes, file = "Stable Genes.rda")
      
      res$Stable <- !(res$padj < 0.05 & abs(res$log2fc) > log2(1.3))
      write.csv(res, file = "Supplementary Table - DE Between Nuclear and Whole Cells.csv")
      
      
      
## Deconvolve
  est.stable <- vector("list", length = length(stable.genes))
  names(est.stable) <- names(stable.genes)
  est.stable <- lapply(est.stable, function(x) list(CIB = list(), DRS = list(), DTA = list()))
  
  # deconvolve
  for (i in names(stable.genes)) {
    for (j in names(sigsBrain)) {
      # setup
      print(Sys.time())
      
      # select stable signature, based on i and j
      stable.sig <- sigsBrain[[j]][stable.genes[[i]][which(stable.genes[[i]] %in% rownames(sigsBrain[[j]]))],]
      
      est.stable[[i]]$DRS[[j]] <- run.DRS(zhang, stable.sig)
      # rownames(est.stable[[i]]$DRS[[j]]) <- colnames(zhang)

      est.stable[[i]]$CIB[[j]] <- run.CIB(from.file = FALSE,
                                     sigObject = stable.sig,
                                     mixString = "Zhang_PureSamples.txt")

       est.stable[[i]]$DTA[[j]] <- run.dtangle(mixture = zhang, 
                                  signature = stable.sig,
                                  alg = "diff",
                                  q = 0.01)
    }  
  }
  
  
  save(est.stable, file = "Zhang Estimates (Stable) (Update 2).rda")
  load("Zhang Estimates (Stable) (Update 2).rda")
    
## Plot 
  est.stable <- lapply(est.stable$notDE, function(x) {
      x <- lapply(x, function(y) {
        y <- lapply(ct, function(z) {
          get.purity(search = z, dat = y)
        }) 
        y <- do.call("rbind", y)
      })
    })
  
  pdf(file = "Zhang Scatterplot (Stable) (Update) (notDE).pdf", height = 3, width = 8)
  purity.plot(est.stable, "CIB") + labs(title = "")
  purity.plot(est.stable, "DRS")
  purity.plot(est.stable, "DTA")
  dev.off()  
  
## Check percentage that are above 80% hits
  above.8 <- lapply(est, function(x) {
    sapply(x, function(y) length(which(y$Proportion > 0.8)) / nrow(y))  
  })
  above.8 <- do.call("rbind", above.8)
  
  above.8.st <- lapply(est.stable, function(x) {
    sapply(x, function(y) length(which(y$Proportion > 0.8)) / nrow(y))  
  })
  above.8.st <- do.call("rbind", above.8.st)
  
  
# ## New metric: effect of stable genes
#  
#   
#   a <- melt(est$CIB)
#   b <- melt(est.stable$CIB)
# 
#   
#   plot.data <- data.frame(Celltype = a$Celltype,
#                           Signature = a$L1,
#                           notDE = b$value - a$value) # proportion in stable.notDE deconv - in allgenes deconv, thus the gain in accuracy
#                           
#   
#   plot.data <- melt(plot.data)
#     
#   plot.data$Celltype <- substr(plot.data$Celltype, 1, 3)
#   plot.data$Celltype <- factor(plot.data$Celltype)
#   levels(plot.data$Celltype) <- paste0(levels(plot.data$Celltype), " (n=", c(12,3,1,5), ")")
#   
#   plot.data <- plot.data[which(plot.data$Signature != "LT"),]
#   plot.data$Signature <- factor(plot.data$Signature, levels = c("IP", "SC", "CA", "LK", "VL", "NG"))
#   
#   levels(plot.data$Signature) <- c("Whole Cell Signature (Matched)", "Whole Cell Signature (SC)", "Nuclear Signature (CA)", "Nuclear Signature (LK)",
#                                    "Nuclear Signature (VL)", "Nuclear Signature (NG)")
#   colnames(plot.data)[3] <- "GeneSelection"
#   
#   pdf(file = "Zhang Gain From Stability (Lister DE Only).pdf", height = 3, width = 8)
#   ggplot(plot.data, aes(y = value, x = Celltype, colour = Signature)) +
#     geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2, jitter.height = 0), alpha = 0.3) +
#     # facet_wrap(~Celltype, nrow = 1) +
#     stat_summary(geom = "point", position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2, jitter.height = 0), fun = function(x) mean(x), shape = "-", size = 8, show.legend = FALSE) +
#     theme_bw() +
#     labs(y = "Change In Proportion\nWhen Selecting Stable Genes") +
#     geom_hline(yintercept = 0) +
#         scale_colour_manual(values = c("firebrick1", "darkorange1", "dodgerblue", "deepskyblue", "darkcyan", "darkolivegreen4", "limegreen")) +
#     theme(panel.border = invis, axis.line.y = element_line())
#   dev.off()
#   
#   
# ## New plot: pooled
#   # plot.data <- melt(est)
#   # # plot.data$value <- 1 - plot.data$value # converts proportion to error, assuming that the true proportion is 1
#   # plot.data$Celltype <- substr(plot.data$Celltype, 1, 3)
#   # plot.data$Celltype <- factor(plot.data$Celltype)
#   # levels(plot.data$Celltype) <- paste0(levels(plot.data$Celltype), "\n(n=", c(12,3,1,5), ")")
#   # 
#   # plot.data$L2 <- factor(plot.data$L2, levels = c("IP", "SC", "CA", "LK", "VL", "NG", "LT"))
#   # levels(plot.data$L2) <- c("Whole Cell (IP) (Matched)", "Whole Cell (SC)", "Nuclear (CA)", "Nuclear (LK)",
#   #                           "Nuclear (VL)", "Nuclear (NG)", "Nuclear (LT)")
#   # 
#   # purity.plot <- function(alg) {
#   #   ggplot(plot.data[which(plot.data$L1 == alg),], aes(x = L2, y = value, fill = L2, colour = L2)) +
#   #     stat_summary(geom = "point", fun = function(x) {x}, position = position_jitter(width = 0.2, height = 0), alpha = 0.4, size = 1.5) +
#   #     stat_summary(geom = "point", fun = function(x) mean(x), shape = "-", size = 8, show.legend = FALSE) +
#   #     facet_wrap(~Celltype, strip.position = "bottom", nrow = 1) +
#   #     theme_bw() +
#   #     labs(y = "Estimated Proportion\nin Pure Whole Cell Samples", title = alg) +
#   #     # theme(panel.border = invis, axis.line.y = invis) +
#   #     # geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 2, colour = "grey50") +
#   #     scale_colour_manual(values = c("firebrick1", "darkorange1", "dodgerblue", "deepskyblue", "darkcyan", "darkolivegreen4", "limegreen")) +
#   #     scale_y_continuous(expand = c(0,0), limits = c(0, 1.01)) +
#   #     theme(axis.text.x = invis, axis.title.x = invis, axis.ticks.x = invis, panel.grid.major.x = invis, legend.title = invis,
#   #           panel.grid.minor.y = invis) 
#   # }
#   
#  
#   
# ## Attempt as a scatterplot!
#   # get non-stable deconvolution
#   a <- lapply(est, function(x) {
#     x <- x[-which(names(x) %in% c("SC", "LT"))]
#     sapply(x, function(y) y[,1])
#   })
#   
#   # get stable deconvolution
#   b <- lapply(est.stable, function(x) {
#     x <- x[-which(names(x) %in% c("SC", "LT"))]
#     x <- sapply(x, function(y) y[,1])
#     colnames(x) <- paste0(colnames(x), ".stable")
#     return(x[,-1]) # IP
#   })
#   
#   
#   
#   # function
#   plot.zhang.scatter <- function(sig, lims = 0.2) {
#     cor <- cors[sig]
#     
#     ggplot(c, aes_string(x = "IP", y = sig, colour = "Celltype")) +
#       geom_point() +
#       geom_abline(intercept = 0, slope = 1) +
#       scale_x_continuous(limits = c(lims, 1)) +
#       scale_y_continuous(limits = c(lims, 1)) +
#       NoLegend() +
#       labs(x = paste0("IP\nr=", cors[sig]))
#     
#   }
#   
#   ## Plot CIB
#   for (j in c("CIB", "DRS", "DTA")) {
#     c <- as.data.frame(cbind(a[[j]], b[[j]]))
#     cors <- round(cor(c)[-1,1],2)
#     c$Celltype <- est$CIB$IP$Celltype
#     
#     plot.list <- list(plot.zhang.scatter("LK"), plot.zhang.scatter("LK.stable"),
#                       plot.zhang.scatter("CA"), plot.zhang.scatter("CA.stable"),
#                       plot.zhang.scatter("VL"), plot.zhang.scatter("VL.stable"),
#                       plot.zhang.scatter("NG"), plot.zhang.scatter("NG.stable"))
#     pdf(file = paste0("Zhang Scatterplot, Irina's Request ", j, ".pdf"), height = 11, width = 8)
#     print(plot_grid(plotlist = plot.list, ncol = 2))
#     dev.off()
#     
#   }
#     
#   
#   
#   
#   
#   
#   
#   for(j in names(cors)) {
#         
#   }
#   
#   
#   
#   
#   
#   # # first, get paired information for all the standard (all gene) deconvolution
#   # plot.ag <- lapply(est, function(x) {
#   #   y <- x[[1]][,1]
#   #   x <- x[-which(names(x) %in% c("IP", "SC", "LT"))]
#   #   
#   #   x <- lapply(x, function(z) {
#   #     cbind(z[,1], y)
#   #   })
#   #   x <- do.call("rbind", x)
#   # })
#   # 
#   # # repeat for stable deconvolution
#   # plot.st <- lapply(est.stable, function(x) {
#   #   x <- x[-which(names(x) %in% c("IP", "SC", "LT"))]
#   #   y <- x[[1]][,1]
#   #   x <- lapply(x, function(z) {
#   #     cbind(z[,1], y)
#   #   })
#   #   x <- do.call("rbind", x)
#   # })
#   # 
#   # # combine 
#   # plot.zhang.scatter <- function(alg, y.min = 0.6) {
#   #   plot.data <- as.data.frame(plot.ag[[alg]])
#   #   colnames(plot.data) <- c("IP", "Nuclear")
#   #   plot.data$NuclearStable <- plot.st[[alg]][,1]
#   #   plot.data$Celltype <- est$CIB$IP$Celltype # autorepeats
#   #   
#   #   # add signature information
#   #   plot.data$Signature <- rep(c("LK", "CA", "VL", "NG"), each = nrow(est$CIB$IP))
#   #   plot.data$Signature <- factor(plot.data$Signature)
#   #   plot.data$Signature2 <- plot.data$Signature
#   #   for (j in levels(plot.data$Signature)) {
#   #     use <- plot.data[which(plot.data$Signature == j),]
#   #     cor <- round(cor(use$IP, use$Nuclear), 2)
#   #     levels(plot.data$Signature)[which(levels(plot.data$Signature) == j)] <- paste0(j, " (r=", cor, ")")
#   #     
#   #     cor <- round(cor(use$IP, use$NuclearStable), 2)
#   #     levels(plot.data$Signature2)[which(levels(plot.data$Signature2) == j)] <- paste0(j, " (r=", cor, ")")
#   #   }
#   #   
#   #   pA <- ggplot(plot.data, aes(x = IP, y = Nuclear, colour = Celltype)) +
#   #     geom_point() +
#   #     facet_wrap(~Signature, ncol = 1) +
#   #     scale_y_continuous(limits = c(y.min, 1)) +
#   #     # theme(axis.title.y = invis) +
#   #     labs(title = alg, y = "All Genes' Estimate", x = "IP Estimate")
#   #   
#   #   pB <- ggplot(plot.data, aes(x = IP, y = NuclearStable, colour = Celltype)) +
#   #     geom_point() +
#   #     facet_wrap(~Signature2, ncol = 1) +
#   #     scale_y_continuous(limits = c(y.min, 1)) +
#   #     # theme(axis.title.y = invis) +
#   #     labs(title = "", y = "Stable Deconvolution's Estimate", x = "IP Estimate")
#   #   
#   #   print(plot_grid(pA + NoLegend(), pB, ncol = 2, rel_widths = c(1, 2)))
#   # }
#   # 
#   # plot.zhang.scatter("CIB", y.min = 0.2)
#   # plot.zhang.scatter("DRS", y.min = 0.2)
#   # plot.zhang.scatter("DTA", y.min = 0.2)
#   # 
#   
 
################################################################################################################################ #
## Output supplementary table ----
 
make.supp.table <- function(e) {
  dat <- lapply(e, function(x) {
    x <- x[-which(names(x) == "LT")]
    names(x) <- gsub("SC", "DM", names(x))
    
    label <- x[[1]]$Celltype
    x <- sapply(x, function(y) round(y[,1], 3))
    rownames(x) <- label
    
    # x
    
    return(as.data.frame(x))
  })
  
  dat <- do.call("cbind", dat)
  dat["Mean",] <- round(colMeans(dat), 3)
  dat["Above_0.8",] <- round(apply(dat[-nrow(dat),], 2, function(x) length(which(x > 0.8))   / length(x) ), 3)
  
  return(dat)
}  
   
dat <- make.supp.table(est)
dat.stab <- make.supp.table(est.stable)

write.csv(dat, "Supplementary Table - Zhang Estimates.csv")
write.csv(dat.stab, "Supplementary Table - Zhang Stable Estimates.csv")
  

################################################################################################################################ #
## Apply findings to the Lister data ----
  
  
## Here, I am going to redeconvolve the Lister data. The whole cell data will be deconvolved by IP, the nuclear by individual signatures
  
## Load signatures
  load("SN/Lister Signature, Library-specific, Nagy Annot.rda")  
  ind.sigs <- lapply(ind.sigs, addENSID)
  ind.sigs <- ind.sigs[-grep("55", names(ind.sigs))]   # remove C55
  
## Load expression data
  load("Normalised RNAseq.rda")
  
## Main deconvolution: Zhang signature on c samples, Individual-specific signatures on n samples...
  lister.est <- list()
  
  # setup signatures
  ct <- c("Neu", "Ast", "Oli") # we have to use the union of available ct, which is these three
  
  ip.sig.for.lister <- sigsBrain$IP
  colnames(ip.sig.for.lister) <- substr(colnames(ip.sig.for.lister), 1, 3)
  ip.sig.for.lister <- ip.sig.for.lister[,ct]
  
  ind.sigs <- lapply(ind.sigs, function(x) {
    x[,ct]
  }) 
    
  
  # ...deconvolve whole cells using IP signature
  lister.est$c <- run.CIB(from.file = FALSE,
                           sigObject = ip.sig.for.lister,
                           mixString = "Crushed_all_rpkm.txt")
  
  lister.est$cs <- run.CIB(from.file = FALSE,
                                  sigObject = ip.sig.for.lister[which(rownames(ip.sig.for.lister) %in% stable.genes$notDE),],
                                  mixString = "Crushed_all_rpkm.txt")
    
  
  # ...deconvolve nuclei using individual signatures
    lister.est$ns <- lister.est$n <- list()
  
    for (j in names(ind.sigs)) {
      print(j)
      Sys.time()
      
      k <- substr(j, start = 2, stop = 10)
      
      lister.est$n[[j]] <- run.CIB(from.file = FALSE,
                                   sigObject = ind.sigs[[j]],
                                   mixString = paste0("Crushed_", k, "_cpm.txt"))
      
      lister.est$ns[[j]] <- run.CIB(from.file = FALSE,
                                    sigObject = ind.sigs[[j]][which(rownames(ind.sigs[[j]])  %in% stable.genes$notDE),],
                                    mixString = paste0("Crushed_", k, "_cpm.txt"))
    }
    
    lister.est$n  <- do.call("rbind", lister.est$n)
    rownames(lister.est$n) <- sapply(strsplit(rownames(lister.est$n), "\\."), "[", 2)
    
    lister.est$ns  <- do.call("rbind", lister.est$ns)
    rownames(lister.est$ns) <- sapply(strsplit(rownames(lister.est$ns), "\\."), "[", 2)
    
  # combine the two deconvolutions
    lister.est$combined <- rbind(lister.est$c[grep("_c", rownames(lister.est$c)),], 
                                 lister.est$n[grep("_n", rownames(lister.est$n)),])
    lister.est$combined.stable <- rbind(lister.est$cs[grep("_c", rownames(lister.est$cs)),], 
                                        lister.est$ns[grep("_n", rownames(lister.est$ns)),])
    lister.est$combined.stableNucOnly <- rbind(lister.est$c[grep("_c", rownames(lister.est$c)),], 
                                               lister.est$ns[grep("_n", rownames(lister.est$ns)),])
    lister.est$combined.stableCellOnly <- rbind(lister.est$cs[grep("_c", rownames(lister.est$cs)),], 
                                                lister.est$n[grep("_n", rownames(lister.est$n)),])
 
# ## Deconvolution with other signatures
#   ## Deconvolve
#     lister.est$Nagy <- run.CIB(from.file = FALSE,
#                                sigObject = sigsBrain$NG,
#                                mixString = paste0("Crushed_all_rpkm.txt"))
#     lister.est$Nagy.stable <- run.CIB(from.file = FALSE,
#                                sigObject = sigsBrain$NG[which(rownames(sigsBrain$NG) %in% stable.genes$notDE),],
#                                mixString = paste0("Crushed_all_rpkm.txt"))
#     lister.est$HCA <- run.CIB(from.file = FALSE,
#                                sigObject = sigsBrain$CA,
#                                mixString = paste0("Crushed_all_rpkm.txt"))
#     lister.est$HCA.stable <- run.CIB(from.file = FALSE,
#                                sigObject = sigsBrain$CA[which(rownames(sigsBrain$CA) %in% stable.genes$notDE),],
#                                mixString = paste0("Crushed_all_rpkm.txt"))
#     lister.est$Vel <- run.CIB(from.file = FALSE,
#                                sigObject = sigsBrain$VL,
#                                mixString = paste0("Crushed_all_rpkm.txt"))
#     lister.est$Vel.stable <- run.CIB(from.file = FALSE,
#                                sigObject = sigsBrain$VL[which(rownames(sigsBrain$VL) %in% stable.genes$notDE),],
#                                mixString = paste0("Crushed_all_rpkm.txt"))
#     lister.est$Zhang <- run.CIB(from.file = FALSE, # here, I note that Zhang and "c" both use the same signature, but the former also includes microglia
#                                sigObject = sigsBrain$IP,
#                                mixString = paste0("Crushed_all_rpkm.txt"))
#     lister.est$Zhang.stable <- run.CIB(from.file = FALSE,
#                                sigObject = sigsBrain$IP[which(rownames(sigsBrain$IP) %in% stable.genes$notDE),],
#                                mixString = paste0("Crushed_all_rpkm.txt"))
#     
#     
#     save(lister.est, file = "Zhang, Lister Reanalyses, Estimates.rda")
    
  
  ## Now scatterplot this...
    a <- t(lister.est$c[grep("_c", rownames(lister.est$c)),]) # Cells by Cellular Signature
    b <- t(lister.est$n[grep("_c", rownames(lister.est$n)),]) # Cells by Nuclear Signature
    c <- t(lister.est$ns[grep("_c", rownames(lister.est$ns)),]) # Cells by Stable-Nuclear Signature
    
    a <- melt(a); b = melt(b); c = melt(c)
    a$Deconvolution <- "Cellular Signature"
    b$Deconvolution <- "Nuclear Signature"
    c$Deconvolution <- "Stable Nuclear Signature"
    
    rownames(a) <- paste0(a$Var1, a$Var2)
    rownames(b) <- paste0(b$Var1, b$Var2)
    rownames(c) <- paste0(c$Var1, c$Var2)
    
    a <- a[rownames(b),]
    d <- data.frame(Celltype = a$Var1,
                    Individual = substr(a$Var2, 1, 4),
                    Cellular = a$value,
                    Nuclear = b$value,
                    NuclearStable = c$value)
    
    pdf(file = "Zhang, Lister Reanalyses, Scatterplot.pdf", height = 3.5, width = 4)
    cor <- round(cor(d$Cellular, d$Nuclear), 2)
    ggplot(d, aes(x = Cellular, y = Nuclear, colour = Celltype, shape = Individual)) +
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      geom_point(size = 3) + 
      # theme_bw() +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(limits = c(0, 1)) +
      theme(panel.border = element_blank()) +
      labs(x = paste0("Whole-cell Signature\nr=", cor), y = "Nuclear Signature")
    
    cor <- round(cor(d$Cellular, d$NuclearStable), 2)
    
    ggplot(d, aes(x = Cellular, y = NuclearStable, colour = Celltype, shape = Individual)) +
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      geom_point(size = 3) + 
      # theme_bw() +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(limits = c(0, 1)) +
      theme(panel.border = element_blank()) +
      labs(x = paste0("Whole-cell Signature\nr=", cor), y = "Stable Nuclear Signature")
    dev.off()
  
  