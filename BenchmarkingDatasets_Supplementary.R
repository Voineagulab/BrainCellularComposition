## This script has two separate becnhmarking datasets, in which experimentally-derived mixtures with known proportion are deconvolved

                                                  ## Part 1
                                #### Using Zhang's purified cells as mixtures. ####

############################################################################################################################### #
## Setup ----

## Start!
  rm(list = ls()); gc()

## Directory
  # root.dir <- "INSERT.YOUR.ROOT.DIRECTORY.PLEASE/"
  root.dir <- "/Volumes/Data1/PROJECTS/BrainCellularComposition/BrainCellularComposition_GitHub/"
  setwd(paste0(root.dir, "Results/BenchmarkingDatasets/Supplementary/"))

## Libraries and functions
  source("../../../Scripts/Fun_Composition.R")
  source("../../../Scripts/Fun_Preprocessing.R")
  load("../../../Data/Preprocessed/geneInfo.rda")
  load("../../../Data/Preprocessed/exonicLength.rda")
  library(Seurat)
  library(rcartocolor)
  library(ggplot2)
  invis <- element_blank()
  

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
zhang <- readxl::read_xlsx("../../../Data/Raw/1-s2.0-S0896627315010193-mmc3.xlsx", sheet = 2, col_names = FALSE)
zhang <- as.data.frame(zhang)

 annot <- as.character(zhang[1,])
  annot[1] <- "Symbol"
  for (j in 1:length(annot)) {
    if (is.na(annot[j])) { # if NA
      closest <- which(!(is.na(annot[1:j]))) %>% max() # find the closest non NA to the left
      annot[j] <- annot[closest] # and assign it
    }
  }

  colnames(zhang) <- annot

## Other formatting
  # remove rows 1:3: annotation
  zhang <- zhang[-c(1:3),]
  
  # set rownames to be $Symbol (column 1)
  dup <- which(duplicated(zhang$Symbol))
  zhang <- zhang[-dup,]
  
  rownames(zhang) <- zhang$Symbol
  zhang <- zhang[,-1]

  # convert to EnsID
  zhang <- addENSID(zhang)
  
  # make expression numeric
  labs <- rownames(zhang)
  zhang <- apply(zhang, 2, as.numeric) %>% as.data.frame()
  rownames(zhang) <- labs

# filter samples
  remove <- grep("Endo|whole|GBM|tumor|fetal|hippo", colnames(zhang))
  zhang <- zhang[,-remove]
  
  colnames(zhang) <- gsub("Human ", "", colnames(zhang))
  colnames(zhang) <- gsub("mature a", "A", colnames(zhang))
  colnames(zhang) <- gsub("/Macrophage", "", colnames(zhang))

# apply minimum threshold of 1 fpkm in any one cell-type
exp_thresh <- 1
zhang <- zhang[which(apply(zhang, 1, max) > exp_thresh),] # 12k genes is reasonable. only necessary in one sample as there's only one neuronal sample...

# write to disk
write.CIB(zhang, "../../../Data/Preprocessed/CIB/Zhang_PureSamples.txt")


################################################################################################################################ #
## Deconvolve ----

## Load signatures
  ct <- c("Astrocytes", "Microglia", "Neurons", "Oligodendrocytes")
  
# load
  load("../../../Data/Preprocessed/Signatures - Brain.rda") 

  # filter to key signatures and cell-types
  sigsBrain <- lapply(sigsBrain[c(2,3,5,6,7,9)], function(x) x[,ct])
  

  
 
  
## Deconvolve
  # setup
  est <- list()
  est$DRS <- est$CIB <- est$DTA <- list()
  
  # deconvolve(~15min per interation)
  for (j in names(sigsBrain)) {
    print(Sys.time())
    
    est$DRS[[j]] <- run.DRS(zhang, sigsBrain[[j]])
    # rownames(est$DRS[[j]]) <- colnames(zhang)

    est$CIB[[j]] <- run.CIB(from.file = FALSE,
                                sigObject = sigsBrain[[j]],
                                mixString = "Zhang_PureSamples.txt")

    est$DTA[[j]] <- run.DTA(mixture = zhang, 
                                    signature = sigsBrain[[j]],
                                    alg = "diff",
                                    q = 0.01)
  }
  
  save(est, file = "Zhang Estimates.rda")

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
 
  
######################################################### FIN PART 1 #############################################################


  
    

                                                  ## Part 2
                                #### Using in vitro RNA mixtures. ####
    
    
    
################################################################################################################################ #
## Setup ----

## Start!
  rm(list = ls()); gc()

## Directory
  ## Directory
  # root.dir <- "INSERT.YOUR.ROOT.DIRECTORY.PLEASE/"
  root.dir <- "/Volumes/Data1/PROJECTS/BrainCellularComposition/BrainCellularComposition_GitHub/"
  wd1 <- paste0(root.dir, "Data/") # working directory 1, for processing
  wd2 <- paste0(root.dir, "Results/BenchmarkingDatasets/RNA/") # and for analysis
  setwd(wd1)
  
## Functions and packages
  source("../Scripts/Fun_Composition.R")
  source("../Scripts/Fun_Preprocessing.R")
  source("../Scripts/Fun_Parameters.R")

## Annotation files
  load(file = "Preprocessed/geneInfo.rda")
  load("Preprocessed/exonicLength.rda")



############################################################################################################################### #
## This section reformats the RPKM data from three neuronal/astrocytic mixtures ----

## First, process the three mixture samples
    load("/Volumes/Data1/PROJECTS/Neurons_Mixture_experiment/RESULTS/STAR_Output/RPKM.rda") # Using premade STAR output... 
    rme <- rpkm[,c(2,5,6,7)] # These columns are my three mixtures!
    rme <- rme[,-1]
  
  ## Filter to protein coding genes
    # Load a file that classifies each annotation by the type of transcript
    load("/Volumes/Data0/PROJECTS/Gavin_2017/Deconvolution/Data/genes.gtf.rda") # This classifies each annotation
  
    # Extracting out the useful information...
    gtf_pcRNAselector <- gtf[grep("protein_coding", gtf[,2]),]
    rm(gtf) # To save memory
    geneInfo <- strsplit(as.character(gtf_pcRNAselector[,9]), split = ";")
    
    # Now I want the EnsIDs of all the pcRNA genes
    pcRNA.list <- sapply(geneInfo, `[`, 4) 
    pcRNA.list <- pcRNA.list[grep("gene_id", pcRNA.list)]
    pcRNA.list <- unique(pcRNA.list)
    pcRNA.list <- strsplit(pcRNA.list, split = " ")
    pcRNA.list <- sapply(pcRNA.list, `[`, 3)
    
    m <- which(rownames(rme) %in% pcRNA.list)
    
    rme.pcRNA <- rme[m,]
    
  ## Expression threshold
    # Min threshold: where the expression must be above 2rpkm in at least on reference
    rme.pcRNA.Tmin <- rme.pcRNA[which(apply(rme.pcRNA, 1, max) > 2),]
  
  ## Saving references
    save(rme.pcRNA.Tmin, file = "Raw/rme_pcRNA_Tmin.rda")
  
## Next, process the pure neuron and astrocyte samples
  ## 2-week differentiated neurons
    load("/Volumes/Data1/PROJECTS/Neurons_Mixture_experiment/RESULTS/STAR_Output/RPKM.rda") # Using premade STAR output... 
    IVdiffNeurons = rpkm[,c(1,2,3)]
    colnames(IVdiffNeurons)[3] = "Differentiated_neurons"
    IVdiffNeurons = as.data.frame(IVdiffNeurons)
  
  ## Cultured astrocytes, transfected with an empty vector. Mean RPKM across 3 biological replicates
    load("/Volumes/Data0/PROJECTS/circRNAs/RBFOX1_OE_EXPERIMENT_2/RESULTS/1_GeneCounts/STAR_geneCounts/rpkm.rda")
    IVAstro = rpkm[,c(6,7,8)]
    IVAstro = cbind(rpkm[,c(1,2)], rowMeans(IVAstro))
    colnames(IVAstro)[3] = "Astrocyte_mean"
    IVAstro = as.data.frame(IVAstro)
  
  ## Quick test to see if the order is the same between the Voineagu datasets
    m = match(rownames(IVAstro), rownames(IVdiffNeurons)); plot(m) # Yup, line is a slope of 1. The bit down the bottom represents all the NAs matching to the first NA   
  
  ## Building a single data frame that combines both references
    refIV = cbind(rownames(IVAstro), IVdiffNeurons$Differentiated_neurons, IVAstro$Astrocyte_mean)
    refIV = as.data.frame(refIV) 
    rownames(refIV) = refIV[,1]
    refIV = refIV[,-1]
    colnames(refIV) = c("Differentiated_neurons", "Astrocyte_mean")
    refIV$Differentiated_neurons = as.numeric(refIV$Differentiated_neurons)
    refIV$Astrocyte_mean = as.numeric(refIV$Astrocyte_mean)
  
  ## How to keep protein-coding genes only
    # Load a file that classifies each annotation by the type of transcript
    load("/Volumes/Data0/PROJECTS/Gavin_2017/Deconvolution/Data/genes.gtf.rda") # This classifies each annotation
    
    # Extracting out the useful information...
    gtf_pcRNAselector = gtf[grep("protein_coding", gtf[,2]),]
    rm(gtf) # To save memory
    geneInfo = strsplit(as.character(gtf_pcRNAselector[,9]), split = ";")
    
    # Now I want the EnsIDs of all the pcRNA genes
    pcRNA.list = sapply(geneInfo, `[`, 4) 
    pcRNA.list = pcRNA.list[grep("gene_id", pcRNA.list)]
    pcRNA.list = unique(pcRNA.list)
    pcRNA.list = strsplit(pcRNA.list, split = " ")
    pcRNA.list = sapply(pcRNA.list, `[`, 3)
    
    m = which(rownames(refIV)%in%pcRNA.list)
    refIV.pcRNA = refIV[m,]

  ## Save
    save(refIV.pcRNA, file = "Raw/rme_refIV_pcRNA.rda")

## Combine into a benchmarking dataset
  ## Expression threshold
    exp_thresh <- 1
    
  ## Load the RNA Mixture Experiment, in a partially-preprocessed state
    load("Raw/rme_pcRNA_Tmin.rda")
    load("Raw/rme_refIV_pcRNA.rda")
    
  ## Process
    common <- intersect(rownames(rme.pcRNA.Tmin), rownames(refIV.pcRNA))
    rme <- cbind(rme.pcRNA.Tmin[common,], refIV.pcRNA[common,])
    colnames(rme) <- c("n40", "n45", "n50", "n100", "n0")
    rme <- rme[,c("n0", "n40", "n45", "n50", "n100")]
    rme <- rme[which(apply(rme, 1, max) > exp_thresh),]
    
    rm(refIV.pcRNA, rme.pcRNA.Tmin)
    
  ## Dataframe of truth
    true_RME <- as.data.frame(matrix(nrow = 5, ncol = 2))
    rownames(true_RME) <- colnames(rme)
    colnames(true_RME) <- c("Neurons", "Astrocytes")
    true_RME$Neurons <- c(0, 0.4, 0.45, 0.5, 1)
    true_RME$Astrocytes <- 1 - true_RME$Neurons
  
  ## Reannotated version
    symbols <- addSymbol(rme)
  
  ## Save
    # rda
    save(rme, symbols, nFeatures, true_RME, file = "Preprocessed/rme.rda")
    
    # cibersort-compatible
    write.CIB(rme, "Preprocessed/CIB/rme.txt")


################################################################################################################################ #
## Create signatures  ----

## First, load signatures
# load("Raw/rme_refIV_pcRNA.rda")
# sigsRME$IH <- refIV.pcRNA; rm(refIV.pcRNA)
# colnames(sigsRME$IH) <- c("Neurons", "Astrocytes")  
# sigsRME$IH <- sigsRME$IH[which(apply(sigsRME$IH, 1, max) > exp_thresh),]
# # sigsBrain$IH <- sigsRME$IH[,c("Neurons", "Astrocytes")]
    
    
load("../../Data/Preprocessed/Signatures - Brain.rda")
sigsRNA <- c(sigsRNA["IH"], sigsBrain)
sigsRNA <- lapply(sigsRNA, function(x) {
  x <- x[,c("Neurons", "Astrocytes")]
  x[which(apply(x, 1, max) > 1),]
})

## Load IH






################################################################################################################################ #
## Begin deconvolution  ----

## I note here that MuSiC will not be run; this is because there is no corresponding single-cell/-nucleus data for the mixtures

## Set up lists 
  # estimates' list
  stats <- list(DRS = list(),
                CIB = list(),
                DTA = list(),
                xCell = list(),
                Blender = list(),
                Linseed = list())
  
  # stats list
  stats <- est


## Run algorithms with custom signatures
  for (j in names(sigsRNA)) { 
    print(paste0(j, ":", Sys.time()))
    est$DRS[[j]] <- run.DRS(rme, sigsRNA[[j]])
    est$CIB[[j]] <- run.CIB(from.file = FALSE,
                            sigObject = sigsRNA[[j]],
                            mixString = "rme.txt")
    est$DTA[[j]] <- run.DTA(rme, sigsRNA[[j]], alg = "diff", q = 0.01)
  }

  # CIB / LK fails, as the algorithm cannot find a "nu" above its threshold. below, I manually set it to a zero matrix for compatibility with scripts
  # est$CIB$LK <- est$CIB$F5; est$CIB$LK[,] <- 0

## Run algorithms with in-built signatures
  est$xCell <- run.xCell(symbols)
  est$Blender <- run.Blender(symbols)
  
  # filter to relevant cell-types in enrichment algorithms
  est$xCell <- lapply(est$xCell, function(x) x[,c("Neurons", "Astrocytes")])
  est$Blender <- lapply(est$Blender, function(x) x[,c("Neurons", "Astrocytes")])

## Run full deconvolution by Linseed
  est$Linseed <- linseed <- list()   

  # on all five samples
  pdf(file = "Linseed Plots.pdf", height = 2.5, width = 2.5)
  linseed$Full <- run.linseed(mixture = rme, nCelltypes = 2, write.plots = TRUE, write.data = TRUE)
  dev.off()
  
  est$Linseed$Full <- linseed$Full$Transformed
  colnames(est$Linseed$Full) <- c("Neurons", "Astrocytes") # manual relabelling of cell-types based on correlation
  
  pdf(file = "Linseed Plots (SVD).pdf", height = 2.5, width = 4)
  linseed$Full$Data$svdPlot() + labs(y = "Cumulative Variance Explained", x = "Number of Dimensions")
  dev.off()
  
  
  
## Save
  save(est, file = "RNA Mixtures Composition Estimates.rda") # load("RNA Mixtures Composition Estimates.rda")
  save(linseed, file = "Raw Linseed Data.rda")
  
  


################################################################################################################################ #
## Statistics  ----

## Stats
  # compute for deconvolution algorithms
  for (j in c("DRS", "DTA", "CIB")) {
    stats[[j]] <- lapply(est[[j]], function (x) { write.stats(true_RME, x, alg = j, error = TRUE) } )
  }
  
  # compute for xCell
  stats$xCell <- list()
  stats$xCell$Raw <- write.stats(true_RME, est$xCell$Raw, alg = "xCell.Raw", error = FALSE) 
  stats$xCell$Trans <- write.stats(true_RME, est$xCell$Transformed, alg = "xCell.Trans", error = FALSE) 
  
  # compute for Blender
  stats$Blender <- list()
  stats$Blender$AverageIndex <- write.stats(true_RME, est$Blender$AverageIndex, alg = "Blender.Average", error = FALSE) 
  stats$Blender$DarmanisIndex <- write.stats(true_RME, est$Blender$DarmanisIndex, alg = "Blender.Darmanis", error = FALSE) 
  
  # compute for Linseed  (full version)
  stats$Linseed <- list()
  stats$Linseed$Raw <- write.stats(true_RME, est$Linseed$Full$Raw, alg = "Linseed.Raw", error = TRUE) 
  stats$Linseed$Transformed <- write.stats(true_RME, est$Linseed$Full$Transformed, alg = "Linseed.Trans", error = TRUE) 
  
  # save
  save(stats, file = "Statistics (Revised).rda")
  
  # condense lists to single matrices (for easy export)
  for (j in names(stats)) {
    stats[[j]] <- do.call("cbind", data.frame(stats[[j]]))
    rownames(stats[[j]]) <- c("rho", "r", "rmse", "nrmse", "mae", "nmae")
  }
  
  for(j in names(stats)) write.csv(stats[[j]], file = paste0("Statistics - ", j, ".csv"))


################################################################################################################################ #
## Plot using matched / default signature ----  
  
## Partial deconvolution algorithms
  plot.list <- list()
  plot.list$CIB <- plot.scatter(t = true_RME, e = est$CIB$IH, ct = "Neurons", calcCor = "r", calcError = "nmae", colour = ct.colours[["Neurons"]], abline.colour = "black") + 
    labs(y = "Estimated NP", x = "True NP") +
    annotate("text", x = 0.88, y = 0.1, label = "CIB")
  plot.list$DRS <- plot.scatter(t = true_RME, e = est$DRS$IH, ct = "Neurons", calcCor = "r", calcError = "nmae", colour = ct.colours[["Neurons"]], abline.colour = "black") + 
    labs(y = "Estimated NP", x = "True NP") + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
    annotate("text", x = 0.88, y = 0.1, label = "DRS")
  plot.list$dtangle <- plot.scatter(t = true_RME, e = est$dtangle$IH, ct = "Neurons", calcCor = "r", calcError = "nmae", colour = ct.colours[["Neurons"]], abline.colour = "black") + 
    labs(y = "Estimated NP", x = "True NP") + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
    annotate("text", x = 0.82, y = 0.1, label = "dtangle")
  # plot.list$blank <- plot.empty
  pdf(file = "Signature-optimised Scatterplots, Deconvolution.pdf", height = 2, width = 4.5)
  plot_grid(plotlist = plot.list, ncol = 3, rel_widths = c(1, 0.8, 0.8))
  dev.off()
  
## Blender and xCell
  plot.list <- list()
  plot.list$bn <- plot.scatter(t = true_RME, e = est$Blender$AverageIndex, ct = "Neurons", calcCor = "r", calcError = FALSE, colour = ct.colours[["Neurons"]], abline = FALSE) + 
    scale_y_continuous(limits = c(NA, NA)) + 
    labs(x = "True NP", y = "Neuronal Enrichment") +
    annotate("text", x = 0.82, y = -0.5, label = "Blender")
  plot.list$xn <- plot.scatter(t = true_RME, e = est$xCell$Transformed, ct = "Neurons", calcCor = "r", calcError = FALSE, colour = ct.colours[["Neurons"]], abline = FALSE) + 
    labs(x = "True NP") +
    annotate("text", x = 0.88, y = 0.005, label = "xCell") +
    theme(axis.title.y = element_blank())
  
  plot.list$ba <- plot.scatter(t = true_RME, e = est$Blender$AverageIndex, ct = "Astrocytes", calcCor = "r", calcError = FALSE, colour = ct.colours[["Astrocytes"]], abline = FALSE) + 
    scale_y_continuous(limits = c(NA, NA)) + 
    scale_y_continuous(limits = c(NA, NA)) + 
    labs(x = "True AP", y = "Astrocytic Enrichment") +
    annotate("text", x = 0.25, y = -0.8, label = "Blender")
  plot.list$xa <- plot.scatter(t = true_RME, e = est$xCell$Transformed, ct = "Astrocytes", calcCor = "r", calcError = FALSE, colour = ct.colours[["Astrocytes"]], abline = FALSE) + 
    labs(x = "True AP") +
    annotate("text", x = 0.88, y = 1e-18, label = "xCell") +
    theme(axis.title.y = element_blank()) +
    scale_y_continuous(breaks = c(0, 5e-18, 1e-17, 1.5e-17, 2e-17), limits = c(0,2.1e-17))
  
  pdf(file = "Signature-optimised Scatterplots, Enrichment.pdf", height = 2, width = 7.5)
  plot_grid(plotlist = plot.list, ncol = 4, rel_widths = c(1, 0.9, 1, 0.9))
  dev.off()
  
## Linseed
  ## Full run
pdf(file = "Linseed Scatterplots.pdf", height = 3.5, width = 3.5)
plot.scatter(t = true_RME, e = est$Linseed$Full$Transformed, ct = "Neurons", calcCor = "r", calcError = "nmae", colour = ct.colours[["Neurons"]]) +
  labs(x = "True NP", y = "Linseed Cell-type 1")
dev.off()

  ## The run on the mixed samples 
  pdf(file = "Linseed Scatterplots (Mixed Samples).pdf", height = 3.5, width = 3.5)
  plot.scatter(t = true_RME[2:4,], e = est$Linseed$Mixed$Transformed, ct = "Neurons", calcCor = "r", calcError = "nmae", colour = ct.colours[["Neurons"]]) +
    labs(x = "True NP", y = "Linseed Cell-type 1") +
    scale_y_continuous(limits = c(0.3, 0.55)) +
    scale_x_continuous(limits = c(0.3, 0.55)) 
  dev.off()

  
################################################################################################################################ #
## Plot deconvolution using mismatched signatures  ----  

## Scatterplot   
  plot.list <- list()
  for(k in c("CIB", "DRS", "DTA")) {
    for(j in names(sigsRNA)) {
      
      annot.pos <- c(Inf, -Inf, 1, -0.5) 
      e <- est[[k]][[j]]
      
      plot.list[[j]] <- plot.scatter(t = true_RME, e = e, ct = "Neurons", ylab = "Estimated Proportion", calcCor = "r", 
                                     calcError = "nmae", colour = "black", abline.colour = "black", abline = TRUE, annot.pos = annot.pos)  +
        labs(title = j, x = "True Proportion") +
        scale_y_continuous(limits = c(0,1)) +
        scale_x_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) +
        theme(plot.title = element_text(hjust = 0.5), panel.grid.minor = element_blank())
      if (!(j %in% c("IH", "VL"))) plot.list[[j]] <- plot.list[[j]] + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
      if (!(j %in% names(sigsRNA)[6:10])) plot.list[[j]] <- plot.list[[j]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
      if (j == "IH") plot.list[[j]] <- plot.list[[j]] + labs(title = "Matched (IH)")
    }
    
    pdf(file = paste0("Origin Test - ", k, ", Scatterplot (Revised).pdf"), height = 4.5, width = 8)
    print(plot_grid(plotlist = plot.list, ncol = 5, rel_widths = c(1.2,1,1,1,1), rel_heights = c(1,1.2)))
    dev.off()
    
  }
 

######################################################### FIN PART 2 #############################################################
