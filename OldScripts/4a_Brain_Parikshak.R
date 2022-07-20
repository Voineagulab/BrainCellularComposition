################################################################################################################################ #
## Setup ----     


## Generic
rm(list = ls())
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Parikshak_revised/")
options(stringsAsFactors = FALSE)


## Load in Parikshak data
load("../../../Data/Preprocessed/Parikshak.rda")

## Filter as per original authors' analyses
  filtered <- parikshak.meta$Network.analysis..CTX | parikshak.meta$Network.Analysis..CB # list of samples in the filtered list (247 of 251)
  parikshak <- parikshak[,filtered]
  parikshak.counts <- parikshak.counts[,filtered]
  parikshak.meta <- parikshak.meta[filtered,]
  

## Load list of nuclear-stable genes
  load("../Lister_CrushedBrains/Stable Genes.rda")
  stable.genes <- stable.genes$notDE

  
## Load in signatures
  # load
  load("../../../Data/Preprocessed/Signatures - Brain.rda")
  load("../../../Data/Preprocessed/Signatures - Brain (Full Ct and Subct).rda")
  load("../../../Data/Preprocessed/Signatures - Brain (MuSiC).rda")
  load("../../../Data/Preprocessed/Multibrain.rda")
  
  ## Reprocess the brain signatures into various alternates
    # first, remove Endothelia from CA and LK, due to low numbers of representative cells
    sigsBrain$CA <- sigsBrain$CA[,-grep("End", colnames(sigsBrain$CA))]
    sigsBrain$LK <- sigsBrain$LK[,-grep("End", colnames(sigsBrain$LK))]
    
    # next, remove OPCs from Darmanis
    sigsBrain$DM <- sigsBrain$DM[,-grep("OPCs", colnames(sigsBrain$DM))]
  
    # for signatures with many celltypes, create simpler and complex version
    reduce.ct <- function(sig, x = sigsBrain, core.ct = c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelia")) {
      full.sig <- paste0(sig, "_")
      
      # create new level of signature with all ct 
      x[[full.sig]] <- x[[sig]][,-grep("Neurons|OPCs", colnames(x[[sig]]))] # remove neurons, as otherwise there are columns for exc/inh/neu
      
      # reduce the original sig to core ct
      x[[sig]] <- x[[sig]][,which(colnames(x[[sig]]) %in% core.ct)]
      
      return(x)
    }
    
    sigs.to.reduce <- c("VL", "NG", "CA", "TS", "LK")
    
    for (j in sigs.to.reduce) sigsBrain <- reduce.ct(j)

  
  ## For nuclear signatures, filter to stable genes
    stable.filter <- function(sig, x = sigsBrain) {
      x[[paste0(sig, ".")]] <- x[[sig]][which(rownames(x[[sig]]) %in% stable.genes),]   
      return(x)
    }
    
    sigs.to.stabilise <- c("VL", "NG", "CA", "LK", "VL_", "NG_", "CA_", "LK_")
  
    for (j in sigs.to.stabilise) sigsBrain <- stable.filter(j)
    
  ## Reprocess the full signatures
    names(sigsFull) <- paste0(names(sigsFull), "z")
    
    # create stable versions
    sigsFull <- stable.filter("VLz", x = sigsFull)
    sigsFull <- stable.filter("NGz", x = sigsFull)
    sigsFull <- stable.filter("LKz", x = sigsFull)
    sigsFull <- stable.filter("CAz", x = sigsFull)
    
    # add to sigsBrain
    sigsBrain <- c(sigsBrain, sigsFull)
  
    
  ## Reprocess the MuSiC signatures
    # this is now performed in the deconvolution loop
    
  ## Add multibrain
    sigsBrain$MB <- multibrain
    sigsBrain <- reduce.ct("MB")
    sigsBrain <- stable.filter("MB")
    
    sigsBrain$MB2 <- multibrain2
    
## Functions & Packages
source("../../../Scripts/Fun_Composition.R")
  
    
## Plotting parameters
  signature.plotting.order <- c("CA", "NG", "LK", "VL", "TS", "MM", "DM", "IP", "F5", "MB")    
  signature.annotation <- c("CA",
                            "NG",
                            "LK",
                            "VL",
                            "TS", 
                            "MM",
                            "DM",
                            "IP",
                            "F5",
                            "MB")
  
  # signature.colours <- c("steelblue4", "royalblue1", "darkturquoise", "dodgerblue",
  #                        "firebrick1", "tomato", "darkorange1", "tan1",
  #                        "grey50", "black")
  # names(signature.colours) <- signature.annotation
  
  signature.colours <- c("darkorange1", pal_lancet()(9))
  names(signature.colours) <- signature.plotting.order
      
## Filtering variables
  ctx <- which(parikshak.meta$RegionID != "vermis")
  deconv <- c("DRS", "DTA", "CIB", "MUS")
  enrich <- c("Blender", "Coex", "xCell")
  main.sigs <- c("CA", "NG", "LK", "VL", "TS", "MM", "DM", "IP", "F5", "MB2")    
  
  
################################################################################################################################ #
## Estimate composition ----  
  
## Setup
est <- list(DRS = list(),
            DTA = list(),
            CIB = list(),
            MUS = list(),
            xCell = list(),
            Blender = list(),
            Coex = list(),
            Linseed = list())


## Apply to all samples
  ## DeconRNASeq
    for (j in names(sigsBrain)) {
      est$DRS[[j]] <- run.DRS(parikshak, sigsBrain[[j]])
    }
  
  ## DTA
    for (j in names(sigsBrain)) {
      est$DTA[[j]] <- run.DTA(parikshak, sigsBrain[[j]], alg = "diff", q = 0.01)
    }  
  
  
  ## CIBERSORT
    for(j in names(sigsBrain)) {
      # report, as runtime is ~hours
      run.no <- which(names(sigsBrain) == j)
      print(paste0(j, ", run ", run.no, " of ", length(sigsBrain), ": ", Sys.time()))
      
      est$CIB[[j]] <- run.CIB(from.file = FALSE,
                              mixString = "Parikshak.txt",
                              sigObject =  sigsBrain[[j]])
      
      
      gc()
      
      save(est, file = "Composition Estimates (partial).rda")
    }
  
  ## MUSIC
    for(j in names(sigsMuSiC)) {
      print(j)
      Sys.time()
      
      ## Deconvolution with merged subtypes
          # celltypes to remove from the standard deconvolution
            remove <- which(sigsMuSiC[[j]]$meta$brain.ct2 == "OPCs") # remove OPCs
            
            if (j %in% c("CA", "LK")) {
              remove <- c(remove, grep("End", sigsMuSiC[[j]]$meta$brain.ct)) # remove endothelia from CA and LK signatures
            }
          
          
          
          # deconvolve with the core.ct
          gc()
          est$MUS[[j]] <- run.music(mixture = parikshak.counts, 
                                    use.meta.column = "brain.ct2", 
                                    signature = sigsMuSiC[[j]]$counts[,-remove], 
                                    signature.meta = sigsMuSiC[[j]]$meta[-remove,]) 
          
          # deconvolve with core celltypes, subsetted to stable genes
          gc()
          est$MUS[[paste0(j, ".")]] <- run.music(mixture = parikshak.counts, 
                                                 use.meta.column = "brain.ct2", 
                                                 signature = sigsMuSiC[[j]]$counts[which(rownames(sigsMuSiC[[j]]$counts) %in% stable.genes),remove], # removes non-stable genes
                                                 signature.meta = sigsMuSiC[[j]]$meta[-remove,]) 
          
          
          # deconvolve with subtypes
          gc()
          est$MUS[[paste0(j, "_")]] <- run.music(mixture = parikshak.counts, 
                                                 use.meta.column = "brain.ct",
                                                 signature = sigsMuSiC[[j]]$counts[,-remove], 
                                                 signature.meta = sigsMuSiC[[j]]$meta[-remove,])
          
          # deconvolve with subtypes, subsetted to stable genes
          est$MUS[[paste0(j, "_.")]] <- run.music(mixture = parikshak.counts, 
                                                 use.meta.column = "brain.ct",
                                                 signature = sigsMuSiC[[j]]$counts[which(rownames(sigsMuSiC[[j]]$counts) %in% stable.genes),-remove], 
                                                 signature.meta = sigsMuSiC[[j]]$meta[-remove,])
          
        ## Deconvolution with the full signature
          if (j %in% c("TS")) next # skip these signatures
          
          if (j == "CA") {
            keep <- which(sigsMuSiC$CA$meta$orig.celltype %in% colnames(sigsFull$CAz))
          } else {
            keep <- 1:nrow(sigsMuSiC[[j]]$meta)
          }
           # deconvolve with all celltypes in the publication
          gc()
          est$MUS[[paste0(j, "z")]] <- run.music(mixture = parikshak.counts, 
                                                 use.meta.column = "orig.celltype",
                                                 signature = sigsMuSiC[[j]]$counts[,keep], 
                                                 signature.meta = sigsMuSiC[[j]]$meta[keep,])
          
          # deconvolve with all celltypes in the publication, subsetted to stable genes
          est$MUS[[paste0(j, "z.")]] <- run.music(mixture = parikshak.counts, 
                                                 use.meta.column = "orig.celltype",
                                                 signature = sigsMuSiC[[j]]$counts[which(rownames(sigsMuSiC[[j]]$counts) %in% stable.genes),keep], 
                                                 signature.meta = sigsMuSiC[[j]]$meta[keep,])
          
      # free up RAM
      sigsMuSiC[[j]] <- sigsMuSiC[[j]][-1]
      
      # end loop
    }
  
  ## xCell
    est$xCell <- run.xCell(symbols)
  
  ## Blender
    est$Blender <- run.Blender(symbols)
  
  ## Linseed
    pdf(file = "Linseed Plots (full).pdf", height = 4, width = 4)
    est$Linseed <- run.linseed(mixture = gtex, 
                                      nCelltypes = 5, 
                                      write.plots = TRUE)
    dev.off()
  
  ## Co-expression
    coex <- list()
    coex$Full <- run.coex(mixture = gtex, 
                          signature = sigsBrain$SC.CA.IP, 
                          only.threshold = FALSE, 
                          sft = 18, # these data  do not reach scale-free topology by a reasonable power; set at 18, where r2 = 0.692, and median.k = 77
                          output.all = TRUE)
    y <- define.celltypes(x = coex$Full, algorithm.type = "Coex", celltype.signature = sigsBrain$SC.CA.IP, return.confidence = TRUE, nMarks.x = 10000, nMarks.sig = 100, bg = rownames(gtex))
    colnames(y$comp) <- y$assignments
    est$Coex <- y$comp

  
## Create a separate list for CTX estimates
  # collect all variance-independent methods!
  est.ctx <- lapply(est[1:4], function(x) {
    lapply(x, function(y) {
      y[ctx,]
    })
  })
  
  
  # run blender
  est.ctx$Blender <- run.Blender(symbols[,ctx])
  
  # run xCell (though it's not always affected by variance)
  est.ctx$xCell <- run.xCell(symbols[,ctx])

  # run Linseed on cortical samples
  pdf(file = "Linseed Plots (CTX).pdf", height = 4, width = 4)
  est.ctx$Linseed <- run.linseed(mixture = gtex[,ctx], 
                                 nCelltypes = 5, 
                                 write.plots = TRUE)
  dev.off()
  
  # run coex on cortical samples
  coex$CTX <- run.coex(mixture = gtex[,ctx], 
                       signature = sigsBrain$SC, 
                       only.threshold = FALSE, 
                       sft = "auto", 
                       output.all = TRUE)
   
  y <- define.celltypes(x = coex$CTX, algorithm.type = "Coex", celltype.signature = sigsBrain$SC.CA.IP, return.confidence = TRUE, nMarks.x = 10000, nMarks.sig = 100, bg = rownames(gtex))
  colnames(y$comp) <- y$assignments
  est.ctx$Coex <- y$comp
  
## Save
save(est, est.ctx, file = "Composition Estimates (Final).rda") # load("Composition Estimates.rda")
save(coex, file = "Raw co-expression data.rda")

## Output as a supplementary table
x <- est[-which(names(est) == "Linseed")]
x$xCell <- x$xCell$Transformed[,c("Neurons", "Astrocytes")]
x$Blender <- x$Blender$AverageIndex
x$Linseed <- x$Linseed$Transformed
x[deconv] <- lapply(x[deconv], function(x) {
  x <- x[which(names(x) %in% main.sigs)]
  names(x) <- gsub("2", "", names(x))
  return(x)
})
y <- sapply(x, function(x) do.call("cbind", x))
y <- do.call("cbind", y)
write.csv(y, file = "ST2_EstimatedComposition_SourceFiles_Parik.csv")


################################################################################################################################ #
## Evaluate using Goodness of Fit ----  


## Calculate GoF when data are log transformed, for equivalence to PsychEncode analyses
  gof <- list()
  for(j in deconv) {
    
    gof[[j]] <- list()
    
    for (k in names(est[[j]])) {
      
      if (k == "DM_") next # skip this run
      if (substr(k, 4, 4) == ".") next # skip deconvolution using subneurons + stable
      if (substr(k, 3, 3) == "z") next # skip deconvolution using all ct
      
      l <- gsub("\\.", "", k) # this means that if the stable deconvolution is used, GoF is calculated using a non-stable signature,
      
      print(paste0(Sys.time(),": ", j, "_", k))
      gof[[j]][[k]] <- write.gof.v2(measuredExp = parikshak, 
                                 estimatedComp = est[[j]][[k]], 
                                 signatureUsed = sigsBrain[[l]])
    }
  }  


## Save
save(gof, file = "Goodness of Fit (Final) (V2 with l).rda")  # load("Goodness of Fit (Final) (V2).rda")


## Plot distributions
  ## Setup and annotate data  
    plot.data <- melt(gof)
    
    # rename columns
    colnames(plot.data) <- c("Metric", "value", "Signature", "Algorithm")
    
    # just correlation as the GoF metric
    plot.data <- plot.data[which(plot.data$Metric == "r"),] # focus on correlation
    
    # add region information
    plot.data$Region <- parikshak.meta$RegionID # does this work?
    plot.data$Region[grep("ba", plot.data$Region)] <- "Cortex"
    plot.data$Region[grep("vermis", plot.data$Region)] <- "Cerebellum"
    plot.data$Region <- factor(plot.data$Region, levels = c("Cortex", "Cerebellum"))
    
    # remove stable/allct and stable/full combined deconvolutions
    # plot.data <- plot.data[-grep("_\\.", plot.data$Signature),]
    
    # add deconvolution filtering information
  plot.data$Filters <- "Standard"
  plot.data$Filters[which(substr(plot.data$Signature, 3, 3) == ".")] <- "Stable"
  plot.data$Filters[which(substr(plot.data$Signature, 3, 3) == "_")] <- "Subneurons"
  # plot.data$Filters[which(substr(plot.data$Signature, 3, 4) == "_.")] <- "Stable Subneurons"
  # plot.data$Filters[which(substr(plot.data$Signature, 3, 3) == "z")] <- "All Subtypes"
  # plot.data$Filters[which(substr(plot.data$Signature, 3, 4) == "z.")] <- "Stable All Subtypes"
  # plot.data$Filters <- factor(plot.data$Filters, levels = c("None", "Stable", 
  #                                                           "Subneurons", "Stable Subneurons", 
  #                                                           "All Subtypes", "Stable All Subtypes")) 
    plot.data$Filters <- factor(plot.data$Filters, levels = c("Standard", "Stable", 
                                                            "Subneurons")) 
    
    
     plot.data <- plot.data[-which(plot.data$Signature %in% c("MB", "MB_", "MB.")),]
    plot.data$Signature <- gsub("2", "", plot.data$Signature)
    plot.data$Signature <- substr(plot.data$Signature, 1, 2) # can remove this information from this column now
    
    # add a tag to identify paired deconvolution runs
    plot.data$id <- paste0(plot.data$Signature, "_", plot.data$Algorithm)
    
    # add colour information
    plot.data$Signature <- factor(plot.data$Signature, levels = signature.plotting.order)
    levels(plot.data$Signature) <- signature.annotation
     
 ## Plot!
    violins <- function(alg, region = "Cortex", y = c(0.3, 0.8)) {
      
      ggplot(plot.data[which(plot.data$Filters == "Standard" & plot.data$Algorithm == alg & plot.data$Region == region),], aes(x = Signature, y = value, colour = Signature, fill = Signature)) +
        geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.9) +
        geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
        theme_bw() +
        scale_colour_manual(values = signature.colours) +
        scale_fill_manual(values = signature.colours) +
        scale_y_continuous(limits = y, expand = c(0,0)) +
        theme(panel.grid = invis, legend.position = "none") +
        geom_hline(yintercept = c(0.5, 0.7), linetype = 2) +
        labs(y = "Goodness of Fit (r)")
    }
    
  # the first plot is a comparison of the goodness of fit under basic conditions (core ct, all genes)
  # it contrasts performance across signatures and regions within an algorithm 
  pdf(file = "Final Goodness of Fit (Main).pdf", height = 3, width = 4)
  violins("CIB")
  dev.off()
  
  pdf(file = "Final Goodness of Fit (Supplementary, Other Algs).pdf", height = 3, width = 8)
  plot_grid(violins("DTA") + labs(title = "DTA"),
            violins("DRS") + labs(title = "DRS") + theme(axis.title.y = invis, axis.text.y = invis), 
            violins("MUS") + labs(title = "MUS") + theme(axis.title.y = invis, axis.text.y = invis), 
            nrow = 1, rel_widths = c(1,0.8,0.6))
  dev.off()
  
  pdf(file = "Final Goodness of Fit (Supplementary, All Algs and Regions).pdf", height = 6, width = 8)
  ggplot(plot.data[which(plot.data$Filters == "Standard"),], aes(x = Signature, y = value, colour = Signature, fill = Signature)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.9) +
    geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
    theme_bw() +
    facet_grid(Region~Algorithm) +
    scale_colour_manual(values = signature.colours) +
    scale_fill_manual(values = signature.colours) +
    theme(panel.grid = invis, legend.position = "bottom") +
    geom_hline(yintercept = c(0.5, 0.7), linetype = 2) +
    labs(y = "Goodness of Fit (r)")
  dev.off()
    
  
    # this next plot explores how stable gene selection impacts performance, comparing otherwise like-for-like deconvolutions
  p <- plot.data[which(substr(plot.data$Signature, 1, 2) %in% sigs.to.stabilise),]

  pdf(file = "Final Goodness of Fit (Alternate Deconvolution Strategies).pdf", height = 6, width = 8)
  ggplot(p, aes(x = substr(Signature,1,2), y = value, colour = Filters, fill = Filters)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.8) +
    geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
    facet_grid(Region~Algorithm) +
    theme_bw() +
    # scale_colour_manual(values = c("black", "darkorange1")) +
    # scale_fill_manual(values = c("black", "darkorange1")) +
    theme(panel.background = element_rect(fill = "grey90"), panel.grid = invis) +
    geom_hline(yintercept = c(0.5, 0.7), linetype = 2) +
    theme(legend.position = "bottom") +
    labs(y = "Goodness of Fit (r)", x = "Signature")
  dev.off()
    
   

## Plot ranks!
    ## Summary plot
      plot.ranks <- function(subset, title, text.colour = "black", subneurons = FALSE) {
        ranks <- lapply(gof, function(x) { sapply(x, function(y) median(y$r[subset])) })
        
        if (subneurons) {
          ranks <- sapply(ranks, function(x) x[substr(names(x), 3, 4) == "_"])
          
        } else {
          ranks <- sapply(ranks, function(x) x[-grep("\\.|_|z|MB$", names(x))])
        }
        
        
        ranks <- do.call("c", ranks)
        names(ranks) <- gsub("_", "", names(ranks))
        names(ranks) <- gsub("2", "", names(ranks))
        
        plot.data <- data.frame(Median = ranks,
                            Algorithm = substr(names(ranks), 1, 3),
                            Signature = substr(names(ranks), 5, 6))
        
        plot.data$Rank <- rank(-plot.data$Median)
        plot.data$Label <- paste0(signif(plot.data$Median, 3), "\n(", plot.data$Rank, ")")
        
        plot.data$Signature <- factor(plot.data$Signature, levels = signature.plotting.order)

        return(ggplot(plot.data, aes(x = Signature, y = Algorithm, fill = as.factor(Rank), label = Label)) +
                 geom_tile(colour = "black") +
                 scale_fill_viridis(discrete = TRUE, direction = -1, alpha = 0.8) +
                 geom_text(size = 2.5, colour = text.colour) +
                 theme_bw() +
                 labs(title = title) +
                 theme(legend.position = "none", axis.title = element_blank(), panel.border = element_blank(),
                       plot.title = element_text(hjust = 0.5), panel.grid = element_blank()))
      }


      pdf(file = "Final Goodness of Fit (Ranks) (r).pdf", height = 3, width = 3.75)
      plot.ranks(subset = ctx, title = "Cortex")
      plot.ranks(subset = c(1:ncol(parikshak))[-ctx], title = "Cerebellum")
      dev.off()

      
## Statistical tests on the effect of stable genes
  ## Function
    compare.stable.gof <- function(alg, sig) {
      x <- gof[[alg]]
      
      ctx.improve <- sign(mean(x[[sig]]$r[ctx]) - mean(x[[paste0(sig, "_")]]$r[ctx])) < 0
      ctx.effect <- mean(x[[sig]]$r[ctx]) - mean(x[[paste0(sig, "_")]]$r[ctx])
      ctx.p <- wilcox.test(x[[sig]]$r[ctx], x[[paste0(sig, "_")]]$r[ctx], paired = TRUE)$p.value
      
      cb.improve <- sign(mean(x[[sig]]$r[-ctx]) - mean(x[[paste0(sig, "_")]]$r[-ctx])) < 0
      cb.effect <- mean(x[[sig]]$r[-ctx]) - mean(x[[paste0(sig, "_")]]$r[-ctx])
      cb.p <- wilcox.test(x[[sig]]$r[-ctx], x[[paste0(sig, "_")]]$r[-ctx], paired = TRUE)$p.value
      
      output <- data.frame(Improved = c(ctx.improve, cb.improve),
                           MeanImprovement = c(ctx.effect, cb.effect),
                           P = c(ctx.p, cb.p))
      rownames(output) <- c("ctx", "cb") 
      
      return(output)
    }
  
  ## Apply function
    stats.gof <- list()
    for (j in names(gof)) {
      stats.gof[[j]] <- list()
      for (k in sigs.to.stabilise[1:4]) {
        stats.gof[[j]][[k]] <- compare.stable.gof(alg = j, sig = k)
      }
      stats.gof[[j]] <- do.call("rbind", stats.gof[[j]] )
    }
  
  ## Output  
    stats.gof <- do.call("rbind", stats.gof)
    write.csv(stats.gof, file = "Paired Wilcox Tests, Stable Genes.csv")
  
  
  
  
  
 
  
################################################################################################################################ #
## Trends in composition ----  

## NP vs. PCA
 # ctx samples
  pca <- princomp(log2(gtex[,ctx] + 0.5), cor = TRUE)
  e <- lapply(est.ctx, function(x) do.call("cbind", x))
  # e <- sapply(e, function(x) x[,grep("Neurons", colnames(x))])
  e <- do.call("cbind", e)
  
  cors <- list()
  for(j in colnames(est$CIB$IP)) { # that is, for each of the five celltypes
    a <- e[,grep(j, colnames(e))]
    b <- data.frame(PC1 = as.numeric(cor(pca$loadings[,1], a, method = "s")),
                    PC2 = as.numeric(cor(pca$loadings[,2], a, method = "s")))
    rownames(b) <- gsub(pattern = paste0(".", j), replacement = "", x = colnames(a))
    cors[[j]] <- b
  }
  
  
  all.runs <- rownames(cors$Neurons)
  
  
  cors <- lapply(cors, function(x) {
    dat <- as.data.frame(matrix(nrow = length(all.runs), ncol = 2))
    colnames(dat) <- c("PC1", "PC2")
    rownames(dat) <- all.runs
    
    m <- match(rownames(dat), rownames(x))
    
    dat$PC1 <- x$PC1[m]
    dat$PC2 <- x$PC2[m]
    dat
  })
  
  cors <- do.call("cbind", cors)
  write.csv(cors, file = "../../SuppTables/STY_GTExCTX_allCTvsPCA.csv")
  
  # # plot CIB IP
  # plot.data <- data.frame(est.ctx$CIB$IP, PC1 = pca$loadings[,1])
  # plot.data <- melt(plot.data, id.vars = "PC1")
  # ggplot(plot.data, aes(x = PC1, y = value, col = variable)) +
  #   geom_point() +
  #   geom_smooth(method = "lm")
  # 
  # plot.data <- data.frame(est.ctx$CIB$IP, PC2 = pca$loadings[,2])
  # plot.data <- melt(plot.data, id.vars = "PC2")
  # ggplot(plot.data, aes(x = PC2, y = value, col = variable)) +
  #   geom_point() +
  #   geom_smooth(method = "lm", se = FALSE) +
  #   labs(y = "Celltype Proportion") +
  #   scale_colour_manual(values = ct.colours)

## Boxplot proportions
  ## Part 1: All deconvolution algorithm combinations
  plot.data <- lapply(est[deconv], function(x) {
    x <- x[main.sigs]
    lapply(x, function(y) { 
      y$Region <- "Cerebellum" 
      y$Region[ctx] <- "Cortex"
      return(y)
    })
  })
  
  plot.data <- melt(plot.data)
  plot.data <- plot.data[,-5]
  
  colnames(plot.data) <- c("Region", "Celltype", "Proportion", "Signature", "Algorithm")
  plot.data$Proportion <- as.numeric(plot.data$Proportion)
  plot.data$Signature <- gsub("2", "", plot.data$Signature)
  
  pdf(file = "Composition Distributions, Deconv, V2.pdf", height = 5.5, width = 8)
  for (j in names(table(plot.data$Celltype))) {
    print(ggplot(plot.data[which(plot.data$Celltype == j),], aes(x = Signature, y = Proportion, colour = Signature, fill = Signature)) +
            geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.7) +
            geom_boxplot(fill = "white", outlier.shape = NA, width = 0.3, position = position_dodge(width = 0.8)) + # note that, on other plots, width = 0.2
            theme_bw() +
            theme(panel.background = element_rect(fill = "grey90"), axis.text.x = element_text(angle = 90),
                  legend.position = "none") +
            geom_hline(yintercept = c(0.5), linetype = 2) +
            # facet_wrap(~Algorithm, ncol = 1) +
            facet_grid(Region~Algorithm) +
            scale_y_continuous(limits = c(0,1)) +
            scale_colour_manual(values = signature.colours) +
            scale_fill_manual(values = signature.colours) +
            labs(y = "Estimated Proportion", title = j))
  }
  dev.off()
    
    # 
    # pdf(file = "Composition Distributions, Deconv.pdf", height = 11, width = 7.5)
    # ggplot(plot.data, aes(x = Combination, y = Proportion, colour = Region, fill = Region)) +
    #   geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.7) +
    #   geom_boxplot(fill = "white", outlier.shape = NA, width = 0.3, position = position_dodge(width = 0.8)) + # note that, on other plots, width = 0.2
    #   theme_bw() +
    #   theme(panel.background = element_rect(fill = "grey90"), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90),
    #         legend.position = "bottom") +
    #   geom_hline(yintercept = c(0.5), linetype = 2) +
    #   facet_wrap(~Celltype, ncol = 1) +
    #   scale_y_continuous(limits = c(0,1)) +
    #   scale_colour_manual(values = reg.colours) +
    #   scale_fill_manual(values = reg.colours) +
    #   labs(y = "Estimated Proportion")
    # dev.off()
  
  ## Part 2: All enrichment algorithms
   #  plot.data <- list(Blender = est$Blender$AverageIndex,
   #                    xCell = est$xCell$Transformed[c("Neurons", "Astrocytes")],
   #                    Coex = est$Coex)
   #  
   #  plot.data <- lapply(plot.data, function(x) {
   #    x$Region <- parikshak.meta$RegionID
   #    return(x)
   #  })
   #  
   #  plot.data <- melt(plot.data)
   #  colnames(plot.data) <- c("Region", "Celltype", "Enrichment", "Algorithm")
   # 
   #   plot.data$Region[grep("ba", plot.data$Region)] <- "CTX"
   #  plot.data$Region[grep("vermis", plot.data$Region)] <- "CB"
   #  
   # pdf(file = "Composition Distributions, Enrichment, V2.pdf", height = 2, width = 7.5)
   # for (j in names(table(plot.data$Celltype))) { 
   #  print(ggplot(plot.data[which(plot.data$Celltype == j),], aes(x = Celltype, colour = Region, fill = Region, y = Enrichment)) +
   #    geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.7) +
   #    geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
   #    theme_bw() +
   #      labs(title = j) +
   #    theme(panel.background = element_rect(fill = "grey90"), axis.title.x = element_blank(), axis.text.x = element_blank(),
   #          legend.position = "none") +
   #    facet_wrap(~Algorithm, ncol = 3, scale = "free_y"))
   # }
   #  dev.off()
   #  
  
  # ## Part 3: Boxplot the best performing algorithm
  #   # from best[length(best)], plot CIB_IP.CA, but since that's not being considered for space reason, go for the second best: CIB_SC.CA.IP
  #   plot.data <- est$CIB$SC.CA.IP
  #   plot.data$Region <- gtexMeta$BroadRegion
  #   plot.data <- melt(plot.data)
  #   colnames(plot.data) <- c("Region", "Celltype", "Proportion")
  #   
  #   pdf(file = "Composition Distributions, Best Algorithm.pdf", height = 2.5, width = 7.5)
  #   ggplot(plot.data, aes(x = Celltype, y = Proportion, colour = Region, fill = Region)) +
  #     geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.7) +
  #     geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
  #     theme_bw() +
  #     theme(panel.background = element_rect(fill = "grey90"), axis.title.x = element_blank(), axis.text.x = element_text(colour = ct.colours, face = "bold")) +
  #     scale_y_continuous(limits = c(0,1)) +
  #     scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
  #     scale_colour_manual(values = reg.colours) +
  #     scale_fill_manual(values = reg.colours) +
  #     labs(y = "Estimated Proportion")
  #   dev.off()
   
  
## Correlations in cell-type proportion across runs!
  # function
  plot.heatmap <- function(celltype, h = 5.5, w = 8) {
    # collect composition estimates from the ctx for the desired celltype
    cors <- lapply(est.ctx, function(x) {
      if("TS" %in% names(x)) {
        x <- x[which(names(x) %in% main.sigs)]
      }
      x <- do.call("cbind", x)
      x <- x[,grep(celltype, colnames(x))]
      return(x)
    })
    
    cors <- do.call("cbind", cors)
    colnames(cors) <- sapply(strsplit(colnames(cors), paste0("\\.", celltype)), "[", 1)
    
    # filter out certain runs
    remove <- c(grep("xCell.Raw", colnames(cors)),
                grep("Blender.DarmanisIndex", colnames(cors)))
    if(celltype == "Endothelia") remove <- c(remove, grep("xCell", colnames(cors)))
    cors <- cors[-remove, -remove]
    
    # rename run names
    colnames(cors) <- gsub(".Transformed", "", colnames(cors))
    colnames(cors) <- gsub(".AverageIndex", "", colnames(cors))
    colnames(cors) <- gsub("2", "", colnames(cors))
    colnames(cors) <- gsub("\\.", " / ", colnames(cors))
    
    # convert to correlation matrix
    cors <- cor(cors, method = "s")
    
    # check for an mark NAs
    n <- which(is.na(cors[,1]))
    if (length(n) > 0) { cors[n,] <- cors[,n] <- -1.5 }
    
    # plot
    colours <- c("black", # for NA, when there's no variance in cell-type proportion
                 "midnightblue", # for -1 < cor < 0
                 viridis_pal()(10))
    
    pdf(file = paste0("Correlation Heatmap V2 ", celltype,".pdf"), height = h, width = w)
    heatmap.2(cors, dendrogram = "row", trace = "none",  density.info = "none", 
          key.title = "NA", key = FALSE, keysize = 2.5, key.xlab = "Spearman", 
          main = NULL, col = colours, breaks = c(-2,-1,seq(0, 1, 0.1)), na.color = "red", na.rm = FALSE,
          cexRow = 0.8, cexCol = 0.8, lwid = c(2,10), lhei = c(2,20), mar = c(4,4))
    dev.off()
  }
  
  
  plot.heatmap("Astrocytes")
  plot.heatmap("Neurons")
  plot.heatmap("Endothelia")
  plot.heatmap("Microglia")
  plot.heatmap("Oligodendrocytes")
  

  
################################################################################################################################ #
## Evaluate signature-free methods ----  

## Evaluate Linseed
  pdf(file = "Linseed Correlation Plot, Deconvolved CTX.pdf", height = 2.5, width = 3.5)
    plot.data <- est.ctx$Linseed$Transformed
    colnames(plot.data) <- paste0("Linseed", 1:5)
    plot.data <- as.data.frame(cor(est.ctx$CIB$IP, plot.data, method = "p"))
    plot.data <- round(plot.data, 2)
    plot.data$ct <- substr(rownames(plot.data), start = 1, stop = 3)
    plot.data <- melt(plot.data)
    colnames(plot.data) <- c("ct", "linseed", "r")
    plot.data$ct <- factor(plot.data$ct, levels = c("Neu", "Ast", "Oli", "Mic", "End"))

    # plot
    print(ggplot(plot.data, aes(x = ct, y = linseed)) +
      geom_tile(aes(fill = r), colour = "black", width = 0.95, height = 0.95) +
      geom_text(aes(label = r), colour = "black", na.rm = TRUE, size = 3) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1,1)) +
      theme_bw() +
      scale_x_discrete(expand = c(0,0)) +
      scale_y_discrete(expand = c(0,0)) +
      # labs(x = "Brain Cell-type Proportion", y = "Linseed Cell-type Proportion") +
      theme(panel.border = element_blank(), axis.title = element_blank(), axis.text.x = element_text(colour = ct.colours, face = "bold")) )
  dev.off()

  pdf(file = "Linseed Correlation Plot, Deconvolved All.pdf", height = 2.5, width = 3.5)
  plot.data <- est$Linseed$Transformed[ctx,]
  colnames(plot.data) <- paste0("Linseed", 1:5)
  plot.data <- as.data.frame(cor(est.ctx$CIB$IP, plot.data, method = "p"))
  plot.data <- round(plot.data, 2)
  plot.data$ct <- substr(rownames(plot.data), start = 1, stop = 3)
  plot.data <- melt(plot.data)
  colnames(plot.data) <- c("ct", "linseed", "r")
  plot.data$ct <- factor(plot.data$ct, levels = c("Neu", "Ast", "Oli", "Mic", "End"))
  
  # plot
  print(ggplot(plot.data, aes(x = ct, y = linseed)) +
          geom_tile(aes(fill = r), colour = "black", width = 0.95, height = 0.95) +
          geom_text(aes(label = r), colour = "black", na.rm = TRUE, size = 3) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1,1)) +
          theme_bw() +
          scale_x_discrete(expand = c(0,0)) +
          scale_y_discrete(expand = c(0,0)) +
          # labs(x = "Brain Cell-type Proportion", y = "Linseed Cell-type Proportion") +
          theme(panel.border = element_blank(), axis.title = element_blank(), axis.text.x = element_text(colour = ct.colours, face = "bold")) )
  dev.off()
  
