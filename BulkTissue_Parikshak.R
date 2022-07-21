################################################################################################################################ #
## Setup ----     


## Start!
  rm(list = ls()); gc()

## Directory
  # root.dir <- "INSERT.YOUR.ROOT.DIRECTORY.PLEASE/"
  root.dir <- "/Volumes/Data1/PROJECTS/BrainCellularComposition/BrainCellularComposition_GitHub/"
  wd1 <- paste0(root.dir, "Data/") # working directory 1, for processing
  wd2 <- paste0(root.dir, "Results/BulkTissue/Parikshak/") # and for analysis
  setwd(wd1)
  
## Functions and packages
  source("../Scripts/Fun_Composition.R")
  source("../Scripts/Fun_Preprocessing.R")
  source("../Scripts/Fun_Parameters.R")

  # please download these files from the GitHub, in ./Files
  load("Preprocessed/exonicLength.rda") # exonic lengths of genes, used for gene length normalisation
  load("Preprocessed/geneInfo.rda") # information about gene symbols, ids, and type
  
## Load signatures
  load("Preprocessed/Signatures - Brain.rda")
  load("Preprocessed/Signatures - Brain (MuSiC).rda")
  load("Preprocessed/Multibrain.rda")
  
################################################################################################################################ #
## Process expression data ----
  
  
## Data from Parikshak et al, 2016
  load("Raw/ASD_RNAseq_ExpressionData.Rdata")
  
## Data from Parikshak et al, 2016
  parikshak <- datExpr.HTSC.unionexon
  parikshak.meta <- datMeta
  rm(datExpr.Cufflinks, datExpr.HTSC.unionexon, datExpr.HTSC.wholegene, datLB.gene, datMeta, datOK.gene)
  
## RPKM normalisation
  parikshak.counts <- parikshak
  parikshak <- rpkm(parikshak)
  
## Expression thresholding
  exp_thresh <- 1 # minimum rpkm cutoff
  NAs <- which(apply(parikshak, 1, anyNA))
  if (length(NAs) > 0) parikshak <- parikshak[-which(apply(parikshak, 1, anyNA)),] # first, remove any rows with an NA
  n <- round(min(table(parikshak.meta$RegionID)) / 2) # a gene must be expressed above exp_thresh in 40 samples (half the number of samples in the least-represented region)
  parikshak <- parikshak[rowSums(parikshak >= exp_thresh) >= n , ] # here, pRPKM >= exp_thresh is a logical matrix, which we need to be TRUE in >= 41 sample

## Create version with symbol-based annotation
  parikshak.symbol <- addSymbol(parikshak)

## Create version of count data with the same rows
  parikshak.counts <- parikshak.counts[rownames(parikshak),]
  
## Augment metadata. 
  # calculate summaries of the sequencing statistics per Parikshak's analyses. Code taken from dhglab's GitHub
  seqInfo <- data.matrix(parikshak.meta[,c(25:43)])
  seqInfo <- seqInfo[,c(6:19,c(5))]
  seqInfo <- t(scale(seqInfo,scale=F))
  PC.seqInfo <- prcomp(seqInfo);
  varexp <- (PC.seqInfo$sdev)^2 / sum(PC.seqInfo$sdev^2)
  topPC.seqInfo <- PC.seqInfo$rotation[,1:2];
  colnames(topPC.seqInfo) <- c("SeqSV1","SeqSV2") ## Recompute since those in datMeta were with additional samples included
  parikshak.meta$seqStatPC1 <- as.numeric(topPC.seqInfo[,1])
  parikshak.meta$seqStatPC2 <- as.numeric(topPC.seqInfo[,2])
  
## Filter as per original authors' analyses
  filtered <- parikshak.meta$Network.analysis..CTX | parikshak.meta$Network.Analysis..CB # list of samples in the filtered list (247 of 251)
  parikshak <- parikshak[,filtered]
  parikshak.counts <- parikshak.counts[,filtered]
  parikshak.meta <- parikshak.meta[filtered,]
  
## Save
  # rda
  save(parikshak, parikshak.symbol, parikshak.counts, parikshak.meta, file = "Preprocessed/Parikshak.rda")
  
  # cibersort-compatible files
  write.CIB(data = parikshak, dir = "Preprocessed/CIB/Parikshak.txt") 

  
################################################################################################################################ #
## Setup for deconvolution ----

  
## Load in signatures
  # already done above

## Reprocess the brain signatures into various alternates
  # for signatures with many celltypes, create simpler and complex version
  reduce.ct <- function(sig, x = sigsBrain, core.ct = c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelia")) {
    
    # reduce the original sig to core ct
    x[[sig]] <- x[[sig]][,which(colnames(x[[sig]]) %in% core.ct)]
    
    return(x)
  }
  
  sigs.to.reduce <- c("DM", "VL", "NG", "CA", "TS", "LK")
  
  for (j in sigs.to.reduce) sigsBrain <- reduce.ct(j)


## Reprocess the MuSiC signatures
  # this is now performed in the deconvolution loop
  
## Add multibrain
  sigsBrain$MB <- multibrain
  

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
  
  library(ggsci)
  signature.colours <- c("darkorange1", pal_lancet()(9))
  names(signature.colours) <- signature.plotting.order
      
## Filtering variables
  ctx <- which(parikshak.meta$RegionID != "vermis")
  deconv <- c("DRS", "DTA", "CIB", "MUS")
  enrich <- c("Blender", "Coex", "xCell")
  main.sigs <- c("CA", "NG", "LK", "VL", "TS", "MM", "DM", "IP", "F5", "MB")    
  
  
################################################################################################################################ #
## Estimate composition ----  
  
  
setwd(wd2)
  
## Setup
  est.parik <- list(DRS = list(),
                    DTA = list(),
                    CIB = list(),
                    MUS = list(),
                    xCell = list(),
                    Blender = list(),
                    Coex = list())


## Apply to all samples
  ## DeconRNASeq
    for (j in names(sigsBrain)) {
      est.parik$DRS[[j]] <- run.DRS(parikshak, sigsBrain[[j]])
    }
  
  ## DTA
    for (j in names(sigsBrain)) {
      est.parik$DTA[[j]] <- run.DTA(parikshak, sigsBrain[[j]], alg = "diff", q = 0.01)
    }  
  
  
  ## CIBERSORT
    for(j in names(sigsBrain)) {
      # report, as runtime is ~hours
      run.no <- which(names(sigsBrain) == j)
      print(paste0(j, ", run ", run.no, " of ", length(sigsBrain), ": ", Sys.time()))
      
      est.parik$CIB[[j]] <- run.CIB(from.file = FALSE,
                              mixString = "Parikshak.txt",
                              sigObject =  sigsBrain[[j]])
      
      save(est.parik, file = "Temporary Estimates.rda") # in case of crash
    }
  
  ## MuSiC
    for(j in names(sigsMuSiC)) {
      print(j)
      Sys.time()
      
      ## Deconvolution with merged subtypes
          # celltypes to remove from the standard deconvolution
          remove <- which(sigsMuSiC[[j]]$meta$brain.ct2 == "OPC") # remove OPCs
            

          gc()
          est.parik$MUS[[j]] <- run.music(mixture = parikshak.counts, 
                                    use.meta.column = "brain.ct2", 
                                    signature = sigsMuSiC[[j]]$counts[,-remove], 
                                    signature.meta = sigsMuSiC[[j]]$meta[-remove,]) 
      

      # end loop
    }
  
  ## xCell
    est.parik$xCell <- run.xCell(symbols)
  
  ## Blender
    est.parik$Blender <- run.Blender(symbols)
  

  ## Co-expression
    coex <- list()
    coex$Full <- run.coex(mixture = gtex, 
                          signature = sigsBrain$Multibrain, 
                          only.threshold = FALSE, 
                          sft = 18, # these data  do not reach scale-free topology by a reasonable power; set at 18, where r2 = 0.692, and median.k = 77
                          output.all = TRUE)
    y <- define.celltypes(x = coex$Full, algorithm.type = "Coex", celltype.signature = sigsBrain$Multibrain, return.confidence = TRUE, nMarks.x = 10000, nMarks.sig = 100, bg = rownames(gtex))
    colnames(y$comp) <- y$assignments
    est.parik$Coex <- y$comp


    
    
################################################################################################################################ #
## Estimate composition in cortical samples ----  
    
## For some algorithms, the composition estimates rely on variance across samples, thus subsetting samples will create different estimates
    
## First, collect all variance-independent methods!
  est.parik.ctx <- lapply(est.parik[1:4], function(x) {
    lapply(x, function(y) {
      y[ctx,]
    })
  })
  
## Rerun on variance-dependent methods: Blender, xCell, and Coex
  # run blender
  est.parik.ctx$Blender <- run.Blender(symbols[,ctx])
  
  # run xCell (though it's not always affected by variance)
  est.parik.ctx$xCell <- run.xCell(symbols[,ctx])

  # run coex on cortical samples
  coex$CTX <- run.coex(mixture = gtex[,ctx], 
                       signature = sigsBrain$SC, 
                       only.threshold = FALSE, 
                       sft = "auto", 
                       output.all = TRUE)
   
  y <- define.celltypes(x = coex$CTX, algorithm.type = "Coex", celltype.signature = sigsBrain$SC.CA.IP, return.confidence = TRUE, nMarks.x = 10000, nMarks.sig = 100, bg = rownames(gtex))
  colnames(y$comp) <- y$assignments
  est.parik.ctx$Coex <- y$comp
  
################################################################################################################################ #
## Save ----  
  
save(est.parik, est.parik.ctx, file = "Composition Estimates.rda") 


################################################################################################################################ #
## Evaluate using Goodness of Fit ----  


## Calculate GoF when data are log transformed, for equivalence to PsychEncode analyses
  gof.parik <- list()
  for(j in deconv) {
    
    gof.parik[[j]] <- list()
    
    for (k in names(est.parik[[j]])) {
      
      print(paste0(Sys.time(),": ", j, "_", k))
      gof.parik[[j]][[k]] <- write.gof.v2(measuredExp = parikshak, 
                                          estimatedComp = est.parik[[j]][[k]], 
                                          signatureUsed = sigsBrain[[k]])
    }
  }  


## Save
  save(gof.parik, file = "Goodness of Fit.rda")  


## Plot distributions
  ## Setup and annotate data  
    plot.data <- melt(gof.parik)
    
    # rename columns
    colnames(plot.data) <- c("Metric", "value", "Signature", "Algorithm")
    
    # just correlation as the GoF metric
    plot.data <- plot.data[which(plot.data$Metric == "r"),] # focus on correlation
    
    # add region information
    plot.data$Region <- parikshak.meta$RegionID 
    plot.data$Region[grep("ba", plot.data$Region)] <- "Cortex"
    plot.data$Region[grep("vermis", plot.data$Region)] <- "Cerebellum"
    plot.data$Region <- factor(plot.data$Region, levels = c("Cortex", "Cerebellum"))
    
    # add a tag to identify paired deconvolution runs
    plot.data$id <- paste0(plot.data$Signature, "_", plot.data$Algorithm)
    
    # add colour information
    plot.data$Signature <- factor(plot.data$Signature, levels = signature.plotting.order)
    levels(plot.data$Signature) <- signature.annotation
     
 ## Plot!
    violins <- function(alg, region = "Cortex", y = c(0.3, 0.8)) {
      
      ggplot(plot.data[which(plot.data$Algorithm == alg & plot.data$Region == region),], aes(x = Signature, y = value, colour = Signature, fill = Signature)) +
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
  pdf(file = "Goodness of Fit.pdf", height = 3, width = 4)
  violins("CIB")
  dev.off()

  
  pdf(file = "Final Goodness of Fit (Supplementary, All Algs and Regions).pdf", height = 6, width = 8)
  ggplot(plot.data, aes(x = Signature, y = value, colour = Signature, fill = Signature)) +
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
    
  

## Plot ranks!
    ## Summary plot
      plot.ranks <- function(subset, title, text.colour = "black", subneurons = FALSE) {
        ranks <- lapply(gof.parik, function(x) { sapply(x, function(y) median(y$r[subset])) })
        
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


      pdf(file = "Final Goodness of Fit (Ranks).pdf", height = 3, width = 3.75)
      plot.ranks(subset = ctx, title = "Cortex")
      plot.ranks(subset = c(1:ncol(parikshak))[-ctx], title = "Cerebellum")
      dev.off()

  

  
################################################################################################################################ #
## Trends in composition ----  



## Boxplot proportions
  ## Part 1: All deconvolution algorithm combinations
  plot.data <- lapply(est.parik[deconv], function(x) {
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
    
  
## Correlations in cell-type proportion across runs!
  # function
  plot.heatmap <- function(celltype, h = 5.5, w = 8) {
    # collect composition estimates from the ctx for the desired celltype
    cors <- lapply(est.parik, function(x) {
      if(class(x) == "data.frame") {
        return(x[,grep(celltype, colnames(x))])
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
    
    # pdf(file = paste0("Correlation Heatmap V2 ", celltype,".pdf"), height = h, width = w)
    heatmap.2(cors, dendrogram = "row", trace = "none",  density.info = "none", 
          key.title = "NA", key = FALSE, keysize = 2.5, key.xlab = "Spearman", 
          main = NULL, col = colours, breaks = c(-2,-1,seq(0, 1, 0.1)), na.color = "red", na.rm = FALSE,
          cexRow = 0.8, cexCol = 0.8, lwid = c(2,10), lhei = c(2,20), mar = c(4,4))
    # dev.off()
  }
  
  
  plot.heatmap("Astrocytes")
  plot.heatmap("Neurons")
  plot.heatmap("Endothelia")
  plot.heatmap("Microglia")
  plot.heatmap("Oligodendrocytes")