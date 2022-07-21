################################################################################################################################ #
## Setup ----     


## Start!
  rm(list = ls()); gc()

## Directory
  # root.dir <- "INSERT.YOUR.ROOT.DIRECTORY.PLEASE/"
  root.dir <- "/Volumes/Data1/PROJECTS/BrainCellularComposition/BrainCellularComposition_GitHub/"
  wd1 <- paste0(root.dir, "Data/") # working directory 1, for processing
  wd2 <- paste0(root.dir, "Results/BulkTissue/GTEx/") # and for analysis
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
  
  
## Read in count-level data
  # the downloaded file is automatically named: GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz
  # please unzip it first!
  gtex <- read.delim("Raw/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct", skip = 2) 
  
  # the above command may take many minutes, and requires ~3GB of memory

  # move gene annotation to rownames
  rownames(gtex) <- sapply(strsplit(as.character(gtex$Name), split = "\\."), `[`, 1)
  gtex <- gtex[,-c(1,2)]

  # colnames
  colnames(gtex) <- gsub("\\.", "-", colnames(gtex))
  
## Read in metadata
  # gtexMeta <- read.csv("Raw/GTEx_v7_Annotations_SampleAttributesDS.csv")
  gtexMeta <- read.delim("Raw/GTEx_v7_Annotations_SampleAttributesDS.txt")
  gtexMeta <- gtexMeta[which(gtexMeta$SAMPID %in% colnames(gtex)),]

## Exclude non-brain samples
  m <- match(colnames(gtex), gtexMeta$SAMPID)
  gtexMeta <- gtexMeta[m,]
  brains <- which(gtexMeta$SMTS == "Brain")
  gtex <- gtex[,brains]
  gtexMeta <- gtexMeta[brains,]

## Convert to numeric
  for (j in 1:ncol(gtex)) gtex[,j] <- as.numeric(as.character(gtex[,j]))

## RPKM
  gtexCounts <- gtex
  gtex <- rpkm(gtex)

## Threshold
  gtex <- gtex[rowSums(gtex >= 1) >= 88,] # 88 is the size of the smallest group
  gtexCounts <- gtexCounts[rownames(gtex),]

## Furnish the sequencing metadata with patient information
  # match
  patientInfo <- read.delim("Raw/GTEx_v7_Annotations_SubjectPhenotypesDS.txt")
  IDs <- strsplit(as.character(gtexMeta$SAMPID), split = "-")
  IDs <- paste0(sapply(IDs, `[`, 1), "-", sapply(IDs, `[`, 2))
  m <- match(IDs, patientInfo$SUBJID)

  # add age
  gtexMeta$Age <- patientInfo$AGE[m]

  # add gender
  gtexMeta$Sex <- patientInfo$SEX[m]
  
  # add broad region
  gtexMeta$BroadRegion <- "sCTX" # sub-cortical
  gtexMeta$BroadRegion[grep("Cerebell", gtexMeta$SMTSD, ignore.case = TRUE)] <- "CB" # cerebellar
  gtexMeta$BroadRegion[grep("Cortex", gtexMeta$SMTSD, ignore.case = TRUE)] <- "CTX" # cortex
  gtexMeta$BroadRegion[which(gtexMeta$SMTSD == "Brain - Spinal cord (cervical c-1)")] <- "SP" # spinal cord
  
## Add gene symbol
  gtexSymbols <- addSymbol(gtex)

  
## Save
  save(gtex, gtexMeta, gtexCounts, gtexSymbols, file = "Preprocessed/GTEx.rda")
  write.CIB(data = gtex, dir = "Preprocessed/CIB/GTEx.txt")

  
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
  gtexMeta$BroadRegion <- "Subcortex" # sub-cortical
  gtexMeta$BroadRegion[grep("Cerebell", gtexMeta$SMTSD, ignore.case = TRUE)] <- "Cerebellum" # cerebellar
  gtexMeta$BroadRegion[grep("Cortex", gtexMeta$SMTSD, ignore.case = TRUE)] <- "Cortex" # cortex
  gtexMeta$BroadRegion[which(gtexMeta$SMTSD == "Brain - Spinal cord (cervical c-1)")] <- "Spinal Cord" # spinal cord
  
  ctx <- which(gtexMeta$BroadRegion == "Cortex")  deconv <- c("DRS", "DTA", "CIB", "MUS")
  enrich <- c("Blender", "Coex", "xCell")
  main.sigs <- c("CA", "NG", "LK", "VL", "TS", "MM", "DM", "IP", "F5", "MB")    
  
  
################################################################################################################################ #
## Estimate composition ----  
  
## Setup
est.gtex <- list(DRS = list(),
            DTA = list(),
            CIB = list(),
            MUS = list(),
            xCell = list(),
            Blender = list(),
            Coex = list(),
            Linseed = list())


## Apply to all samples
  ## DeconRNASeq
    
  # for (j in names(sigsBrain)[grep("z", names(sigsBrain))]) {
    for (j in names(sigsBrain)) {
      est.gtex$DRS[[j]] <- run.DRS(gtex, sigsBrain[[j]])
    }
  
  ## DTA
    for (j in names(sigsBrain)) {
    # for (j in names(sigsBrain)[grep("z", names(sigsBrain))]) {
      est.gtex$DTA[[j]] <- run.DTA(gtex, sigsBrain[[j]], alg = "diff", q = 0.01)
    }  
  
  
  ## CIBERSORT
    for(j in names(sigsBrain)) {
      # report, as runtime is ~hours
      run.no <- which(names(sigsBrain) == j)
      print(paste0(j, ", run ", run.no, " of ", length(sigsBrain), ": ", Sys.time()))
      
      est.gtex$CIB[[j]] <- run.CIB(from.file = FALSE,
                             mixString = "GTEx.txt",
                             sigObject =  sigsBrain[[j]])
      
      gc()
    
    }
  
  
  ## MUSIC
  for(j in names(sigsMuSiC)) {
      print(j)
      Sys.time()
      
        # celltypes to remove from the standard deconvolution
        remove <- which(sigsMuSiC[[j]]$meta$brain.ct2 == "OPCs") # remove OPCs
        
        if (j %in% c("CA", "LK")) {
          remove <- c(remove, grep("End", sigsMuSiC[[j]]$meta$brain.ct)) # remove endothelia from CA and LK signatures
        }
        
        # deconvolve with the core.ct
        gc()
        est.gtex$MUS[[j]] <- run.music(mixture = gtexCounts,
                                  use.meta.column = "brain.ct2",
                                  signature = sigsMuSiC[[j]]$counts[,-remove],
                                  signature.meta = sigsMuSiC[[j]]$meta[-remove,])
        
      # end loop
    }
    
  ## xCell
    est.gtex$xCell <- run.xCell(gtexSymbols)
  
  ## Blender
    est.gtex$Blender <- run.Blender(gtexSymbols)
  
  
  ## Co-expression
    coex <- list()
    coex$Full <- run.coex(mixture = gtex, 
                          signature = sigsBrain$SC.CA.IP, 
                          only.threshold = FALSE, 
                          sft = 18, # these data  do not reach scale-free topology by a reasonable power; set at 18, where r2 = 0.692, and median.k = 77
                          output.all = TRUE)
    y <- define.celltypes(x = coex$Full, algorithm.type = "Coex", celltype.signature = sigsBrain$SC.CA.IP, return.confidence = TRUE, nMarks.x = 10000, nMarks.sig = 100, bg = rownames(gtex))
    colnames(y$comp) <- y$assignments
    est.gtex$Coex <- y$comp

    
################################################################################################################################ #
## Estimate composition in cortical samples ----  
  
## Create a separate list for CTX estimates
  # collect all variance-independent methods!
  est.gtex.ctx <- lapply(est.gtex[1:4], function(x) {
    lapply(x, function(y) {
      y[ctx,]
    })
  })
  
  # run blender
  est.gtex.ctx$Blender <- run.Blender(gtexSymbols[,ctx])
  
  # run xCell (though it's not always affected by variance)
  est.gtex.ctx$xCell <- run.xCell(gtexSymbols[,ctx])

  # run Linseed on cortical samples
  pdf(file = "Linseed Plots (CTX).pdf", height = 4, width = 4)
  est.gtex.ctx$Linseed <- run.linseed(mixture = gtex[,ctx], 
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
  est.gtex.ctx$Coex <- y$comp
  
## Save
save(est.gtex, est.gtex.ctx, file = "Composition Estimates (Final) (New MUS).rda") # load("Composition Estimates (Final) (New MUS).rda")





################################################################################################################################ #
## Evaluate using Goodness of Fit ----  

## Setup 
gof <- list()

## Calculate GoF when data are log transformed, for equivalence to PsychEncode analyses
for(j in deconv) {
  
  gof[[j]] <- list()
  
  for (k in names(est.gtex[[j]])) {
    
    if (k == "DM_") next # skip this run
    if (substr(k, 4, 4) == ".") next # skip deconvolution using subneurons + stable
    if (substr(k, 3, 3) == "z") next # skip deconvolution using all ct
    
    l <- gsub("\\.", "", k) # this means that if the stable deconvolution is used, GoF is calculated using a non-stable signature,
    
    print(paste0(Sys.time(),": ", j, "_", k))
    gof[[j]][[k]] <- write.gof.v2(measuredExp = gtex, 
                               estimatedComp = est.gtex[[j]][[k]], 
                               signatureUsed = sigsBrain[[l]])
  }
}  


## Save
save(gof, file = "Goodness of Fit (Final) (V2 with l).rda")  # load("Goodness of Fit (Final) (V2 with l).rda")


## Setup plotting data
  ## Setup and annotate data  
  plot.data <- melt(gof)
  
  # rename columns
  colnames(plot.data) <- c("Metric", "value", "Signature", "Algorithm")
  
  # just correlation as the GoF metric
  plot.data <- plot.data[which(plot.data$Metric == "r"),] # focus on correlation
  
  # add region information
  plot.data$Region <- gtexMeta$BroadRegion 
  plot.data$Region <- factor(plot.data$Region, levels = c("Cortex", "Cerebellum", "Subcortex", "Spinal Cord"))
  
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
  

  
## Plot! this is for standard (no filter) deconvolutions
  ## Function
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
  
  pdf(file = "Final Goodness of Fit (Supplementary, All Algs and Regions).pdf", height = 10.5, width = 8)
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

  pdf(file = "Final Goodness of Fit (Alternate Deconvolution Strategies).pdf", height = 10, width = 8)
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

  
  
  pdf(file = "Final Goodness of Fit (Ranks).pdf", height = 3, width = 3.75)
  plot.ranks(subset = ctx, title = "Cortex")
  plot.ranks(subset = which(gtexMeta$BroadRegion == "Cerebellum"), title = "Cerebellum")
  plot.ranks(subset = which(gtexMeta$BroadRegion == "Subcortex"), title = "Subcortex")
  plot.ranks(subset = which(gtexMeta$BroadRegion == "Spinal Cord"), title = "Spinal Cord")
  dev.off()
  
  
  

  

################################################################################################################################ #
## Trends in composition ----  


## Boxplot proportions
  ## Part 1: All deconvolution algorithm combinations
  plot.data <- lapply(est[deconv], function(x) {
    x <- x[main.sigs]
    lapply(x, function(y) { 
      y$Region <- gtexMeta$BroadRegion
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
  ## Function
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
  