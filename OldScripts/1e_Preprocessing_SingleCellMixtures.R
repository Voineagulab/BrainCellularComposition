## In this script, the scme is generated

################################################################################################################################ #
## Setup ----

## Generic
rm(list=ls())
options(stringsAsFactors = FALSE)

## Set directory
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/")

## Functions and libraries
  source("../Scripts/Fun_Preprocessing.R")
  
## Annotation files
  load(file = "Preprocessed/geneInfo.rda")
  load("Preprocessed/exonicLength.rda")
  

################################################################################################################################ #
## Via random sampling  ----    

## Preprocess data from Darmanis et al.
  # read in
  darmanis <- read.csv("Raw/Darmanis2015_scRNAseq.csv")
  rownames(darmanis) <- darmanis[,1]
  darmanis <- darmanis[,-1]
  
  # add ensID
  darmanis <- addENSID(darmanis)
  
  # annotate
  darmanis_meta <- read.table("Raw/Darmanis2015_Meta.txt", sep = "\t", header = 1)
  darmanis_meta <- darmanis_meta[-grep("hybrid", darmanis_meta$cell_type_s),] # "hybrid" refers to transcriptomes that likely contain >1 cell
  rownames(darmanis_meta) <- darmanis_meta$Sample_Name_s
  darmanis <- darmanis[,rownames(darmanis_meta)]
  darmanis_meta$cell_type_s <- as.factor(darmanis_meta$cell_type_s)
  levels(darmanis_meta$cell_type_s) <- c("Astrocytes", "Endothelia", "FoetalQuiescent", "FoetalRep", "Microglia", "Neurons", "Oligodendrocytes", "OPC")
  darmanis_meta$cell_type_s <- as.character(darmanis_meta$cell_type_s)
  # colnames(darmanis) <- as.character(darmanis_meta$cell_type_s[match(colnames(darmanis), darmanis_meta$Sample_Name_s)])
  # colnames(darmanis) <- sapply(strsplit(colnames(darmanis), "\\."), "[", 1)

  # filter to five key cell-types
  keep <- which(darmanis_meta$cell_type_s %in% c("Astrocytes", "Endothelia", "Microglia", "Neurons", "Oligodendrocytes"))
  darmanis <- darmanis[,keep]
  darmanis_meta <- darmanis_meta[keep,]
  
  
## Create signature!
  # dataframe
  sigsSCME <- as.data.frame(matrix(nrow = nrow(darmanis), ncol = 5))
  colnames(sigsSCME) <- levels(as.factor(darmanis_meta$cell_type_s))
  rownames(sigsSCME) <- rownames(darmanis)

  # aggregation
  for(j in colnames(sigsSCME)) { sigsSCME[,j] <- rowSums(darmanis[,which(darmanis_meta$cell_type_s == j)]) }
  
  # cpm
  sigsSCME <- as.data.frame(apply(sigsSCME, 2, function(x) {
    lib.size <- 10^6 / sum(x)
    x <- x * lib.size  
    return(x)
  }))

  # expression threshold based on signature
  remove <- which(apply(sigsSCME, 1, max) < 1)
  sigsSCME <- sigsSCME[-remove,]
  darmanis <- darmanis[-remove,]
  
    
## Mixture creation
  nMix <- 100
  nReps <- 100
  
  scme.v2 <- create.snmeRand(x = darmanis,
                              nMix = nMix,
                              nReps = nReps,
                              meta = darmanis_meta,
                              meta.columns = "cell_type_s",
                              rpkm = FALSE)

## Save
  # rda
  save(scme.v2, sigsSCME, file = "Preprocessed/scme.rda")

  # cibersort-compatible scmes
  write.CIB(data = scme.v2$mixture$cpm, dir = "Preprocessed/CIB/scme.v2.txt")
  
################################################################################################################################ #
## Via gradient sampling  ----    

## Create a gradient!
  colnames(darmanis) <- darmanis_meta$cell_type_s
  ct <- names(table(darmanis_meta$cell_type_s))
  
  ## Get abundances
    norm.abundance <- aggregate(colSums(darmanis), list(darmanis_meta$cell_type_s), mean)
    # norm.abundance$x <- max(norm.abundance$x) / norm.abundance$x
    norm.abundance$OverMax <- max(norm.abundance$x) / norm.abundance$x
    norm.abundance$OverMin <- norm.abundance$x / min(norm.abundance$x)
    
    base.n <- (nMix) / sum(norm.abundance$OverMax)
    norm.n <- floor(norm.abundance$OverMax * base.n)
    
    # if this number is greater than the total available, set it to the latter
    is.greater <- norm.n > table(colnames(darmanis)) 
    norm.n[is.greater] <- table(colnames(darmanis))[is.greater] 
    
    norm.n <- lapply(norm.n, function(y) { data.frame(c(0, y, 0, y)) })
    norm.n <- do.call("cbind", norm.n)
    colnames(norm.n) <- norm.abundance$Group.1
    props <- norm.n
    
  
  # run
  scme.gradients <- confound.proportion(dat = darmanis,
                                               ct = ct, 
                                               setProps = props, 
                                               nPerGroup = 50, 
                                               nPerMixture = nMix,
                                               keepPureExpression = FALSE)
  
  
  
  # rpkm normalise
  scme.gradients$confound.p <- rpkm(scme.gradients$confound.p)
  
  ## Save
  save(scme.gradients, file = "Preprocessed/scme.gradients.rda")  
  

############################################################# FIN ################################################################
