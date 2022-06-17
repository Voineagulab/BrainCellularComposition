## In this file, I process data for the RNA Mixture Experiment

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

############################################################################################################################### #
## This section reformats the RPKM data from three neuronal/astrocytic mixtures ----

## Load
load("/Volumes/Data1/PROJECTS/Neurons_Mixture_experiment/RESULTS/STAR_Output/RPKM.rda") # Using premade STAR output... 
  rme = rpkm[,c(2,5,6,7)] # These columns are my three mixtures!
  rme = rme[,-1]

## Filter to protein coding genes
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
  
  m = which(rownames(rme)%in%pcRNA.list)
  
  rme.pcRNA = rme[m,]
  
## Expression threshold
  # Min threshold: where the expression must be above 2rpkm in at least on reference
  rme.pcRNA.Tmin = rme.pcRNA[which(apply(rme.pcRNA, 1, max) > 2),]

  

## Saving references
  save(rme.pcRNA.Tmin, file = "Raw/rme_pcRNA_Tmin.rda")
  
############################################################################################################################### #
## This section creates a dataframe for the standardised neurons and astrocyte IH reference  ----


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


save(refIV.pcRNA, file = "Raw/rme_refIV_pcRNA.rda")


  
################################################################################################################################ #
## Create benchmarking mixture dataset: RNA from pure astrocyte and pure neuronal culture ----    
  
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
