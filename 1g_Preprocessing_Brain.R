## In this script,  bulk-brain-derived datasets from Parikshak et al. and GTEx are preprocessed

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
  
## Parameters
  exp_thresh <- 1 # minimum rpkm cutoff



################################################################################################################################ #
## Preprocessing of Parikshak et al. transcriptomes for broader analyses of the brain ----

## Data from Parikshak et al, 2016
  load("Raw/Parikshak2016_RNASeq.Rdata")
  parikshak <- datExpr.HTSC.unionexon
  parikshak.meta <- datMeta
  rm(datExpr.Cufflinks, datExpr.HTSC.unionexon, datExpr.HTSC.wholegene, datLB.gene, datMeta, datOK.gene)
  
## RPKM normalisation
  parikshak.counts <- parikshak
  parikshak <- rpkm(parikshak)
  
## Expression thresholding
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
  
## Save
  # rda
  save(parikshak, parikshak.symbol, parikshak.counts, parikshak.meta, file = "Preprocessed/Parikshak.rda")
  
  # cibersort-compatible files
  write.CIB(data = parikshak, dir = "Preprocessed/CIB/Parikshak.txt") 
  
  
## List of samples used in analyses
  samplesUsed <- data.frame(SampleID = rownames(parikshak.meta),
                            Region = parikshak.meta$RegionID,
                            Age = parikshak.meta$Age,
                            Diagnosis = parikshak.meta$ASD.CTL,
                            composition.estimated = parikshak.meta$Network.analysis..CTX | parikshak.meta$Network.Analysis..CB,
                            used.in.DE = parikshak.meta$ASD.vs.CTL..CTX,
                            used.in.WGCNA = parikshak.meta$Network.analysis..CTX)
  
  write.csv(samplesUsed, file = "Preprocessed/Samples Used In Parikshak.csv")
  

  
################################################################################################################################ #
## GTEx ----       
  
## Read in count-level data
  gtex <- read.table("Raw/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.txt")

  # adjust column names
  colnames(gtex) <- unlist(gtex[1,])
  gtex <- gtex[-1,]

  # move gene annotation to rownames
  rownames(gtex) <- sapply(strsplit(as.character(gtex$Name), split = "\\."), `[`, 1)
  gtex <- gtex[,-c(1,2)]

## Read in metadata
  gtexMeta <- read.csv("Raw/GTEx_v7_Annotations_SampleAttributesDS.csv")
  gtexMeta <- gtexMeta[which(gtexMeta$SAMPID%in%colnames(gtex)),]

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
  symbols <- addSymbol(gtex)

  
## Save
  save(gtex, gtexMeta, gtexCounts, symbols, file = "Preprocessed/GTEx.rda")
  write.CIB(data = gtex, dir = "Preprocessed/CIB/GTEx.txt")


## Samples used
  samplesUsed <- data.frame(SampleID = (gtexMeta$SAMPID),
                            Region = gtexMeta$SMTSD,
                            BroadRegion = gtexMeta$BroadRegion,
                            Age = gtexMeta$Age)
  
  write.csv(samplesUsed, file = "Preprocessed/Samples Used In GTEx.csv")
  write.csv(samplesUsed, file = "../Results/SuppTables/ST1_SampleDescriptions_SourceFiles/Samples Used In GTEx.csv")


################################################################################################################################ #
## Preprocessing of Parikshak et al. transcriptomes for analyses of ASD ----

## Data from Parikshak et al, 2016
  load("Raw/Parikshak2016_RNASeq.Rdata")
  pCounts <- datExpr.HTSC.unionexon[,-grep("vermis", colnames(datExpr.HTSC.unionexon))] # removing samples from the cerebellar vermins
  pMeta <- datMeta[-grep("vermis", rownames(datMeta)),] # removing samples from the cerebellar vermis
  rm(datExpr.Cufflinks, datExpr.HTSC.unionexon, datExpr.HTSC.wholegene, datLB.gene, datMeta, datOK.gene)
  
## RPKM normalisation
  pRPKM <- rpkm(pCounts)
  
## Expression thresholding
  NAs <- which(apply(pRPKM, 1, anyNA))
  if (length(NAs) > 0) pRPKM <- pRPKM[-which(apply(pRPKM, 1, anyNA)),] # first, remove any rows with an NA
  n <- round(min(table(pMeta$ASD.CTL)) / 2) # a gene must be expressed above exp_thresh in 41 samples (half the size of the smallest group)
  pRPKM <- pRPKM[rowSums(pRPKM >= exp_thresh) >= n , ] # here, pRPKM >= exp_thresh is a logical matrix, which we need to be TRUE in >= 41 sample
  
## Outlier removal: for standardisation, remove the same outliers as Parikshak et al in their original analyses (of which there are 4)
  outliers <- which(pMeta$Network.analysis..CTX == FALSE)
  pRPKM <- pRPKM[,-outliers]
  pMeta <- pMeta[-outliers,]
    
## Apply same standards to the count data
  pCounts <- pCounts[which(rownames(pCounts)%in%rownames(pRPKM)),
                     which(colnames(pCounts)%in%colnames(pRPKM))] 
    
## Features list
  pFeatures <- cbind(rownames(pRPKM), exonicLength[match(rownames(pRPKM), rownames(exonicLength))])
  pFeatures <- as.data.frame(pFeatures)
  colnames(pFeatures) <- c("EnsID", "ExonicLength")
  pFeatures$EnsID <- as.character(pFeatures$EnsID)
  pFeatures$ExonicLength <- as.numeric(pFeatures$ExonicLength)
  
  # add total gene length
  geneLength <- read.csv("/Volumes/Data1/PROJECTS/Nicole_TTseq/ANNOTATION/Ensembl_gene_length.csv")
  geneLength["Total_length"] <- NA 
  geneLength$Total_length <- geneLength$Gene.end..bp. - geneLength$Gene.start..bp.
  m <- match(pFeatures$EnsID, geneLength$EnsID)
  pFeatures$Total_length <- geneLength$Gene.end..bp.[m] - geneLength$Gene.start..bp.[m]
    
  pFeatures$Gene.Symbol <- pFeatures$Biotype <- "-" # placeholder
  pFeatures$Gene.Symbol <- geneInfo$Gene.Symbol[match(pFeatures$EnsID, geneInfo$ensID)]
  pFeatures$Biotype <- geneInfo$Biotype[match(pFeatures$EnsID, geneInfo$ensID)]
  
## Augment metadata. 
  # calculate summaries of the sequencing statistics per Parikshak's analyses. Code taken from dhglab's GitHub
  seqInfo <- data.matrix(pMeta[,c(25:43)])
  seqInfo <- seqInfo[,c(6:19,c(5))]
  seqInfo <- t(scale(seqInfo,scale=F))
  PC.seqInfo <- prcomp(seqInfo);
  varexp <- (PC.seqInfo$sdev)^2 / sum(PC.seqInfo$sdev^2)
  topPC.seqInfo <- PC.seqInfo$rotation[,1:2];
  colnames(topPC.seqInfo) <- c("SeqSV1","SeqSV2") # recompute since those in datMeta were with additional samples included
  pMeta$seqStatPC1 <- as.numeric(topPC.seqInfo[,1])
  pMeta$seqStatPC2 <- as.numeric(topPC.seqInfo[,2])
  
## Generate a version of pRPKM annotated with geneSymbol
  pSymbol <- addSymbol(pRPKM)
  
## Filter to the subset analysed in Parikshak et al.
  # get filter
  subset <- which(pMeta$ASD.vs.CTL..CTX) # these 106 samples are those which Parikshak et al used in their DE analyses; it excludes young and/or Dup15q samples

  # apply filter
  pRPKM <- pRPKM[,subset]
  pSymbol <- pSymbol[,subset]
  pMeta <- pMeta[subset,]
  pCounts <- pCounts[,subset]
  
## Save
  # rda
  save(pCounts, pRPKM, pSymbol, pMeta, pFeatures, file = "Preprocessed/ASD.rda")
  
  # cibersort-compatible files
  write.CIB(data = pRPKM, dir = "Preprocessed/CIB/pRPKM.txt") 

  
############################################################# FIN ################################################################