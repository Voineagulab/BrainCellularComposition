## In this script, all signatures are preprocessed

## Note that the script Preprocessing_a_SingleNucleusData must be run first
## Specifically, this file is needed: Data/Preprocessed/Seurat.rda 

################################################################################################################################ #
## Setup ----

## Generic
rm(list=ls())
options(stringsAsFactors = FALSE)

## Directory
  # root.dir <- "INSERT.YOUR.ROOT.DIRECTORY.PLEASE/"
  root.dir <- "/Volumes/Data1/PROJECTS/BrainCellularComposition/BrainCellularComposition_GitHub/"
  wd1 <- paste0(root.dir, "Data/") # working directory 1, for processing
  setwd(wd1)


## Functions and libraries
source("../Scripts/Fun_Preprocessing.R")
library(tidyr)

## Parameters and settings
  # expression threshold
  exp_thresh <- 1
  
  # lists to hold the different flavours of signature.
  sigsMuSiC <- sigsBrain <- sigsSNME <- list()
  # sigsBrain: signatures containing up to five cell-types, for application to the scme and healthy brain transcriptomes
  # sigsSNME: signatures for use in the snme. Not loaded in the above line, but the previous load
  # sigsMuSiC: specially-formatted signatures required for use with the MuSiC algorith
  
  # annotation files
  load(file = "Preprocessed/geneInfo.rda")
  load("Preprocessed/exonicLength.rda")

################################################################################################################################ #
## F5 ----    

## Signature 1: A CAGE expression resource of (cultured) neuron and astrocyte from the FANTOM5 consortium

# read in data
F5.sig <- read.csv("Raw/Fantom5Consortium_CAGE_Filtered.csv") 
F5.meta <- read.csv("Raw/Fantom5Consortium_Meta.csv") 

# locate neurons and astrocytes
n <- grep("Neurons",F5.meta$Sample)
a <- grep("Astr",F5.meta$Sample)
nL <- F5.meta$Library[n[1:3]] 
nS <- F5.sig[,colnames(F5.sig)%in%nL] # Expression data for the cortical neuronal samples
aL <- F5.meta$Library[a[4:6]]
aS <- F5.sig[,colnames(F5.sig)%in%aL] # Expression data for the cortical astrocyte samples (i.e., ignoring cerebellar astrocytes)

# basic expression profiles of astrocytes and neurons  
ids <- F5.sig$X
F5.sig <- cbind(rowMeans(nS),rowMeans(aS)) # Keeping only the average of the three astrocyte or neuron samples
rownames(F5.sig) <- ids
colnames(F5.sig) <- c("Neurons", "Astrocytes")
F5.sig <- as.data.frame(F5.sig)


# adding EnsIDs to the annotation!
F5.sig <- addENSID(F5.sig)

# min threshold: where the expression must be above 1 tpm in at least on cell-type
F5.sig <- F5.sig[which(apply(F5.sig, 1, max) > exp_thresh),]

# add to the signature list
sigsBrain$F5 <- F5.sig[,c("Neurons", "Astrocytes")]
sigsSNME$F5 <- addSymbol(F5.sig[,c("Neurons", "Astrocytes")]) # though not CPM, it's the only available normalisation


################################################################################################################################ #
## IP ----    

## Signature 2: An RNA-seq expression resource of purified cell-types isolated from the human brain, without culture. 

# download table S4 from Zhang et al., 2016

## Read in 
  IP.sig <- readxl::read_xlsx("Raw/1-s2.0-S0896627315010193-mmc3.xlsx", sheet = 2, col_names = FALSE)
  IP.sig <- as.data.frame(IP.sig)

## Set cell-type annotation to columns
  # as we are reading a .xlsx file, the column name structure is non-standard. 
  # it annotates cell-type, and columns with NA are actually the annotation to its left.
  annot <- as.character(IP.sig[1,])
  annot[1] <- "Symbol"
  for (j in 1:length(annot)) {
    if (is.na(annot[j])) { # if NA
      closest <- which(!(is.na(annot[1:j]))) %>% max() # find the closest non NA to the left
      annot[j] <- annot[closest] # and assign it
    }
  }

  colnames(IP.sig) <- annot

## Other formatting
  # remove rows 1:3: annotation
  IP.sig <- IP.sig[-c(1:3),]
  
  # set rownames to be $Symbol (column 1)
  dup <- which(duplicated(IP.sig$Symbol))
  IP.sig <- IP.sig[-dup,]
  
  rownames(IP.sig) <- IP.sig$Symbol
  IP.sig <- IP.sig[,-1]

  # convert to EnsID
  IP.sig <- addENSID(IP.sig)
  
  # make expression numeric
  labs <- rownames(IP.sig)
  IP.sig <- apply(IP.sig, 2, as.numeric) %>% as.data.frame()
  rownames(IP.sig) <- labs

# aggregate samples within a cell-type
  n <- colnames(IP.sig)[grep("Neurons", colnames(IP.sig))]
  a <- colnames(IP.sig)[grep("mature astrocytes", colnames(IP.sig))]
  o <- colnames(IP.sig)[grep("Oligodendrocytes", colnames(IP.sig))]
  m <- colnames(IP.sig)[grep("Microglia", colnames(IP.sig))]
  e <- colnames(IP.sig)[grep("Endothelial", colnames(IP.sig))]

IP.sig <- data.frame(Neurons = IP.sig[,n],
                     Astrocytes = rowMeans(IP.sig[,a]),
                     Oligodendrocytes = rowMeans(IP.sig[,o]),
                     Microglia = rowMeans(IP.sig[,m]),
                     Endothelia = rowMeans(IP.sig[,e]))

# apply minimum threshold of 1 fpkm in any one cell-type
IP.sig <- IP.sig[which(apply(IP.sig, 1, max) > exp_thresh),]

# add to the signature list
sigsBrain$IP <- IP.sig
sigsSNME$IP <- addSymbol(IP.sig) # though not CPM, it's the only normalisation


################################################################################################################################ #
## DM ----    

## Signature 3: A single-cell RNAseq expression resource from surgically-removed temporal lobe tissue. Sourced from Darmanis et al 2015.
  # from https://github.com/VCCRI/CIDR-examples/tree/master/Brain, download: brainTags.csv and SraRunTable.txt. 
  # to their filenames I have appended _Darmanis2015 for labelling clarity
  
  # tag tables were downloaded from NCBI Gene Expression Omnibus (GSE67835). They were combined into one, and hybrid cells were excluded.
  DM.sig <- read.csv("Raw/brainTags_Darmanis2015.txt")
  rownames(DM.sig) <- DM.sig[,1]
  DM.sig <- DM.sig[,-1]
  
  # read in annotation
  DM.meta <- read.csv("Raw/SraRunTable_Darmanis2015.txt", sep = "\t")
  DM.meta <- DM.meta[-grep("hybrid", DM.meta$cell_type_s),]
  
  # add ensID
  DM.sig <- addENSID(DM.sig)
  DM.mus <- DM.sig
  DM.cpm <- apply(DM.sig, 2, function(x) {
      lib.size <- 10^6 / sum(x)
      x <- x * lib.size  
      return(x)  
  })
  DM.sig <- rpkm(DM.sig)
  
  # aggregate samples within a cell-type
  n <- colnames(DM.sig)[which(DM.meta$cell_type_s == "neurons")]
  a <- colnames(DM.sig)[which(DM.meta$cell_type_s == "astrocytes")]
  o <- colnames(DM.sig)[which(DM.meta$cell_type_s == "oligodendrocytes")]
  m <- colnames(DM.sig)[which(DM.meta$cell_type_s == "microglia")]
  e <- colnames(DM.sig)[which(DM.meta$cell_type_s == "endothelial")]
  p <- colnames(DM.sig)[which(DM.meta$cell_type_s == "OPC")]
  
  DM.sig <- data.frame(Neurons = rowMeans(DM.sig[,n]),
                       Astrocytes = rowMeans(DM.sig[,a]),
                       Oligodendrocytes = rowMeans(DM.sig[,o]),
                       Microglia = rowMeans(DM.sig[,m]),
                       Endothelia = rowMeans(DM.sig[,e]),
                       OPCs = rowMeans(DM.sig[,p]))
  
  DM.cpm <- data.frame(Neurons = rowMeans(DM.cpm[,n]),
                       Astrocytes = rowMeans(DM.cpm[,a]),
                       Oligodendrocytes = rowMeans(DM.cpm[,o]),
                       Microglia = rowMeans(DM.cpm[,m]),
                       Endothelia = rowMeans(DM.cpm[,e]),
                       OPCs = rowMeans(DM.cpm[,p]))
  
  # threshold
  DM.sig <- DM.sig[which(apply(DM.sig, 1, max) > exp_thresh),]
  DM.cpm <- DM.cpm[which(apply(DM.cpm, 1, max) > exp_thresh),]
  
  # add to signature list
  sigsBrain$DM <- DM.sig
  sigsSNME$DM <- addSymbol(DM.cpm)
  

## Create Music signature
  ## Process metadata
    DM.meta <- DM.meta[-grep("fetal", DM.meta$cell_type_s),]
    DM.meta$cell_type_s <- factor(DM.meta$cell_type_s)
    levels(DM.meta$cell_type_s) <- c("Astrocytes", "Endothelia", "Microglia", "Neurons", "Oligodendrocytes", "OPCs")
    DM.meta <- data.frame(orig.ident = DM.meta$Sample_Name_s,
                          Individual = DM.meta$age_s, # age is a proxy for individual
                          orig.celltype = DM.meta$cell_type_s,
                          row.names = DM.meta$Sample_Name_s)
  ## Process counts data
    DM.mus <- DM.mus[,rownames(DM.meta)]
    DM.mus <- DM.mus[rownames(DM.sig),] # expression filtering to those in the rpkm version
    
  ## Save
    sigsMuSiC$DM <- list(counts = DM.mus,
                         meta = DM.meta)

################################################################################################################################ #
## MM ----    

## Signature 4: similar derivation technique to sigIP, except from the mouse brain. MM stands for Mus Musculus

# read in 
# MM.sig <- read.csv("Raw/Zhang2014_mouseRNASeq.csv", row.names = 1)
MM.sig <- readxl::read_xlsx("Raw/Zhang2014_mouseRNASeq_originalDownload.xlsx", sheet = 1)
MM.sig <- as.data.frame(MM.sig)
rownames(MM.sig) <- MM.sig$`Gene symbol`


# filter to desired columns
MM.sig <- MM.sig[,c("Astrocytes", "Neuron", "Myelinating Oligodendrocytes", "Microglia", "Endothelial Cells")] # removes annotation columns

# filter to homologous genes
homologues <- read.delim("Raw/HOM_MouseHumanSequence.txt") # downloaded from: http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
MM.sig <- MM.sig[which(rownames(MM.sig) %in% homologues$Symbol),]

# replace mouse gene symbol with human gene symbol
m <- match(rownames(MM.sig), homologues$Symbol)
m <- m + 1 # to reflect the fact that the human homologue of a mouse gene in row m is found in row m + 1
MM.sig$Symbol <- homologues$Symbol[m]
dup <- MM.sig$Symbol[which(duplicated(MM.sig$Symbol))]
MM.sig <- MM.sig[-which(MM.sig$Symbol %in% dup),]

rownames(MM.sig) <- MM.sig$Symbol
MM.sig <- MM.sig[,-6]

# replace human gene symbol with human EnsID
MM.sig <- addENSID(MM.sig)

# rename columns
colnames(MM.sig) <- c("Astrocytes", "Neurons", "Oligodendrocytes", "Microglia", "Endothelia")
MM.sig <- MM.sig[,c("Neurons", "Astrocytes","Oligodendrocytes", "Microglia", "Endothelia")]

# threshold
MM.sig <- MM.sig[which(apply(MM.sig, 1, max) > exp_thresh),]

# add to signature list
sigsBrain$MM <- MM.sig
sigsSNME$MM <- addSymbol(MM.sig) # again, not properly normalised, but only available option...

# ################################################################################################################################ #
# ## IH ----    
# 
# ## Signature 5: in-house generated, it contains RNA from pure astrocyte and pure neuron cell culture
# load("Raw/rme_refIV_pcRNA.rda")
# sigsRME$IH <- refIV.pcRNA; rm(refIV.pcRNA)
# colnames(sigsRME$IH) <- c("Neurons", "Astrocytes")  
# sigsRME$IH <- sigsRME$IH[which(apply(sigsRME$IH, 1, max) > exp_thresh),]
# # sigsBrain$IH <- sigsRME$IH[,c("Neurons", "Astrocytes")]


################################################################################################################################ #
## Single-nucleus signatures from processed Seurat data: VL, NG, CA, TS, and LK ----


## First, for sigsSNME, I note that CA and VL will have a signature separately generated in the snme script, for use only with the snme rather than other datasets


## Load
  load("Preprocessed/SeuratObjects.rda")


## The first step is to reannotate cell-types
  # for the brain, we will be using Neurons, Astrocytes, Microglia, Endothelia, Oligodendrocytes, and OPCs. 
  # excitatory and inhibitory subtypes of neurons will be present in a separate signature...

## Function
rename.ct <- function(x, orig.label, new.label) {
  # create metadata column if not already present
  if (!("brain.ct" %in% colnames(x@meta.data))) x$brain.ct <- "."
  
  # get indices of cells with the original labl
  g <- grep(orig.label, x$orig.celltype)
  
  # rename these
  x$brain.ct[g] <- new.label
  
  return(x)
}

# VL
obj$VL <- rename.ct(obj$VL, "AST", "Astrocytes")
obj$VL <- rename.ct(obj$VL, "^L|^Neu", "Excitatory")
obj$VL <- rename.ct(obj$VL, "^IN", "Inhibitory")
obj$VL <- rename.ct(obj$VL, "Endo", "Endothelia")
obj$VL <- rename.ct(obj$VL, "Microglia", "Microglia")
obj$VL <- rename.ct(obj$VL, "Oligodendrocytes", "Oligodendrocytes")
obj$VL <- rename.ct(obj$VL, "OPC", "OPCs")

# NG
obj$NG <- rename.ct(obj$NG, "Astro", "Astrocytes")
obj$NG <- rename.ct(obj$NG, "Inhib", "Inhibitory")
obj$NG <- rename.ct(obj$NG, "Ex", "Excitatory")
obj$NG <- rename.ct(obj$NG, "Endo", "Endothelia")
obj$NG <- rename.ct(obj$NG, "Micro", "Microglia")
obj$NG <- rename.ct(obj$NG, "Oligo", "Oligodendrocytes")
obj$NG <- rename.ct(obj$NG, "OPC", "OPCs")

# CA
obj$CA <- rename.ct(obj$CA, "Astro", "Astrocytes")
obj$CA <- rename.ct(obj$CA, "Inh", "Inhibitory")
obj$CA <- rename.ct(obj$CA, "Exc", "Excitatory")
obj$CA <- rename.ct(obj$CA, "Endo", "Endothelia")
obj$CA <- rename.ct(obj$CA, "Micro", "Microglia")
obj$CA <- rename.ct(obj$CA, "Oligo", "Oligodendrocytes")
obj$CA <- rename.ct(obj$CA, "OPC", "OPCs")

# TS
obj$TS <- rename.ct(obj$TS, "Astro", "Astrocytes")
obj$TS$brain.ct[which(obj$TS$Class == "GABAergic")] <- "Inhibitory" # a simpler way of getting ct information
obj$TS$brain.ct[which(obj$TS$Class == "Glutamatergic")] <- "Excitatory" # per above
obj$TS <- rename.ct(obj$TS, "Endo", "Endothelia")
obj$TS <- rename.ct(obj$TS, "Micro|^PVM", "Microglia")
obj$TS <- rename.ct(obj$TS, "Oligo", "Oligodendrocytes")
obj$TS <- rename.ct(obj$TS, "OPC", "OPCs")

# LK
obj$LK <- rename.ct(obj$LK, "Ast", "Astrocytes")
obj$LK <- rename.ct(obj$LK, "In", "Inhibitory")
obj$LK <- rename.ct(obj$LK, "Ex", "Excitatory")
obj$LK <- rename.ct(obj$LK, "End", "Endothelia")
obj$LK <- rename.ct(obj$LK, "Mic", "Microglia")
obj$LK <- rename.ct(obj$LK, "Oli", "Oligodendrocytes")
obj$LK <- rename.ct(obj$LK, "OPC", "OPCs")

# counts
ct.counts <- sapply(obj, function(x) table(x$brain.ct))
ct.counts <- do.call("cbind", ct.counts)


## Next, create CPM and RPKM signatures
  ## Function
    create.seurat.signature <- function(w) {
      # print
      print(w@project.name)
      
      # get counts
      counts <- as.data.frame(w@assays$RNA@counts)
      
      # convert to EnsID and remove non-coding genes
      counts <- addENSID(counts)
      
      # get CPM of every cell
      # cpm <- apply(counts, 2, function(x) {
      #   lib.size <- 10^6 / sum(x)
      #   x <- x * lib.size
      #   return(x)
      # })
      # below is an alternate way, easier on RAM 
      cpm <- counts
      for (j in 1:ncol(cpm)) {
          if (j %% 1000 == 0) print(paste0("Processing sample ", j))
          lib.size <- 10^6 / sum(cpm[,j])
          cpm[,j] <- cpm[,j] * lib.size
      }
      
      # cpm <- as.data.frame(cpm)
      
      # get RPKM of every cell
      rpkm <- length.correct(cpm)
      
      # signatures are the average normalised expression of every member
      output <- list(rpkm = list(), cpm = list())
      for (j in rownames(ct.counts)) {
        # print((j))
        k <- which(w$brain.ct == j)
        
        output$cpm[[j]] <- rowMeans(cpm[,k]) 
        output$rpkm[[j]] <- rowMeans(rpkm[,k]) 
      }
      
      output <- lapply(output, function(x) as.data.frame(do.call("cbind", x)))
      
      # add neurons, which come from pooling exc and inh cells
      neu <- which(w$brain.ct %in% c("Excitatory", "Inhibitory"))
      
      output$cpm$Neurons <- rowMeans(cpm[,neu])
      output$rpkm$Neurons <- rowMeans(rpkm[,neu])
      
      
      # expression threshold: a gene is kept if > 1 unit in at least 1 cell-type
      output <- lapply(output, function(x) {
        keep <- which(apply(x, 1, max, na.rm = TRUE) > 1)
        x <- x[keep,]
        return(x)
      })
      
      # return
      return(output)
    }

## Apply function
  x <- lapply(obj, create.seurat.signature)

## Partition output into rpkm and cpm lists
  # rpkm, for use with scme / bulk brain data
  sigsBrain <- c(sigsBrain, lapply(x, function(y) y$rpkm))
  
  # cpm, for use with the snmes
  sigsSNME$TS <- addSymbol(x$TS$cpm) # sigsSNME already has specially-prepared CA and VL signatures, so add the others
  sigsSNME$LK <- addSymbol(x$LK$cpm)
  sigsSNME$NG <- addSymbol(x$NG$cpm)



## Finally, create MuSiC compatible signatures
  ## Function
    create.seurat.MuSiC <- function(w) {
      # print
      print(w@project.name)
      
      # get counts
      counts <- as.data.frame(w@assays$RNA@counts)
      
      # convert to EnsID and remove non-coding genes
      counts <- addENSID(counts)
      
      # grab metadata
      meta <- w@meta.data
      
      # add another metadata column where exc and inh are combined into neurons
      meta$brain.ct2 <- meta$brain.ct
      meta$brain.ct2[grep("Exc|Inh", meta$brain.ct2)] <- "Neurons"
      
      # return
      output <- list(counts = counts, 
                     meta = meta)
      return(output)  
    }

  ## Apply function
    sigsMuSiC <- lapply(obj, create.seurat.MuSiC)
    
  ## Expression threshold the signatures based on the corresponding brain signature
    for (j in names(sigsMuSiC)) {
      sigsMuSiC[[j]]$counts <- sigsMuSiC[[j]]$counts[rownames(sigsBrain[[j]]),]
    }
    

## Remove ct
  sigsBrain$CA <- sigsBrain$CA[,-which(colnames(sigsBrain$CA) == "Endothelia")]
  sigsBrain$LK <- sigsBrain$LK[,-which(colnames(sigsBrain$LK) == "Endothelia")]
    
## Save
save(sigsBrain, file = "Preprocessed/Signatures - Brain.rda")
save(sigsMuSiC, file = "Preprocessed/Signatures - Brain (MuSiC).rda")
save(sigsSNME, file = "Preprocessed/Signatures - SNME (Incomplete).rda") # marked as incomplete as it is missing CA and VL; these are generated when the corresponding mixtures are generated


################################################################################################################################ #
## Multibrain ----  
  
## Five human signatures, without subtype resolution
  x <- sigsBrain[c("IP", "DM", "VL", "NG", "CA")]
  
  # filter to common ct (columns) and genes (rows)
  common <- Reduce(intersect, lapply(x, rownames)) # this gets the genes (rownames) common to all signatures
  x <- lapply(x, function(y) y[common, colnames(x$IP)])
  
  # quantile normalise
  qn <- do.call("cbind", x)
  # qn <- qn[,-which(colnames(qn) == "CA.Endothelia")] # note that this is a difference to the paper, that version had CA's endothelia. it is removed here for convenience in these replication scripts
  qn <- as.data.frame(normalize.quantiles(as.matrix(qn), copy = FALSE))
  
  # average across signatures
  multibrain <- list()
  for (j in colnames(x$IP)) {
    multibrain[[j]] <- rowMeans(qn[,grep(j, colnames(qn))])
  }
  multibrain <- as.data.frame(do.call("cbind", multibrain))

## Save
  save(multibrain, file = "Preprocessed/Multibrain.rda")


############################################################# FIN ################################################################
