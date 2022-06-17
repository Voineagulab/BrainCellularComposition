##################################################################################################################################
## Libraries

require(data.table)
require(GenomicFeatures)
require(gplots)
require(ggplot2)
require(reshape2)
require(gtools)
require(dtangle)
require(preprocessCore)
require(Seurat)

## And load other files...
load("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/exonicLength.rda") # exonic lengths of genes, used for gene length normalisation
load("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/geneInfo.rda") # information about gene symbols, ids, and type

##################################################################################################################################
## General data preprocessing functions

## Convert counts to RPKM. There is one note in the function.
rpkm <- function(counts) { 
  counts <- counts[which(rownames(counts) %in% rownames(exonicLength)),]  # note 1: counts must have EnsID as rownames   
  
  # Format the exonicLength matrix/list
  length <- transpose(exonicLength)[[1]] # note 2: the file "exonicLength" must be loaded into the global environment,
  names(length) <- rownames(exonicLength)
  
  # Stats
  m <- match(rownames(counts), names(length))
  length <- length[m]/1000
  libsize <- apply(counts, 2, sum) / 10^6
  
  # Normalise for library size, then length
  cpm <- counts 
  for (j in c(1:ncol(cpm))) cpm[,j] <- counts[,j]/libsize[j]
    
  rpkm <- cpm 
  for (j in c(1:ncol(rpkm))) rpkm[,j] <- cpm[,j]/length
  
  return(rpkm)
}

## Convert counts to TPM. There is one note in the function.
tpm <- function(counts) { 
  counts <- counts[which(rownames(counts) %in% rownames(exonicLength)),]  # note 1: counts must have EnsID as rownames   
  
  # Format the exonicLength matrix/list
  length <- transpose(exonicLength)[[1]] # note 2: the file "exonicLength" must be loaded into the global environment,
  names(length) <- rownames(exonicLength)
  
  # Stats on length
  m <- match(rownames(counts), names(length))
  length <- length[m]/1000
  
  # Normalise for length
  length.corrected <- counts
  for (j in c(1:ncol(length.corrected))) length.corrected[,j] <- counts[,j]/length
  
  # Stats on libsize
  libsize <- apply(length.corrected, 2, sum) / 10^6
  
  # Normalise for libsize
  tpm <- length.corrected
  for (j in c(1:ncol(tpm))) tpm[,j] <- length.corrected[,j]/libsize[j]
  
  return(tpm)
}

## Divide expression in a dataframe by length in kilobases
  length.correct <- function(exp) {
    exp <- exp[which(rownames(exp) %in% rownames(exonicLength)),]  
    
    length <- transpose(exonicLength)[[1]] # note 2: the file "exonicLength" must be loaded into the global environment,
    names(length) <- rownames(exonicLength)
    
    m <- match(rownames(exp), names(length))
    length <- length[m]/1000
    
    norm <- apply(exp, 2, function(x) x / length)  
    norm <- as.data.frame(norm)
    return(norm)
  }


## Convert rownames from EnsID to GeneSymbol. 1 note.
addSymbol <- function(dataframe) {
  annot <- geneInfo[which(geneInfo$Biotype == "protein_coding"),] # Note 1: selects protein-coding genes only! 
  
  sharedAnnotation <- annot[which(annot$ensID%in%rownames(dataframe)),]
  matches <- match(sharedAnnotation$ensID, rownames(dataframe))
  
  # Add gene.symbol to a new column
  dataframe$Gene.Symbol <- "-"
  dataframe$Gene.Symbol[matches] <- sharedAnnotation$Gene.Symbol 
  dataframe <- dataframe[which(dataframe$Gene.Symbol != "-"),]      
  dataframe <- dataframe[!(duplicated(dataframe$Gene.Symbol) | duplicated(dataframe$Gene.Symbol, fromLast = TRUE)), ] # remove deprecated entries
  
  # put gene.symbol into rownames
  rownames(dataframe) <- dataframe$Gene.Symbol
  dataframe <- dataframe[,-which(colnames(dataframe) %in% "Gene.Symbol")]
  return(dataframe)
}  
  

## Convert rownames from GeneSymbol to EnsID. 1 note.
addENSID <- function(dataframe, pc = TRUE) {
  if (pc) {
    annot <- geneInfo[which(geneInfo$Biotype == "protein_coding"),] # Note 1: selects protein-coding genes only!  
  } else {
    annot <- geneInfo
  }
  
  # Match symbols in dataframe and annotation file
  sharedAnnotation <- annot[which(annot$Gene.Symbol%in%rownames(dataframe)),]
  matches <- match(sharedAnnotation$Gene.Symbol, rownames(dataframe))
  
  # Add ensID to a new column
  dataframe$ensID <- "-"
  dataframe$ensID[matches] <- sharedAnnotation$ensID # A sanity check was performed, which bound m$Approved.Symbol, and it always matched to rownames!
  dataframe <- dataframe[which(dataframe$ensID != "-"),]
  
  # Put ensID into rownames
  rownames(dataframe) <- dataframe$ensID
  dataframe <- dataframe[,-which(colnames(dataframe) %in% "ensID")]
  return(dataframe)
}


## Write a gene expression table with geneID in rownames as a CIBERSORT compatible file
write.CIB <- function(data, dir) {
  data <- cbind(rownames(data), data)
  colnames(data)[1] <- "EnsID"
  write.table(data, file = dir, sep = "\t", quote = FALSE, row.names = FALSE)
}

## Write a standard "phenotype class" table for CIBERSORT. This follows the most basic form of every cell-type being compared to every other cell-type
write.CIB.pheno <- function(data, dir, basic = TRUE) {
  n <- ncol(data)
  
  if (basic) {
    pheno <- matrix(data = 2, nrow = n, ncol = n)
    diag(pheno) <- 1
    pheno <- as.data.frame(pheno)
    rownames(pheno) <- colnames(data)
    
    write.table(pheno, col.names = FALSE, sep = "\t", quote = FALSE, file = dir)  
  }
  
}

##################################################################################################################################
## Choose n markers per cell-type within a signature

## Collect the top n markers of each cell-type
nTopFeatures <- function(signature, n, alg = "diff") {
  # the signature is a dataframe of aggregated cell-types x genes
  # n is the number of markers per cell-type, giving a total of n*ncol(signature) markers
  
  # setup
    sig <- t(log2(signature + 1))
    
    ps <- list() # ps is required by dtangle's find_markers() function.  
    for (j in rownames(sig)) { ps[[j]] <- match(j, rownames(sig)) }
    
    res <- as.data.frame(matrix(ncol = 2, nrow = n*nrow(sig)))
    colnames(res) <- c("EnsID", "forCelltype")
    
  # find markers, using find_markers() from the dtangle package
  allMarkers <- find_markers(sig, marker_method = alg, data_type = "rna-seq", pure_samples = ps)
  nMarkers <- sapply(allMarkers$V, function(x) { names(x)[1:n] } ) 

  # add to res
  counter <- 1
  for (j in colnames(nMarkers)) {
    for(k in 1:n) {
      res$EnsID[counter] <- nMarkers[k,j]
      res$forCelltype[counter] <- j
      counter <- counter + 1
    }
  }
  
  # alternate situation when n = 1
  if (n == 1) {
    res$EnsID <- nMarkers
    res$forCelltype <- names(nMarkers)
  }

  # output
  return(res)
}

nRandomFeatures <- function(signature, n, seed) {
  # the signature is a dataframe of aggregated cell-types x genes
  # n is the number of markers per cell-type, giving a total of n*ncol(signature) markers
  
  # setup
  genes <- rownames(signature) # the data of interest
  nMarkers <- n*ncol(signature) # number of markers per celltype * number of celltypes = number of markers
  set.seed(seed)
  
  # collect
  res <- genes[sample(x = 1:length(genes), size = nMarkers, replace = FALSE)]
  
  # output
  return(res)
}



##################################################################################################################################
## scme-related functions

## Random sampling
  create.scmeRand <- function(cellTypes, # celltypes to use
                              nMix, # number of cells to mix per sample
                              nReps, # number of samples to make
                              x = darmanis, # a dataframe of single-cells, where the column name is the celltype
                              rpkm = TRUE) {
    # subset single-cell data setup
    x <- x[,which(colnames(x) %in% cellTypes)]
    
    # setup dataframes
    mixed <- as.data.frame(matrix(nrow = nrow(x), ncol = nReps))
    meta <- as.data.frame(matrix(nrow = nReps, ncol = length(cellTypes)))
    colnames(meta) <- cellTypes
    rownames(mixed) <- rownames(x)
    
    # collect information on celltypes
    colnames(x) <- sapply(strsplit(colnames(x), "\\."), `[`, 1)
    ratios <- table(colnames(x))
    ratios <- ratios[cellTypes]
    
    for(i in 1:nReps) {
      # report progress to user
      print(paste0(i / nReps * 100, "%"))
  
      # sampling vector
      sampleVec <- sample(seq(from = 1, to = ncol(x), by = 1), size = nMix, replace = FALSE)
      
      # sample these single cells, and collect:
        # the average if rpkm, the sum if counts
        if (rpkm) {
          mixed[,i] <- rowMeans(x[,sampleVec])  
        } else {
          mixed[,i] <- rowSums(x[,sampleVec])  
        }
        
      # annotate the output metadata with proportions
        # the number in the mixture if counts
        if (rpkm) {
          ct <- table(colnames(x)[sampleVec])
          for(j in names(ct)) meta[i,j] <- ct[j] / nMix  
        } else {
          for(j in cellTypes) {
            ct <- x[,intersect(sampleVec,grep(j, colnames(x)))]
            meta[i,j] <- sum(ct) / sum(mixed[,i])  
          }
        }
        
    }
    
    # function output
    output <- list()
    output$mixed <- mixed
    output$tProps <- meta
    return(output)
  }

## Gradient sampling
  create.scmeGrad <- function(cellTypes, nMix, nReps, gradCellType) {
    # gradCellType is the celltype for which a gradient from 0-100% is desired
    # other arguments per create.scmeRand
  
    # subset single-cell data setup
    x <- darmanis[,which(colnames(darmanis) %in% cellTypes)]
    
    # setup dataframes
    mixed <- as.data.frame(matrix(nrow = nrow(x), ncol = nReps))
    meta <- as.data.frame(matrix(nrow = nReps, ncol = length(cellTypes)))
    colnames(meta) <- cellTypes
    rownames(mixed) <- rownames(x)
    
    # collect information on celltypes
    colnames(x) <- sapply(strsplit(colnames(x), "\\."), `[`, 1)
    ratios <- table(colnames(x))
    ratios <- ratios[cellTypes]
    nGrad <- ratios[names(ratios) == gradCellType]
    if (nGrad > nMix) nGrad <- nMix
    iGrad <- grep(gradCellType, colnames(x))
    iNotGrad <- (1:ncol(x))[-grep(gradCellType, colnames(x))]
    
    for(i in 1:nReps) {
      # report to user
      print(paste0(i / nReps * 100, "%"))
    
      # pick a random number from 1:nGrad for how many gradCellType cells to take
      n <- sample(seq(from = 1, to = nGrad, by = 1), size = 1, replace = FALSE)
    
      # get the index of a random n gradCellTypes
      n <- sample(iGrad, size = n, replace = FALSE)
      
      # randomly sample 100 - length(n) cells of other origins
      sampleVec <- c(n, sample(iNotGrad, size = nMix-length(n), replace = FALSE))
  
      # sample these single cells, and collect average
      mixed[,i] <- rowMeans(x[,sampleVec])
      
      # annotate the output metadata with proportions
      ct <- table(colnames(x)[sampleVec])
      for(j in names(ct)) meta[i,j] <- ct[j] / nMix
      for(j in names(ct)) meta[is.na(meta[,j]),j] <- 0
    }
    
    # function output
    output <- list()
    output$mixed <- mixed
    output$tProps <- meta
    return(output)
  }
  
  # a new version of create.scmeGrad, aimed at reducing collinearity
  create.scmeGrad2 <- function(cellTypes, gradCellType, nGrad, nNotGrad, nReps) {
    # gradCellType is the celltype for which a gradient from 0-100% is desired
    # other arguments per create.scmeRand
  
    # subset single-cell data setup
    x <- darmanis[,which(colnames(darmanis) %in% cellTypes)]
    
    # setup dataframes
    mixed <- as.data.frame(matrix(nrow = nrow(x), ncol = nReps))
    meta <- as.data.frame(matrix(nrow = nReps, ncol = length(cellTypes)))
    colnames(meta) <- cellTypes
    rownames(mixed) <- rownames(x)
    
    # collect information on celltypes
    colnames(x) <- sapply(strsplit(colnames(x), "\\."), `[`, 1)
    ratios <- table(colnames(x))
    ratios <- ratios[cellTypes]
    # nGrad <- ratios[names(ratios) == gradCellType]
    iGrad <- grep(gradCellType, colnames(x))
    iNotGrad <- (1:ncol(x))[-grep(gradCellType, colnames(x))]
    
    for(i in 1:nReps) {
      # report to user
      print(paste0(i / nReps * 100, "%"))
    
      # pick a random number from 1:nGrad for how many gradCellType cells to take
      n <- sample(seq(from = 1, to = nGrad, by = 1), size = 1, replace = FALSE)
      # print(paste0(i / nReps * 100, "%, nGrad = ", n)) # used in troubleshooting
      nTotal <- n + nNotGrad
    
      # get the index of a random n gradCellTypes
      n <- sample(iGrad, size = n, replace = FALSE)
      
      ## the following line is the (major) point of difference between v1 and v2
      
      # randomly sample nNotGrad cells from the cell-types that aren't gradCellType
      sampleVec <- c(n, sample(iNotGrad, size = nNotGrad, replace = FALSE))
  
      # sample these single cells, and collect average
      mixed[,i] <- rowMeans(x[,sampleVec])
      
      # annotate the output metadata with proportions
      ct <- table(colnames(x)[sampleVec])
      for(j in names(ct)) meta[i,j] <- ct[j] / nTotal
      for(j in names(ct)) meta[is.na(meta[,j]),j] <- 0
    }
    
    # function output
    output <- list()
    output$mixed <- mixed
    output$tProps <- meta
    return(output)
  }
  
##################################################################################################################################
## snme-related functions
  
## Random sampling
create.snmeRand <- function(nMix, # number of cells to mix per sample
                            nReps, # number of samples to make
                            x = darmanis, # a dataframe of single-cells
                            meta, # a dataframe of celltype labels
                            meta.columns, # columns in the dataframe to use as true labels during mixing
                            rpkm = TRUE) {
  
  # recalculate library size in the meta data
  meta$libSize <- colSums(x)
  
  # setup mixture dataframe
  mixed <- as.data.frame(matrix(nrow = nrow(x), ncol = nReps))
  rownames(mixed) <- rownames(x)
  
  # setup true-proportion dataframe
  true <- list()
  for(j in meta.columns) {
    ct <- levels(as.factor(meta[,j]))
    true[[j]] <- as.data.frame(matrix(nrow = nReps, ncol = length(ct)))
    colnames(true[[j]]) <- ct
  }
  
  # for loop for mixture generation
  for(i in 1:nReps) {
    # report progress to user
    print(paste0(i / nReps * 100, "%"))

    # sampling vector
    sampleVec <- sample(colnames(x), size = nMix, replace = FALSE)
    
    # sample these single cells, and collect:
     
      if (rpkm) {  # the average if rpkm...
        mixed[,i] <- rowMeans(x[,sampleVec])  
      } else { # and the sum if counts
        mixed[,i] <- rowSums(x[,sampleVec])  
      }
      
    # annotate the output metadata with proportions
      
      if (rpkm) { # the number in the mixture if rpkm
        
        for(j in meta.columns) {
          ct <- table(meta[sampleVec,j])
          ct <- ct / nMix
          ct <- ct[colnames(true[[j]])]
          true[[j]][i,] <- ct
        }
          
      } else { # the ratio of reads if counts
        
        for(j in meta.columns) {
          meta$temp <- meta[,j]
          ratio <- aggregate(libSize~temp, meta[sampleVec,], sum)
          ratio$libSize <- ratio$libSize / sum(ratio$libSize)
          rownames(ratio) <- ratio$temp
          ratio <- ratio[colnames(true[[j]]),]
          
          true[[j]][i,] <- ratio$libSize
          # ct <- x[,intersect(sampleVec,grep(j, colnames(x)))]
          # true[i,j] <- sum(ct) / sum(mixed[,i])  
        }
        
      }
      
  }
  
  # in the true dataframe, the absence of a celltype is coded as NA; convert this to 0
  true <- lapply(true, function(x) {
    x <- apply(x, 2, function(y) {
      z <- which(is.na(y))
      y[z] <- 0
      return(y)
    })
    return(x)
  })
  
  # create a cpm version of the mixture
  cpm <- apply(mixed, 2, function(x) {
    lib.size <- 10^6 / sum(x)
    x <- x * lib.size  
    return(x)
  })
  cpm <- as.data.frame(cpm)
  
  # function output
  output <- list()
  output$mixture <- list(counts = mixed, cpm = cpm)
  output$true <- true
  return(output)
}

##################################################################################################################################
## Differential expression simulation

## I want a function that can mix as many scme cell-types as I'd like...
## ...with each having a defined range of proportions in Group1, and in Group2...
## ...and then I can independently simulate adding a few RPKM to some genes, and taking away from others...
## ...but given that I'm mixing RPKM, DE might have to be performed by Limma...

# confound.proportion <- function(ct, # a character vector of celltypes to use in the mixture
#                                 setProps, # a dataframe of proportion ranges to use any of the celltypes in ct. Ct not specified will be randomly sampled
#                                 nPerGroup = 50, # the number of mixtures to make per group (and there will be two groups...)
#                                 genesUp = FALSE, # a character vector of genes to be upregulated
#                                 genesDown = FALSE,  # a character vector of genes to be downregulated
#                                 absoluteFC = 1.3) { # the desired fold-change to up- or down-regulated genes
#   
#   ## Filter
#   dat <- darmanis[,which(colnames(darmanis) %in% ct)]
#   colnames(dat) <- sapply(strsplit(colnames(dat), "\\."), "[", 1)
#   
#   ## Find existing proportions
#   available <- table(colnames(dat))
#   
#   ## Determine which subset of ct has been supplied setProps ("used")
#   used <- ct[which(ct %in% colnames(setProps))]
#   usedIndex <- which(colnames(dat) %in% used)
# 
#   ## Define dataframes
#   mixed <- as.data.frame(matrix(nrow = nrow(dat), ncol = 2*nPerGroup))
#   meta <- as.data.frame(matrix(nrow = 2*nPerGroup, ncol = length(ct)))
#   colnames(meta) <- ct
#   rownames(mixed) <- rownames(dat)
#   colnames(mixed) <- rownames(meta) <- c(paste0("Alpha_", 1:nPerGroup), paste0("Beta_", 1:nPerGroup))
# 
#   ## Collect indices of celltypes
#   indices <- list()
#   for(j in names(available)) indices[[j]] <- grep(j, colnames(dat))
#   
#   ## Sample cells at different proportions between groups
#   for(i in 1:ncol(mixed)) {
#     # report to user
#     print(paste0(i / (2*nPerGroup) * 100, "%"))
#     
#     # define sampleVec as a vector holding the index of all cells to take
#     sampleVec <- c()
# 
#     # collect proportions for every celltype with a specified (i.e. non-random) range
#     if(i <= nPerGroup) {
#       for (j in used) { 
#         # pick a random number of cells in the specified range for Group1
#         n <- sample(seq(from = setProps[1,j], to = setProps[2,j], by = 1), size = 1, replace = FALSE)
#         
#         # get the index of a random n gradCellTypes
#         sampleVec <- c(sampleVec, sample(indices[[j]], size = n, replace = FALSE))
#       }
#     } else {
#       for (j in used) { 
#         # pick a random number of cells in the specified range for Group1
#         n <- sample(seq(from = setProps[3,j], to = setProps[4,j], by = 1), size = 1, replace = FALSE)
#         
#         # get the index of a random n gradCellTypes
#         sampleVec <- c(sampleVec, sample(indices[[j]], size = n, replace = FALSE))
#       }
#     }
# 
#     # now sample other celltypes
#     sampleVec <- c(sampleVec, sample((1:ncol(dat))[-usedIndex], size = 100 - length(sampleVec), replace = FALSE))
#       
#     # sample these single cells, and collect average
#     mixed[,i] <- rowMeans(dat[,sampleVec])
#     
#     # annotate the output metadata with proportions
#     annot <- table(colnames(dat)[sampleVec])
#     for(j in names(annot)) meta[i,j] <- annot[j] / 100
#     for(j in names(annot)) meta[is.na(meta[,j]),j] <- 0
#   }
#   
#   ## Output
#   res <- list(mixed = mixed, meta = meta)
#   return(res)
# }
# 
# confound.expression <- function(confounded.proportions, genesUp, genesDown, absoluteFC = 1.3, noise = 0.3) {
#   # confounded.proportions: a list generated using confound.proportions
#   # genesUp: a dataframe of genes to be upregulated. $EnsID denotes the gene, and $forCelltype is its annotation
#   # genesDown: per genesUp, but genes here are downregulated
#   # absoluteFC = 1.3: the desired fold-change to up- or down-regulated genes, not in log scale!
#   # noise = 0.3: sets the fc in each sample as between  70% and 130% of the absoluteFC, but mean = absoluteFC
#   
#   # rename for convenience
#   x <- confounded.proportions
#   
#   # auto-detect group size in confounded.proportions "alpha" as the first nPerGroup samples
#   alpha <- grep("Alpha", colnames(x$mixed))  
# 
#   # upregulate genes
#   for (j in genesUp$EnsID) {
#     # generate a vector of fold changes with added noise
#     noisyFC <- runif(n = nPerGroup, min = absoluteFC-noise, max = absoluteFC+noise)
#     
#     # multiply the expression of gene j by absoluteFC * geneNoise
#     x$mixed[j,alpha] <- x$mixed[j,alpha]*noisyFC
#   }
#   
#   # downregulate genes
#   for (j in genesDown$EnsID) {
#     # generate a vector of fold changes with added noise
#     noisyFC <- 1 / runif(n = nPerGroup, min = absoluteFC-noise, max = absoluteFC+noise)
#     
#     # multiply the expression of gene j by absoluteFC * geneNoise
#     x$mixed[j,alpha] <- x$mixed[j,alpha]*noisyFC
#   }
# 
#   
#   # add dataframe revealing which genes were upregulated, and which were downregulated
#   x$de <- rbind(genesUp, genesDown)
#   x$de$Direction <- rep(c("Up", "Down"), each = nrow(genesUp))
#   
#   # add metadata
#   x$de.meta <- data.frame(absoluteFC = absoluteFC, noise = noise)
#   rownames(x$de.meta) <- "Value"
#   
#   # output
#   return(x)
# }

# ## V2!
# confound.proportion <- function(ct, # a character vector of celltypes to use in the mixture
#                                 setProps, # a dataframe of proportion ranges to use any of the celltypes in ct. Ct not specified will be randomly sampled
#                                 dat = darmanis, # the dataframe from which to sample cells, with colnames as the celltype
#                                 nPerGroup = 50) { # the number of mixtures to make per group (and there will be two groups...)
#                                 
#   ## Filter
#   dat <- dat[,which(colnames(dat) %in% ct)]
#   colnames(dat) <- sapply(strsplit(colnames(dat), "\\."), "[", 1)
#   
#   ## Find existing proportions
#   available <- table(colnames(dat))
#   
#   ## Determine which subset of ct has been supplied setProps ("used")
#   used <- ct[which(ct %in% colnames(setProps))]
#   usedIndex <- which(colnames(dat) %in% used)
# 
#   ## Define dataframes
#   mixed <- as.data.frame(matrix(nrow = nrow(dat), ncol = 2*nPerGroup))
#   meta <- as.data.frame(matrix(nrow = 2*nPerGroup, ncol = length(ct)))
#   colnames(meta) <- ct
#   rownames(mixed) <- rownames(dat)
#   colnames(mixed) <- rownames(meta) <- c(paste0("Alpha_", 1:nPerGroup), paste0("Beta_", 1:nPerGroup))
# 
#   ## Collect indices of celltypes
#   indices <- list()
#   for(j in ct) indices[[j]] <- grep(j, colnames(dat))
#   
#   ## Hold cell-type specific expression
#   deconvolved.expression <- list()
#   for (j in ct) deconvolved.expression[[j]] <- mixed
#   
#   ## Sample cells at different proportions between groups
#   for(i in 1:ncol(mixed)) {
#     # report to user
#     print(paste0(i / (2*nPerGroup) * 100, "%"))
#     
#     # define sampleVec as a vector holding the index of all cells to take
#     sampleVec <- c()
# 
#     # collect proportions for every celltype with a specified (i.e. non-random) range
#     if(i <= nPerGroup) {
#       for (j in used) { 
#         # pick a random number of cells in the specified range for Group1
#         n <- sample(seq(from = setProps[1,j], to = setProps[2,j], by = 1), size = 1, replace = FALSE)
#         
#         # get the index of a random n gradCellTypes
#         sampleVec <- c(sampleVec, sample(indices[[j]], size = n, replace = FALSE))
#       }
#     } else {
#       for (j in used) { 
#         # pick a random number of cells in the specified range for Group1
#         n <- sample(seq(from = setProps[3,j], to = setProps[4,j], by = 1), size = 1, replace = FALSE)
#         
#         # get the index of a random n gradCellTypes
#         sampleVec <- c(sampleVec, sample(indices[[j]], size = n, replace = FALSE))
#       }
#     }
# 
#     # now sample other celltypes
#     sampleVec <- c(sampleVec, sample((1:ncol(dat))[-usedIndex], size = 100 - length(sampleVec), replace = FALSE))
#       
#     # sample these single cells, and collect average
#     mixed[,i] <- rowMeans(dat[,sampleVec])
#     
#     # for each cell-type, separately collect the average expression that goes into the mixture
#     for (j in ct) {
#       if (length(which(indices[[j]] %in% sampleVec)) >= 2) {
#         g <- indices[[j]][which(indices[[j]] %in% sampleVec)]
#         deconvolved.expression[[j]][,i] <- rowMeans(dat[,g])  
#       } 
#       if (length(which(indices[[j]] %in% sampleVec)) == 1) {
#         g <- indices[[j]][which(indices[[j]] %in% sampleVec)]
#         deconvolved.expression[[j]][,i] <- dat[,g]
#       }
#       if (length(which(indices[[j]] %in% sampleVec)) == 0) {
#         deconvolved.expression[[j]][,i] <- 0
#       }
#     }
#     
#     # annotate the output metadata with proportions
#     annot <- table(colnames(dat)[sampleVec])
#     for(j in ct) meta[i,j] <- annot[j] / 100
#     for(j in ct) meta[is.na(meta[,j]),j] <- 0
#   }
#   
#   ## Output
#   res <- deconvolved.expression
#   res$mixed <- mixed
#   res$meta.prop <- meta
#   return(res)
# }
# 
# confound.expression <- function(confounded.proportions, cell.type = "mixed", genesUp, genesDown, absoluteFC = 1.3, noise = 0.3, id = "new") {
#   # confounded.proportions: a list generated using confound.proportions
#   # cell.type: sets the level of confounded.proportions to be perturbed. to perturb multiple levels, run the function multiple times
#   # genesUp: a dataframe of genes to be upregulated. $EnsID denotes the gene, and $forCelltype is its annotation
#   # genesDown: per genesUp, but genes here are downregulated
#   # absoluteFC = 1.3: the desired fold-change to up- or down-regulated genes, not in log scale!
#   # noise = 0.3: sets the fc in each sample as between  70% and 130% of the absoluteFC, but mean = absoluteFC
#   # id = "new": this sets the name of the new levels to be added to confounded.proportions
#   
#   # a tally
#   print("New call started")
#   
#   # rename for convenience
#   x <- confounded.proportions
#   new.name <- paste0(id, ".", cell.type)
#   x[[new.name]] <- x[[cell.type]]
#   
#   # auto-detect group size in confounded.proportions "alpha" as the first nPerGroup samples
#   alpha <- grep("Alpha", colnames(x[[cell.type]]))  
# 
#   # upregulate genes
#   for (j in genesUp$EnsID) {
#     # generate a vector of fold changes with added noise
#     noisyFC <- runif(n = nPerGroup, min = absoluteFC-noise, max = absoluteFC+noise)
#     
#     # multiply the expression of gene j by absoluteFC * geneNoise
#     x[[new.name]][j,alpha] <- x[[new.name]][j,alpha]*noisyFC
#   }
#   
#   # downregulate genes
#   for (j in genesDown$EnsID) {
#     # generate a vector of fold changes with added noise
#     noisyFC <- 1 / runif(n = nPerGroup, min = absoluteFC-noise, max = absoluteFC+noise)
#     
#     # multiply the expression of gene j by absoluteFC * geneNoise
#     x[[new.name]][j,alpha] <- x[[new.name]][j,alpha]*noisyFC
#   }
#   
#   if (cell.type != "mixed") { # reconstitute expression, if perturbing anything by the "mixed" matrix (because that's already reconstituted)
#     # create new level in x to hold reconstituted expression
#     newer.name <- paste0(new.name, ".reconstituted")
#     x[[newer.name]] <- x[[new.name]]
#     x[[newer.name]][] <- 0 # a zero matrix
#     
#     
#     for (j in 1:ncol(x$mixed)) { # for each sample
#       for (k in colnames(x$meta)) { # take every cell-type
#         if(k == cell.type) {
#           x[[newer.name]][,j] <- x[[newer.name]][,j] + (x[[new.name]][,j] * x$meta.prop[j,k])
#         } else { # if the cell-type is the pertubed one, take the pertubed matrix... else use the non-pertubed one!
#           x[[newer.name]][,j] <- x[[newer.name]][,j] + (x[[k]][,j] * x$meta.prop[j,k])
#         }
#         # that is, reconstituted expression is the average of (perturbed) expression in each cell-type, weighted by respective cell-type abundance
#         # note that the initial value is 0
#       }
#     }
#   }
# 
#   # leave a trail of what perturbations were made
#       if ("meta.exp" %in% names(x)) {
#         y <- x[["meta.exp"]] # if it already exists (i.e., if the function has been applied to confounded.proportions before), append new metadata
#       } else {
#         y <- list()
#       }
#     
#     y[[new.name]] <- list()
#     
#     # add dataframe revealing which genes were upregulated, and which were downregulated
#     y[[new.name]]$de <- rbind(genesUp, genesDown)
#     y[[new.name]]$de$Direction <- rep(c("Up", "Down"), each = nrow(genesUp))
#     
#     # perturbation in fc and the inherent noise
#     y[[new.name]]$de.meta <- data.frame(absoluteFC = absoluteFC, noise = noise)
#     rownames(y[[new.name]]$de.meta) <- "Value"
# 
#     # append
#     x$meta.exp <- y
#     
#   # output
#   return(x)
# }
  
# ## V3!
#   # this essentially updates the function to streamline it, removing the raw output of non-perturbed celltypes
#   confound.proportion <- function(ct, # a character vector of celltypes to use in the mixture
#                                   setProps, # a dataframe of proportion ranges to use any of the celltypes in ct. Ct not specified will be randomly sampled
#                                   dat, # the dataframe from which to sample cells, with colnames as the celltype
#                                   reduce.filesize = TRUE,
#                                   nPerGroup) { # the number of mixtures to make per group (and there will be two groups...)
#     
#     ## Filter
#     dat <- dat[,which(colnames(dat) %in% ct)]
#     colnames(dat) <- sapply(strsplit(colnames(dat), "\\."), "[", 1)
#     
#     ## Find existing proportions
#     available <- table(colnames(dat))
#     
#     ## Determine which subset of ct has been supplied setProps ("used")
#     used <- ct[which(ct %in% colnames(setProps))]
#     usedIndex <- which(colnames(dat) %in% used)
#     
#     ## Define dataframes
#     mixed <- as.data.frame(matrix(nrow = nrow(dat), ncol = 2*nPerGroup))
#     meta <- as.data.frame(matrix(nrow = 2*nPerGroup, ncol = length(ct)))
#     colnames(meta) <- ct
#     rownames(mixed) <- rownames(dat)
#     colnames(mixed) <- rownames(meta) <- c(paste0("Alpha_", 1:nPerGroup), paste0("Beta_", 1:nPerGroup))
#     
#     ## Collect indices of celltypes
#     indices <- list()
#     for(j in ct) indices[[j]] <- grep(j, colnames(dat))
#     
#     ## Hold cell-type specific expression
#     deconvolved.expression <- list()
#     for (j in ct) deconvolved.expression[[j]] <- mixed
#     if (reduce.filesize) deconvolved.expression <- deconvolved.expression[colnames(setProps)]
#     
#     
#     ## Sample cells at different proportions between groups
#     for(i in 1:ncol(mixed)) {
#       # report to user
#       print(paste0(i / (2*nPerGroup) * 100, "%"))
#       
#       # define sampleVec as a vector holding the index of all cells to take
#       sampleVec <- c()
#       
#       # collect proportions for every celltype with a specified (i.e. non-random) range
#       if(i <= nPerGroup) {
#         for (j in used) { 
#           # pick a random number of cells in the specified range for Group1
#           n <- sample(seq(from = setProps[1,j], to = setProps[2,j], by = 1), size = 1, replace = FALSE)
#           
#           # get the index of a random n gradCellTypes
#           sampleVec <- c(sampleVec, sample(indices[[j]], size = n, replace = FALSE))
#         }
#       } else {
#         for (j in used) { 
#           # pick a random number of cells in the specified range for Group1
#           n <- sample(seq(from = setProps[3,j], to = setProps[4,j], by = 1), size = 1, replace = FALSE)
#           
#           # get the index of a random n gradCellTypes
#           sampleVec <- c(sampleVec, sample(indices[[j]], size = n, replace = FALSE))
#         }
#       }
#       
#       # now sample other celltypes
#       sampleVec <- c(sampleVec, sample((1:ncol(dat))[-usedIndex], size = 100 - length(sampleVec), replace = FALSE))
#       
#       # sample these single cells, and collect average
#       mixed[,i] <- rowMeans(dat[,sampleVec])
#       
#       # for each cell-type, separately collect the average expression that goes into the mixture
#       for (j in names(deconvolved.expression)) {
#         
#         if (length(which(indices[[j]] %in% sampleVec)) >= 2) {
#           g <- indices[[j]][which(indices[[j]] %in% sampleVec)]
#           deconvolved.expression[[j]][,i] <- rowMeans(dat[,g])  
#         } 
#         
#         if (length(which(indices[[j]] %in% sampleVec)) == 1) {
#           g <- indices[[j]][which(indices[[j]] %in% sampleVec)]
#           deconvolved.expression[[j]][,i] <- dat[,g]
#         }
#         
#         if (length(which(indices[[j]] %in% sampleVec)) == 0) {
#           deconvolved.expression[[j]][,i] <- 0
#         }
#         
#       }
#       
#       # annotate the output metadata with proportions
#       annot <- table(colnames(dat)[sampleVec])
#       for(j in ct) meta[i,j] <- annot[j] / 100
#       for(j in ct) meta[is.na(meta[,j]),j] <- 0
#     }
#     
#     ## Output
#     res <- deconvolved.expression
#     res$mixed <- mixed
#     res$meta.prop <- meta
#     return(res)
#   }
#   
#   confound.expression <- function(confounded.proportions, cell.type = "mixed", genesUp, genesDown, absoluteFC = 1.3, noise = 0.3, id = "new") {
#     # confounded.proportions: a list generated using confound.proportions
#     # cell.type: sets the level of confounded.proportions to be perturbed. to perturb multiple levels, run the function multiple times
#     # genesUp: a dataframe of genes to be upregulated. $EnsID denotes the gene, and $forCelltype is its annotation
#     # genesDown: per genesUp, but genes here are downregulated
#     # absoluteFC = 1.3: the desired fold-change to up- or down-regulated genes, not in log scale!
#     # noise = 0.3: sets the fc in each sample as between  70% and 130% of the absoluteFC, but mean = absoluteFC
#     # id = "new": this sets the name of the new levels to be added to confounded.proportions
#     
#     # a tally
#     print("New call started")
#     
#     # rename for convenience
#     x <- confounded.proportions
#     new.name <- paste0(id, ".", cell.type)
#     x[[new.name]] <- x[[cell.type]]
#     
#     # auto-detect group size in confounded.proportions "alpha" as the first nPerGroup samples
#     alpha <- grep("Alpha", colnames(x[[cell.type]]))  
#     
#     # upregulate genes
#     for (j in genesUp$EnsID) {
#       # generate a vector of fold changes with added noise
#       noisyFC <- runif(n = nPerGroup, min = absoluteFC-noise, max = absoluteFC+noise)
#       
#       # multiply the expression of gene j by absoluteFC * geneNoise
#       x[[new.name]][j,alpha] <- x[[new.name]][j,alpha]*noisyFC
#     }
#     
#     # downregulate genes
#     for (j in genesDown$EnsID) {
#       # generate a vector of fold changes with added noise
#       noisyFC <- 1 / runif(n = nPerGroup, min = absoluteFC-noise, max = absoluteFC+noise)
#       
#       # multiply the expression of gene j by absoluteFC * geneNoise
#       x[[new.name]][j,alpha] <- x[[new.name]][j,alpha]*noisyFC
#     }
#     
#     if (cell.type != "mixed") { # reconstitute expression, if perturbing anything but the "mixed" matrix (because that's already reconstituted)
#       # create new level in x to hold reconstituted expression
#       newer.name <- paste0(new.name, ".reconstituted")
#       x[[newer.name]] <- x[[new.name]]
#       x[[newer.name]][] <- 0 # a zero matrix
#       
#       
#       for (j in 1:ncol(x$mixed)) { # for each sample
#         for (k in colnames(x$meta)) { # take every cell-type
#           if(k == cell.type) {
#             x[[newer.name]][,j] <- x[[newer.name]][,j] + (x[[new.name]][,j] * x$meta.prop[j,k])
#           } else { # if the cell-type is the pertubed one, take the pertubed matrix... else use the non-pertubed one!
#             x[[newer.name]][,j] <- x[[newer.name]][,j] + (x[[k]][,j] * x$meta.prop[j,k])
#           }
#           # that is, reconstituted expression is the average of (perturbed) expression in each cell-type, weighted by respective cell-type abundance
#           # note that the initial value is 0
#         }
#       }
#     }
#     
#     # leave a trail of what perturbations were made
#     if ("meta.exp" %in% names(x)) {
#       y <- x[["meta.exp"]] # if it already exists (i.e., if the function has been applied to confounded.proportions before), append new metadata
#     } else {
#       y <- list()
#     }
#     
#     y[[new.name]] <- list()
#     
#     # add dataframe revealing which genes were upregulated, and which were downregulated
#     y[[new.name]]$de <- rbind(genesUp, genesDown)
#     y[[new.name]]$de$Direction <- rep(c("Up", "Down"), each = nrow(genesUp))
#     
#     # perturbation in fc and the inherent noise
#     y[[new.name]]$de.meta <- data.frame(absoluteFC = absoluteFC, noise = noise)
#     rownames(y[[new.name]]$de.meta) <- "Value"
#     
#     # append
#     x$meta.exp <- y
#     
#     # NEW BIT OF CODE TO REDUCE RAM BURDER
#     x <- x[c("mixed", "meta.prop", new.name, newer.name, "meta.exp")]
#     
#     # output
#     return(x)
#   }
  
## V4!
  ## Changelog
    ## Increased the number of cells per mixture to 500 (from 100)
    ## No need to save $mixed, as all models will be run in the perturbed data
  
  
  confound.proportion <- function(ct, # a character vector of celltypes to use in the mixture
                                  setProps, # a dataframe of proportion ranges to use any of the celltypes in ct. Ct not specified will be randomly sampled
                                  dat, # the dataframe from which to sample cells, with colnames as the celltype
                                  nPerMixture = 500,
                                  keepPureExpression = FALSE,
                                  nPerGroup) { # the number of mixtures to make per group (and there will be two groups...)
    
  ## Filter
    dat <- dat[,which(colnames(dat) %in% ct)]
    colnames(dat) <- sapply(strsplit(colnames(dat), "\\."), "[", 1)
    
  ## Find existing proportions
    available <- table(colnames(dat))
    libSizes <- colSums(dat)
    
  ## Determine which subset of ct has been supplied setProps ("used")
    used <- ct[which(ct %in% colnames(setProps))]
    usedIndex <- which(colnames(dat) %in% used)
    
  ## Define dataframes
    mixed <- as.data.frame(matrix(nrow = nrow(dat), ncol = 2*nPerGroup))
    meta <- as.data.frame(matrix(nrow = 2*nPerGroup, ncol = length(ct)))
    colnames(meta) <- ct
    rownames(mixed) <- rownames(dat)
    colnames(mixed) <- rownames(meta) <- c(paste0("Alpha_", 1:nPerGroup), paste0("Beta_", 1:nPerGroup))
    
  ## Collect indices of celltypes
    indices <- list()
    for(j in ct) indices[[j]] <- grep(j, colnames(dat))
    
  ## List to store cell-type-specific expression
    deconvolved.expression <- list()
    for (j in ct) deconvolved.expression[[j]] <- mixed
    
  ## Report on the simulation
    cat("\n")
    print(paste0("Simulating a dataset with ", used, " changed by ", (setProps[3,] - setProps[1,]) / nPerMixture, "%"))
    
  ## Sample cells at different proportions between groups
    for(i in 1:ncol(mixed)) {
      # further report to user
      cat("Simulating sample", i, "of", (2*nPerGroup), "\r")
      flush.console()
      
      # define sampleVec as a vector holding the index of all cells to take
      sampleVec <- c()
      
      # collect proportions for every celltype with a specified (i.e. non-random) range
      if(i <= nPerGroup) {
        for (j in used) { 
          # pick a random number of cells in the specified range for Group1
          n <- sample(seq(from = setProps[1,j], to = setProps[2,j], by = 1), size = 1, replace = FALSE)
          
          # get the index of a random n gradCellTypes
          sampleVec <- c(sampleVec, sample(indices[[j]], size = n, replace = FALSE))
        }
      } else {
        for (j in used) { 
          # pick a random number of cells in the specified range for Group1
          n <- sample(seq(from = setProps[3,j], to = setProps[4,j], by = 1), size = 1, replace = FALSE)
          
          # get the index of a random n gradCellTypes
          sampleVec <- c(sampleVec, sample(indices[[j]], size = n, replace = FALSE))
        }
      }
      
      ## Now sample other celltypes
        if (length(sampleVec) > nMix) {
          
          sampleVec <- sample(sampleVec, size = nMix, replace = FALSE) # if more than the specified number of cells has been sampled, reduce to that number
          
        } else { # add more
          if (length(which(names(available) %in% used)) != length(available)) { # if there are ct without specified abundance ranges, sample from those
            sampleVec <- c(sampleVec, 
                           sample((1:ncol(dat))[-usedIndex], size = nPerMixture - length(sampleVec), replace = FALSE))  
          } else { # if an abundance range is specified for each ct, then sample from any celltype!
            u <- 1:ncol(dat)
            u <- u[-which(u %in% sampleVec)] # this is the index of all cells which are yet to be sampled
            sampleVec <- c(sampleVec, 
                           sample(u, size = nPerMixture - length(sampleVec), replace = FALSE))  
          }
        }
          
      
      
      # sample these single cells, and collect average for a pseudobulk mixed sample
      mixed[,i] <- rowSums(dat[,sampleVec]) # using rowSums instead of rowMeans enables the use of count data
      
      ## For each cell-type, separately collect the average expression that goes into the mixture
        for (j in names(deconvolved.expression)) {
          
          if (length(which(indices[[j]] %in% sampleVec)) >= 2) {
            g <- indices[[j]][which(indices[[j]] %in% sampleVec)]
            deconvolved.expression[[j]][,i] <- rowSums(dat[,g])
          }
          
          if (length(which(indices[[j]] %in% sampleVec)) == 1) {
            g <- indices[[j]][which(indices[[j]] %in% sampleVec)]
            deconvolved.expression[[j]][,i] <- dat[,g]
          }
          
          if (length(which(indices[[j]] %in% sampleVec)) == 0) {
            deconvolved.expression[[j]][,i] <- 0
          }
          
        }
      
      # annotate the output metadata with proportions
      ratio <- aggregate(libSizes[sampleVec]~names(libSizes)[sampleVec], FUN = base::sum)
      colnames(ratio) <- c("ct", "libSize")
      rownames(ratio) <- ratio$ct
      ratio$libSize <- ratio$libSize / sum(ratio$libSize)
      meta[i,] <- 0 # the reason this is coded like this is 
      for (j in rownames(ratio)) { meta[i,j] <- ratio[j,"libSize"] }
    }
    
    ## Output
      output <- list()
      output$confound.p <- mixed
      if (keepPureExpression) output$pure.exp <- deconvolved.expression
      output$meta.prop <- meta
      return(output)
  }
  
  confound.expression <- function(x, starting.level = "confound.p", replace.starting.level = TRUE, ct.specific = "all", genesUp, genesDown, absoluteFC = 1.3, noise = 0.1, id = "confound.pe") {
    # x: a list generated using confound.proportions
    # cell.type: sets the level of confounded.proportions to be perturbed. to perturb multiple levels, run the function multiple times
    # replace.starting.level: get rid of the starting.level from the output, useful to save space
    # ct: perturb expression in which cell-type (either all, or specify a ct!)
    # genesUp: a dataframe of genes to be upregulated. $EnsID denotes the gene, and $forCelltype is its annotation
    # genesDown: per genesUp, but genes here are downregulated
    # absoluteFC = 1.3: the desired fold-change to up- or down-regulated genes, not in log scale!
    # noise = 0.3: sets the fc in each sample as between  70% and 130% of the absoluteFC, but mean = absoluteFC
    # id = "confound.pe": this sets the name of the new levels to be added to confounded.proportions
    
    # setup
    if (ct.specific == "all") {
      cat(" Confounding expression by", absoluteFC, "±", noise, ".")
      x[[id]] <- x[[starting.level]]
    } else {
      cat(" Confounding expression by", absoluteFC, "±", noise, "in", ct.specific, ".")
      x[[id]] <- x$pure.exp
    }
    
    # auto-detect group size in confounded.proportions "alpha" as the first nPerGroup samples
    alpha <- grep("Alpha", colnames(x[[starting.level]]))  
    
    # upregulate genes
    for (j in genesUp$EnsID) {
      # generate a vector of fold changes with added noise
      noisyFC <- runif(n = nPerGroup, min = absoluteFC-noise, max = absoluteFC+noise)
      
      # multiply the expression of gene j by absoluteFC * geneNoise
      if (ct.specific == "all") { 
        x[[id]][j,alpha] <- round(x[[id]][j,alpha] * noisyFC) # round ensures and integer; useful for counts  
      } else { 
        x[[id]][[ct.specific]][j,alpha] <- round(x[[id]][[ct.specific]][j,alpha] * noisyFC)
      }
    }
    
    # downregulate genes
    for (j in genesDown$EnsID) {
      # generate a vector of fold changes with added noise
      noisyFC <- 1 / runif(n = nPerGroup, min = absoluteFC-noise, max = absoluteFC+noise)
      
      # multiply the expression of gene j by absoluteFC * geneNoise
      if (ct.specific == "all") { 
        x[[id]][j,alpha] <- round(x[[id]][j,alpha] * noisyFC) # round ensures and integer; useful for counts  
      } else { 
        x[[id]][[ct.specific]][j,alpha] <- round(x[[id]][[ct.specific]][j,alpha] * noisyFC)
      }
    }
    
    # reconstitute expression, if perturbing anything but the "mixed" matrix (because that's already reconstituted)
      if (ct.specific != "all") { 
        # create new level in x to hold reconstituted expression
        x[[id]] <- x[[id]]$Astrocytes + x[[id]]$Excitatory + x[[id]]$Inhibitory + x[[id]]$Oligodendrocytes + x[[id]]$OPC
      }
    
    
    # leave a trail of what perturbations were made
    if ("meta.exp" %in% names(x)) {
      y <- x[["meta.exp"]] # if it already exists (i.e., if the function has been applied to confounded.proportions before), append new metadata
    } else {
      y <- list()
    }
    
    y[[id]] <- list()
    
    # add dataframe revealing which genes were upregulated, and which were downregulated
    y[[id]]$genes <- rbind(genesUp, genesDown)
    y[[id]]$genes$Direction <- rep(c("Up", "Down"), each = nrow(genesUp))
    
    # perturbation in fc and the inherent noise
    y[[id]]$sim.parameters <- data.frame(absoluteFC = absoluteFC, noise = noise)
    rownames(y[[id]]$sim.parameters) <- "Value"
    
    # append
    x$meta.exp <- y
    
    # trim to reduce ram burden
    if (replace.starting.level) {
      w <- which(names(x) == starting.level)
      x <- x[-w]
    } 
  
    # output
    return(x)
  }
  
 