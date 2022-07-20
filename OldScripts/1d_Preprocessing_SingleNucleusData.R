## In this script, mixtures are made from single-nucleus data

################################################################################################################################ #
## Setup ----

## Generic
  rm(list=ls())
  options(stringsAsFactors = FALSE)

## Set directory
  setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/")

## Functions and libraries
  source("../Scripts/Fun_Preprocessing.R")
  source("../Scripts/Fun_Composition.R")
  load(file = "Preprocessed/geneInfo.rda")
  load("Preprocessed/exonicLength.rda")
  library(Seurat)

## Lists
  # to hold the dataset-level Seurat objects
  obj <- list() # to store data from...
  obj$VL <- list() # Velmeshev 2019
  obj$NG <- list() # Nagy 2020
  obj$CA <- list() # Hodge 2019
  obj$TS <- list() # Tasic 2018
  obj$LK <- list() # Lake 2018
  

## Key parameters
  # for seurat preprocessing
  min.cells <- 3 # during the initial load, a gene is excluded if in < 3 cells 
  min.features <- 200 # during the initial load, a barcode is excluded < 200 features are expressed
  min.depth <- 1000 # a barcode is excluded if nCount_RNA < this value
  max.depth.percentile <- 0.995 # a barcode is excluded if nCount_RNA > this percentile within the dataset
  max.mito <- 5
  min.celltype.n <- 200 # minimum number of members in a celltype for it to be kept. applied to anything used for creating mixtures (at this stage, Vel and HCA, but the former passes this criterion for all celltypes anyway...)
  
  # preprocessing options
  downsample <- FALSE
  downsample.n <- NA; if (downsample) downsample.n <- NA
  use.SCTransform <- FALSE
  
  # mixture options
  nMix <- 500 # number of cells to aggregate per mixture
  nReps <- 100 # number of mixtures to simulate

## Functions
  ## Function for downsampling the dataset to a set number of barcodes
    downsample.fun <- function(x, n = downsample.n) {
      if (ncol(x) <= downsample.n) {
        print("No downsampling performed (Reason: number of cells in dataset is already less than or equal to the downsampling number)")
      } else {
        sample <- sample(colnames(x), size = n, replace = FALSE)
        x <- subset(x, cells = sample)  
      }
      return(x)
    }
  
  ## General function for preprocessing sn data (normalise, filters, and scales)
    get.max.depth <- function(x) {
      max.depth <- quantile(x@meta.data$nCount_RNA, probs = max.depth.percentile)
    }
    
    preprocess.fun <- function(x, run.downsample = downsample, SCTransform = use.SCTransform, max.depth = max.depth) {
      # quantify mitochondrial reads
      x[["percent.mito"]] <- PercentageFeatureSet(object = x, pattern = "^MT-")
      
      # filter to remove outlier nuclei: 
      
      x <- subset(x = x, subset = (nCount_RNA > min.depth) & (nCount_RNA < max.depth) & (percent.mito < max.mito))
      
      # downsample
      if (run.downsample) { x <- downsample.fun(x) }
      
      # normalise expression levels
      x <- NormalizeData(object = x, normalization.method = "LogNormalize", scale.factor = 10000) # standard parameters for Seurat
    
      # find variable genes (i.e. features)
      x <- FindVariableFeatures(object = x, selection.method = "vst", nfeatures = 2000)
      
      
      # further normalisation
      if (use.SCTransform) {
        x <- SCTransform(object = x, vars.to.regress = c("nCount_RNA", "percent.mito")) 
      }
      
      # output
      return(x)
    } 
    
  ## Function for brief UMAP visualisation. x must be the output of preprocess.fun(x)
  UMAP.fun <- function(x, dims = 30) {
    x <- ScaleData(object = x, vars.to.regress = c("nCount_RNA", "percent.mito")) 
    x <- RunPCA(x, npcs = dims)
    x <- FindNeighbors(object = x, dims = 1:dims) 
    x <- FindClusters(object = x, resolution = 1) # "We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets"
    x <- RunUMAP(object = x, dims = 1:dims)
  }
  
  ## Function for integration!
  run.integration <- function(x = obj, SCTransform = use.SCTransform) {
    if (SCTransform) {
      int.features <- SelectIntegrationFeatures(object.list = x, nfeatures = 2000)
      int <- PrepSCTIntegration(object.list = x, anchor.features = int.features)
      int <- FindIntegrationAnchors(object.list = int, normalization.method = "SCT", anchor.features = int.features) 
      int <- IntegrateData(anchorset = int, normalization.method = "SCT")
    } else {
      int <- FindIntegrationAnchors(object.list = x)
      int <- IntegrateData(anchorset = int)
    }
    
    return(int)
  }
  
  ## Function to remove features
  filter.features <- function(x, search = "^MT-|^RPS|^RPL") {
    keep <- grep(pattern = search, rownames(x@assays$RNA@counts), invert = TRUE) # ribosomal and mitochondrial genes  
    x <- subset(x, features = keep)
    return(x)
  }


################################################################################################################################ #
## CA ----
  
## Read in
  dat <- read.csv("Raw/human_MTG_2018-06-14_exon-matrix.csv")
  
  # add gene symbol
  meta <- read.csv("Raw/human_MTG_2018-06-14_genes-rows.csv")
  dat <- dat[,-1] # remove an annotation column
  rownames(dat) <- meta$gene
  
  # create Seurat object
  obj$CA <- CreateSeuratObject(counts = dat,
                                min.cells = round(ncol(dat) / 100),
                                min.features = min.features,
                                project = "HCA")
  
  # augment metadata
  meta <- read.csv("Raw/human_MTG_2018-06-14_samples-columns.csv")
  rownames(meta) <- meta$sample_name
  meta <- meta[colnames(obj$CA),] 
  
  obj$CA$Individual <- meta$donor
  obj$CA$orig.celltype <- meta$cluster
  
## Remove cells with no class
  keep <- which(!(obj$CA$orig.celltype == "no class"))
  obj$CA <- subset(obj$CA, cells = keep)
  
## Preprocess
  max.depth <- get.max.depth(obj$CA)
  obj$CA <- preprocess.fun(obj$CA, max.depth = max.depth)
  
## The below section was commented out, and is now run prior to generating the snme mixtures!
  ## Quick filtering: restrict to cell-types with > min.celltype.n members
    # this includes many exc and inh subtypes, as well as ast and oli outgroups
    # keep <- table(obj$CA$orig.celltype)
    # keep <- names(keep[which(keep > min.celltype.n)])
    # keep <- keep[-grep("no class", keep)]
    # 
    # keep <- which(obj$CA$orig.celltype %in% keep)
    # obj$CA <- subset(obj$CA, cells = keep)
    
## UMAP
  obj$CA <- UMAP.fun(obj$CA)
  
  pdf(file = "QC/HCA Initial UMAP.pdf", height = 8, width = 8)
  DimPlot(object = obj$CA, reduction = "umap", group.by = "orig.celltype", label = TRUE)
  dev.off()

################################################################################################################################ #
## Velmeshev ----

## Load
  dat <- Read10X("Raw/Velmeshev2019/")
  obj$VL <- CreateSeuratObject(counts = dat,
                            min.cells = round(ncol(dat) / 100),
                            min.features = min.features,
                            project = "Velmeshev")
  rm(dat)
  
## Annotate
  meta <- read.table("Raw/Velmeshev2019/meta.txt", sep = "\t", header = TRUE)
  m <- match(colnames(obj$VL), meta$cell)
  
  obj$VL$orig.celltype <- meta$cluster[m]
  obj$VL$Region <- meta$region[m]
  obj$VL$Disorder <- meta$diagnosis[m]
  obj$VL$orig.ident <- paste0("Vel_", substr(meta$cell[m], start = 1, stop = 16))
  obj$VL$Individual <- meta$individual[m]

## Quick filtering: restrict to PFC control samples
  keep <- which(obj$VL$Region == "PFC" & obj$VL$Disorder == "Control") # from ~105K to  ~30K nuclei
  obj$VL <- subset(obj$VL, cells = keep)
  # save(obj, file = "Preprocessed/sn_inCreation.rda")
  
## Preprocess and normalise
  max.depth <- get.max.depth(obj$VL)
  obj$VL <- preprocess.fun(obj$VL, max.depth = max.depth)
  
## UMAP
  obj$VL <- UMAP.fun(obj$VL)

  pdf(file = "QC/Velmeshev Initial UMAP.pdf", height = 8, width = 8)
  DimPlot(object = obj$VL, reduction = "umap", group.by = "orig.celltype", label = TRUE)
  dev.off()


################################################################################################################################ #
## Nagy ----

## Load 
  dat <- Read10X("Raw/Nagy2020/")
  obj$NG <- CreateSeuratObject(counts = dat,
                            min.cells = round(ncol(dat) / 100),
                            min.features = min.features,
                            project = "Nagy")
  rm(dat)

## Augment the metadata
  meta <- strsplit(colnames(obj$NG), "\\.")
  meta <- lapply(meta, function(x) {
    x <- c(x[1], strsplit(x[2], "_")[[1]])
    return(x)
  })
  
  meta <- do.call("rbind", meta)
  meta <- as.data.frame(meta)
  colnames(meta) <- c("Celltype", 
                      "Individual", # inferred from there being 17 CTL and 17 MDD patients in the study, and there are 17 levels that are always Control and 17 always Suicide in $Disorder 
                      "Disorder", # note that suicide is equivalent to MDD in this study
                      "Batch", # likely batch given its correspondence across individuals, and that it has six levels which is the number of reported batches. not used by us, so of no importance
                      "Barcode")
  
  obj$NG$orig.celltype <- meta$Celltype
  obj$NG$orig.ident <- paste0("Nagy_", meta$Barcode)
  obj$NG$Individual <- meta$Individual
  obj$NG$Disorder <- meta$Disorder
  
## Minor cell filtering: remove cells labelled as "Mix", or those from MDD individuals
  keep1 <- which(obj$NG$Disorder == "Control")
  keep2 <- grep("Mix_", obj$NG$orig.celltype, invert = TRUE)
  keep <- intersect(keep1, keep2)
  obj$NG <- subset(obj$NG, cells = keep) 
  
## Preprocess and normalise
  max.depth <- get.max.depth(obj$NG)
  obj$NG <- preprocess.fun(obj$NG, max.depth = max.depth)
 
## Visualise dataset properties
  obj$NG <- UMAP.fun(obj$NG)

  pdf(file = "QC/Nagy Initial UMAP.pdf", height = 8, width = 8)
  DimPlot(object = obj$NG, reduction = "umap", group.by = "orig.celltype", label = TRUE)
  dev.off()
  
  
################################################################################################################################ #
## Tasic ----
  
## Load
  # dat <- read.csv("Raw/Tasic2018_GSE115746_cells_exon_counts.csv")
  load("Raw/Tasic2018_GSE115746_cells_exon_counts.rda") # instead of reading the csv
  
  meta <- read.csv("Raw/Tasic2018_ST10_Full_Metadata.csv")
  rownames(meta) <- meta$sample_name
  
## First pass filtering of cells
  # restrict to cells with metadata
  keep <- intersect(colnames(dat), meta$sample_name)
  dat <- dat[,keep]
  meta <- meta[keep,]
  
  # restrict to ALM (anterior lateral motor cortex), as that's a frontal cortical areas
  keep <- which(meta$brain_region == "ALM")
  dat <- dat[,keep]
  meta <- meta[keep,]
  
  # filter out classless cells or those classified as low quality
  remove <- which(meta$class %in% c("Low Quality", "No Class"))
  dat <- dat[,-remove]
  meta <- meta[-remove,]
  
  
## Convert to human EnsID
  # move gene names to rows
  rownames(dat) <- dat$X
  dat <- dat[,-which(colnames(dat) == "X")]
  
  # filter to homologous genes
  homologues <- read.delim("Raw/HOM_MouseHumanSequence.txt")
  dat <- dat[which(rownames(dat) %in% homologues$Symbol),]
  
  # replace mouse gene symbol with human gene symbol
  m <- match(rownames(dat), homologues$Symbol)
  m <- m + 1
  dat$Symbol <- homologues$Symbol[m]
  
  remove <- which(is.na(dat$Symbol)) # removes NAs
  dat <- dat[-remove,]
  
  remove <- dat$Symbol[which(duplicated(dat$Symbol))] # removes duplicated genes
  dat <- dat[-which(dat$Symbol %in% remove),]
  
  rownames(dat) <- dat$Symbol
  dat <- dat[,-which(colnames(dat) == "Symbol")]
  
## Make into Seurat object!
  obj$TS <- CreateSeuratObject(counts = dat,
                               min.cells = round(ncol(dat) / 100),
                               min.features = min.features,
                               project = "Tasic")
  rm(dat)
  
## Augment the metadata
  obj$TS$Individual <- as.character(meta$donor)
  obj$TS$orig.celltype <- as.character(meta$cluster)
  obj$TS$Class <- as.character(meta$class)
  obj$TS$Subclass <- as.character(meta$subclass)
  
## Remove celltypes
  keep <- grep("Peri|^SMC|^VLMC", obj$TS$orig.celltype, invert = TRUE) # this removes pericytes, smooth muscle cells, and VLMCs
  obj$TS <- subset(obj$TS, cells = keep)  
  
## Preprocess and normalise
  max.depth <- get.max.depth(obj$TS)
  obj$TS <- preprocess.fun(obj$TS, max.depth = max.depth)
  
  ## Visualise dataset properties
  obj$TS <- UMAP.fun(obj$TS)

  
################################################################################################################################ #
## Lake ----

## Load
  dat <- read.table("Raw/Lake2018_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt", sep = "\t")
  
  meta <- data.frame(SampleID = colnames(dat), 
                                 orig.celltype = sapply(strsplit(colnames(dat), "_"), `[`, 1),
                                 Individual = sapply(strsplit(colnames(dat), "_"), `[`, 2))
  
  obj$LK <- CreateSeuratObject(counts = dat,
                               min.cells = round(ncol(dat) / 100),
                               min.features = min.features,
                               project = "Lake")  

## Augment metadata
  obj$LK$orig.celltype <- meta$orig.celltype
  obj$LK$Individual <- meta$Individual
  
## Minor cell filtering: remove Pericytes
  keep <- which(obj$LK$orig.celltype != "Per")
  obj$LK <- subset(obj$LK, cells = keep) 
  
## Preprocess and normalise
  max.depth <- get.max.depth(obj$LK)
  obj$LK <- preprocess.fun(obj$LK, max.depth = max.depth)
  
## Visualise dataset properties
  obj$LK <- UMAP.fun(obj$LK)
  

################################################################################################################################ #
## Save ----
  
## Save
  obj <- lapply(obj, filter.features) # this removes mitochondrial and ribosomal reads
  save(obj, file = "Preprocessed/SeuratObjects.rda")
  
      
################################################################################################################################ #
## Create base mixtures ----
  
  
## First: prefilter to remove any lowly-abundant cell-types from being used in the simulation. 
  # for CA
  keep <- table(obj$CA$orig.celltype)
  keep <- names(keep[which(keep > min.celltype.n)])
  # keep <- keep[-grep("no class", keep)]
  
  keep <- which(obj$CA$orig.celltype %in% keep)
  obj$CA <- subset(obj$CA, cells = keep)
  
  # for VL, no cell-types are below this threshold
  
## Setup lists
  snme <- list()
  cellsForX <- list()

## First, partitition nuclei into two discrete groups
  cellsForX$CA.sig <- sample(colnames(obj$CA), size = ncol(obj$CA) / 2, replace = FALSE)
  cellsForX$VL.sig <- sample(colnames(obj$VL), size = ncol(obj$VL) / 2, replace = FALSE)
  
  cellsForX$CA.mix <- colnames(obj$CA)[-which(colnames(obj$CA) %in% cellsForX$CA.sig)]
  cellsForX$VL.mix <- colnames(obj$VL)[-which(colnames(obj$VL) %in% cellsForX$VL.sig)]
  
## Create randomly-sampled mixtures
  # hca
  dat <- as.data.frame(obj$CA@assays$RNA@counts)
  dat <- dat[,cellsForX$CA.mix]
  snme$CA <- create.snmeRand(x = dat,
                              nMix = nMix,
                              nReps = nReps,
                              meta = obj$CA@meta.data[cellsForX$CA.mix,],
                              meta.columns = "orig.celltype",
                              rpkm = FALSE)

  # velmeshev
  dat <- as.data.frame(obj$VL@assays$RNA@counts)
  dat <- dat[,cellsForX$VL.mix]
  snme$VL <- create.snmeRand(x = dat,
                              nMix = nMix,
                              nReps = nReps,
                              meta = obj$VL@meta.data[cellsForX$VL.mix,],
                              meta.columns = "orig.celltype",
                              rpkm = FALSE)

  
  
## Save
  write.CIB(data = snme$CA$mixture$cpm, dir = "Preprocessed/CIB/snme_HCA.txt")
  write.CIB(data = snme$VL$mixture$cpm, dir = "Preprocessed/CIB/snme_Vel.txt")
  
  save(snme, file = "Preprocessed/snme.rda")  
 
  
################################################################################################################################ #
## Create gradient mixtures for use with complete deconvolution ----
  
## Remove mitochondrial and ribosomal reads
  obj <- lapply(obj, filter.features)
  
## From the Vel data
  ## Extract Vel data
    cells <- as.data.frame(obj$VL@assays$RNA@counts)
    meta <- obj$VL@meta.data
    colnames(cells) <- meta$orig.celltype
    
  ## Note: no need to partition nuclei into two discrete groups, as we're not using the VL signature!!
    
  ## Remove NRGN
    remove <- grep("NRGN", meta$orig.celltype)
    cells <- cells[,-remove]
    meta <- meta[-remove,]
    
  ## Gradient simulation 1: 
  nMix <- 500
  ct <- names(table(meta$grad1))
    
    # here, you want to restrict to just the merge2 cell-types, and change exc proportion greatly
    meta$grad1 <- meta$orig.celltype
    meta$grad1[grep("AST", meta$grad1)] <- "Astrocytes"
    meta$grad1[grep("OPC", meta$grad1)] <- "OPCs"
    meta$grad1[grep("^IN", meta$grad1)] <- "Neurons_Inh"
    meta$grad1[grep("^L|mat", meta$grad1)] <- "Neurons_Exc"
    
    ## Set proportions of excitatory cells
    props <- list()
    
      ## A) excitatory range from 0-50%
      j <- list(meta$grad1)
      j[[1]][-grep("Exc", j[[1]])] <- "Other"
      norm.abundance <- aggregate(meta$nCount_RNA, j, mean)
      rownames(norm.abundance) <- norm.abundance$Group.1
      
      exc.scale.factor <- norm.abundance["Other", "x"] / norm.abundance["Neurons_Exc", "x"] 
    
      props$A <- ceiling(data.frame(Neurons_Exc = c(0, 0.5*exc.scale.factor*nMix, 0, 0.5*exc.scale.factor*nMix)))
      
      ## B): balanced abundance of ct
        norm.abundance <- aggregate(meta$nCount_RNA, list(meta$grad1), mean)
        # norm.abundance$x <- max(norm.abundance$x) / norm.abundance$x
        norm.abundance$OverMax <- max(norm.abundance$x) / norm.abundance$x
        norm.abundance$OverMin <- norm.abundance$x / min(norm.abundance$x)
        
        base.n <- (nMix*3) / sum(norm.abundance$OverMax)
        
        norm.n <- floor(norm.abundance$OverMax * base.n)
        norm.n <- lapply(norm.n, function(y) { data.frame(c(0, y, 0, y)) })
        norm.n <- do.call("cbind", norm.n)
        colnames(norm.n) <- norm.abundance$Group.1
        props$B <- norm.n
        

    # change cellnames for this run
    colnames(cells) <- meta$grad1
    
    # run
    snme.gradients <- list()
    for (j in names(props)) {
      snme.gradients[[j]] <- confound.proportion(dat = cells, 
                                                 ct = ct, 
                                                 setProps = props[[j]], 
                                                 nPerGroup = 50, 
                                                 nPerMixture = nMix,
                                                 keepPureExpression = FALSE)
    }
    
 
    # cpm normalise
    snme.gradients <- lapply(snme.gradients, function(x) {
      x$confound.p <- apply(x$confound.p, 2, function(y) {
        lib.size <- 10^6 / sum(y)
        y <- y * lib.size    
      })
      return(x)
    })
    
  ## Save
  save(snme.gradients, file = "Preprocessed/snme_grad_VL.rda")  
  

## From the CA data
  ## Extract CA data
    cells <- as.data.frame(obj$CA@assays$RNA@counts)
    meta <- obj$CA@meta.data
    colnames(cells) <- meta$orig.celltype
  
  ## Note: no need to partition nuclei into two discrete groups, as we're not using the CA signature!!
  
  ## Filter CA cells as per other simulations
    keep <- names(which(table(meta$orig.celltype) > 200))
    keep <- which(meta$orig.celltype %in% keep)
    
    cells <- cells[,keep]
    meta <- meta[keep,]
    meta$orig.celltype <- as.character(meta$orig.celltype)
    
  ## Gradient simulation
    nMix <- 500
    
  
  ## Merge celltypes
    meta$grad1 <- meta$orig.celltype
    meta$grad1[grep("Ast", meta$grad1)] <- "Astrocytes"
    meta$grad1[grep("OPC", meta$grad1)] <- "OPCs"
    meta$grad1[grep("Oli", meta$grad1)] <- "Oligodendrocytes"
    meta$grad1[grep("Inh", meta$grad1)] <- "Neurons_Inh"
    meta$grad1[grep("Exc", meta$grad1)] <- "Neurons_Exc"
    ct <- names(table(meta$grad1))
    
  ## Set proportions of excitatory cells
    props <- list()
    
    ## A) excitatory range from 0-50%
      j <- list(meta$grad1)
      j[[1]][-grep("Exc", j[[1]])] <- "Other"
      norm.abundance <- aggregate(meta$nCount_RNA, j, mean)
      rownames(norm.abundance) <- norm.abundance$Group.1
      
      exc.scale.factor <- norm.abundance["Other", "x"] / norm.abundance["Neurons_Exc", "x"] 
      
      props$A <- ceiling(data.frame(Neurons_Exc = c(0, 0.5*exc.scale.factor*nMix, 0, 0.5*exc.scale.factor*nMix)))
      
    ## B): balanced abundance of ct
      norm.abundance <- aggregate(meta$nCount_RNA, list(meta$grad1), mean)
      # norm.abundance$x <- max(norm.abundance$x) / norm.abundance$x
      norm.abundance$OverMax <- max(norm.abundance$x) / norm.abundance$x
      norm.abundance$OverMin <- norm.abundance$x / min(norm.abundance$x)
      
      base.n <- (nMix*2) / sum(norm.abundance$OverMax)
      
      norm.n <- floor(norm.abundance$OverMax * base.n)
      norm.n <- lapply(norm.n, function(y) { data.frame(c(0, y, 0, y)) })
      norm.n <- do.call("cbind", norm.n)
      colnames(norm.n) <- norm.abundance$Group.1
      props$B <- norm.n
      
  # change cellnames for this run
  colnames(cells) <- meta$grad1
  
  # run
  snme.gradients <- list()
  for (j in names(props)) {
    snme.gradients[[j]] <- confound.proportion(dat = cells, 
                                               ct = ct, 
                                               setProps = props[[j]], 
                                               nPerGroup = 50, 
                                               nPerMixture = nMix,
                                               keepPureExpression = FALSE)
  }
  
  
  # cpm normalise
  snme.gradients <- lapply(snme.gradients, function(x) {
    x$confound.p <- apply(x$confound.p, 2, function(y) {
      lib.size <- 10^6 / sum(y)
      y <- y * lib.size    
    })
    return(x)
  })
  
  ## Save
  save(snme.gradients, file = "Preprocessed/snme_grad_CA.rda")  
  

################################################################################################################################ #
## Create signatures from the above Seurat objects ----
  
sigsSNME <- list()
merges <- list()
drops <- list()


## In this section, special versions of the CA and VL signatures are made, for deconvolving their respective mixtures


## First, get base signature with all cell-types
  for(j in names(obj)) {
    # setup
    print(paste0("Starting signature ", j))
    x <- list()
    use <- cellsForX[[paste0(j, ".sig")]]

    # collect metadata, for use with MuSiC
    x$Meta <- obj[[j]]@meta.data[use,]

    # collect all cells' count, for use with MuSiC
    x$Full.counts <- as.data.frame(obj[[j]]@assays$RNA@counts[,use])

    # collect average cpm expression of every cell
    x$Full.cpm <- apply(x$Full.counts, 2, function(x) {
      lib.size <- 10^6 / sum(x)
      x <- x * lib.size
      return(x)
    })

    # collect average cpm of original clusters
    y <- list()
    for(k in levels(as.factor(x$Meta$orig.celltype))) {
      ct <- which(x$Meta$orig.celltype == k)
      y[[k]] <- rowMeans(x$Full.cpm[,ct])
    }
    x$orig.celltype <- as.data.frame(do.call("cbind", y))

    # remove $Full.cpm for filesize reason
    x <- x[-grep("Full.cpm", names(x))]

    sigsSNME[[j]] <- x
  }


## Next, create signatures where related cell-types are pooled together
  merges <- list()

  ## Function to merge
    merger <- function(x = sigsSNME, sig, levels, name) {
      y <- list() # to hold expression
      z <- x[[sig]]$Meta # to hold labels
      z[,name] <- NA

      for (j in levels) {
        print(paste0("Merging ", j, " of ", sig, "_", name))

        # collect indices
        ct <- merges[[sig]][[j]]
        representatives <- which(x[[sig]]$Meta$orig.celltype %in% ct)
        z[representatives,name] <- j

        # collect expression
        y[[j]] <- rowMeans(x[[sig]]$Full.cpm[,representatives])
      }

      # combine expression into a single dataframe
      y <- as.data.frame(do.call("cbind", y))

      # edit x
      x[[sig]][[name]] <- y
      x[[sig]]$Meta <- z
      return(x)
    }

  ## Merge celltypes in Velmeshev
    merges$VL <- list()

    # main cell-types
    merges$VL$Neurons <- c("IN-PV", "IN-SST", "IN-SV2C", "IN-VIP", "L2/3", "L4", "L5/6", "L5/6-CC", "Neu-mat", "Neu-NRGN-I", "Neu-NRGN-II")
    merges$VL$Astrocytes <- c("AST-FB", "AST-PP")
    merges$VL$Endothelial <- "Endothelial"
    merges$VL$Microglia <- "Microglia"
    merges$VL$Oligodendrocytes <- "Oligodendrocytes"
    merges$VL$OPC <- "OPC"

    # sub-types of neurons
    merges$VL$Neurons_Inh <- c("IN-PV", "IN-SST", "IN-SV2C", "IN-VIP")
    merges$VL$Neurons_Exc <- c("L2/3", "L4", "L5/6", "L5/6-CC", "Neu-mat")
    merges$VL$Neurons_NRGN <- c("Neu-NRGN-I", "Neu-NRGN-II")

    # finer resolution
    # merges$VL$Neurons_InhA <- c("IN-SST", "IN-VIP")
    # merges$VL$Neurons_InhB <- c("IN-PV","IN-SV2C")
    # merges$VL$Neurons_ExcA <- c("L2/3", "L4", "L5/6-CC", "Neu-mat")
    # merges$VL$Neurons_ExcB <- c("L5/6")

    merges$VL$Neurons_InhA <- c("IN-SST", "IN-PV")
    merges$VL$Neurons_InhB <- c("IN-VIP","IN-SV2C")
    merges$VL$Neurons_ExcA <- c("L2/3", "L4", "L5/6-CC")
    merges$VL$Neurons_ExcB <- c("L5/6")
    merges$VL$Neurons_ExcC <- c("Neu-mat")

    # apply
    sigsSNME <- merger(sig = "Vel",
                        levels = c("Neurons",
                                   "Astrocytes",
                                   "Oligodendrocytes",
                                   "Microglia",
                                   "Endothelial",
                                   "OPC"),
                        name = "merge1")
    sigsSNME <- merger(sig = "Vel",
                        levels = c("Neurons_Exc",
                                   "Neurons_Inh",
                                   "Neurons_NRGN",
                                   "Astrocytes",
                                   "Oligodendrocytes",
                                   "Microglia",
                                   "Endothelial",
                                   "OPC"),
                        name = "merge2")
    sigsSNME <- merger(sig = "Vel",
                        levels = c("Neurons_ExcA",
                                   "Neurons_InhA",
                                   "Neurons_ExcB",
                                   "Neurons_InhB",
                                   "Neurons_ExcC",
                                   "Neurons_NRGN",
                                   "Astrocytes",
                                   "Oligodendrocytes",
                                   "Microglia",
                                   "Endothelial",
                                   "OPC"),
                        name = "merge3")

  ## HCA
    x <- levels(as.factor(sigsSNME$CA$Meta$orig.celltype))
    merges$CA <- list()

    # merge to the level of cell-types
    merges$CA$Neurons <- x[grep("^Ex|^Inh", x)]
    merges$CA$Astrocytes <- x[grep("^Astro", x)]
    merges$CA$Oligodendrocytes <- x[grep("^Oligo", x)]
    merges$CA$OPCs <- x[grep("^OPC", x)]

    # sub-types of neurons
    merges$CA$Neurons_Exc <- x[grep("^Ex", x)]
    merges$CA$Neurons_Inh <- x[grep("^Inh", x)]

    # apply
    sigsSNME <- merger(sig = "HCA",
                        levels = c("Neurons",
                                   "Astrocytes",
                                   "Oligodendrocytes",
                                   "OPCs"),
                        name = "merge1")
    sigsSNME <- merger(sig = "HCA",
                        levels = c("Neurons_Exc",
                                   "Neurons_Inh",
                                   "Astrocytes",
                                   "Oligodendrocytes",
                                   "OPCs"),
                        name = "merge2")


## Signatures in which cell-types are dropped
  ## This is for Velmeshev only...
  ## Drop1: remove all inh neurons, and deconvolve using merged celltypes
    # for DRS/DTA/CIB
    g <- grep("Inh", colnames(sigsSNME$VL$merge2))
    sigsSNME$VL$drop1 <- sigsSNME$VL$merge2[,-g]

    # for Music
    x <- sigsSNME$VL$Meta$merge2
    x[grep("Inh", x)] <- "Drop"
    sigsSNME$VL$Meta$drop1 <- x

  ## Drop2: remove one (abundant) subtype of excitatory neurons, and remerge
    # for DRS/DTA/CIB
    sigsSNME$VL$drop2 <- sigsSNME$VL$merge2

    g <- grep("L4", merges$VL$Neurons_Exc)
    g <- merges$VL$Neurons_Exc[-g]

    sigsSNME$VL$drop2$Neurons_Exc <- rowMeans(sigsSNME$VL$Full.cpm[,which(sigsSNME$VL$Meta$orig.celltype %in% g)]) # cor between left and right is 0.9996919

    # for Music
    x <- sigsSNME$VL$Meta$merge2
    x[grep("L4", sigsSNME$VL$Meta$orig.celltype)] <- "Drop"
    sigsSNME$VL$Meta$drop2 <- x

  ## Drop3: remove one (abundant) subtype of excitatry neurons, but don't remerge
    # for DRS/DTA/CIB
    g <- grep("L4", colnames(sigsSNME$VL$orig.celltype))
    sigsSNME$VL$drop3 <- sigsSNME$VL$orig.celltype[,-g] # removes L4, typically 15% of the mixtures

    # for Music
    x <- sigsSNME$VL$Meta$orig.celltype
    x[grep("L4", x)] <- "Drop"
    sigsSNME$VL$Meta$drop3 <- x

  ## Drop4: remove Oli, deconvolve merged
    # for DRS/DTA/CIB
    g <- grep("Oli", colnames(sigsSNME$VL$merge1))
    sigsSNME$VL$drop4 <- sigsSNME$VL$merge1[,-g]

    # for Music
    x <- sigsSNME$VL$Meta$merge1
    x[grep("Oli", x)] <- "Drop"
    sigsSNME$VL$Meta$drop4 <- x

## Save
    save(sigsSNME, file = "Preprocessed/Signatures - SNME (incomplete).rda")

    
################################################################################################################################ #
## Correlations between and within signatures ---- 
  
## Function
  plot.heatmap <- function(x, y = NA, method = "Spearman", main = NA, combine = FALSE, xName = NA, yName = NA, abundances = NA, h.margin = 10, v.margin = 10) {
    # setup
    colours <- rev(colorRampPalette(brewer.pal(8, "PiYG"))(25))
    met <- "p"; if(method == "Spearman") met <- "s"
    
    
    # create correlation matrix
    if (is.na(y)) {
      corMx <- cor(x, method = met)
      if(!(is.na(abundances))) {
        m <- match(rownames(corMx), names(abundances))
        rownames(corMx) <- paste0(rownames(corMx), " (", abundances[m], ")")
      }
      
    } else {
      if(combine) {
        colnames(x) <- paste0(xName, colnames(x))
        colnames(y) <- paste0(yName, colnames(y))
        
        dat <- cbind(x, y)
        sideColours <- c(rep("black", ncol(x)),
                         rep("grey90", ncol(y)))
        corMx <- cor(dat, method = met)
      } else {
        corMx <- cor(y, x, method = met)  
      }
    }
    
    # heatmap
    if (combine) {
      heatmap.2(corMx, trace = "none", col = colours, density.info = "none", key.title = "NA", key.xlab = method, main = main, ColSideColors = sideColours, RowSideColors = sideColours)  
    } else {
      heatmap.2(corMx, trace = "none", col = colours, density.info = "none", key.title = "NA", key.xlab = method, main = main, margins = c(v.margin,h.margin))
    }
    
  }

## Produce correlation heatmaps within a signature
  pdf(file = "../Results/Revisions/snme/Signature Heatmaps - All Celltypes.pdf", height = 8, width = 8)
  plot.heatmap(sigsSNME$VL$orig.celltype, main = "VL", abundances = table(obj$VL$orig.celltype), h.margin = 10.5, v.margin = 9)  
  # plot.heatmap(sigsSNME$NG$orig.celltype, main = "NG", abundances = table(obj$NG$orig.celltype), h.margin = 8, v.margin = 7)  
  plot.heatmap(sigsSNME$CA$orig.celltype, main = "CA", abundances = table(obj$CA$orig.celltype), h.margin = 15, v.margin = 13)  
  dev.off()
  
  pdf(file = "../Results/Revisions/snme/Signature Heatmaps - All Celltypes (Pearson).pdf", height = 6, width = 8)
  plot.heatmap(sigsSNME$VL$orig.celltype, main = "Vel", method = "Pearson")  
  plot.heatmap(sigsSNME$NG$orig.celltype, main = "Nagy", method = "Pearson")  
  plot.heatmap(sigsSNME$CA$orig.celltype, main = "HCA", method = "Pearson")  
  dev.off()
  
  pdf(file = "../Results/Revisions/snme/Signature Heatmaps - All Celltypes (Pearson) (Core Genes).pdf", height = 6, width = 8)
  plot.heatmap(sigsSNME$VL$orig.celltype[core.genes,], main = "Vel", method = "Pearson")  
  plot.heatmap(sigsSNME$NG$orig.celltype[core.genes,], main = "Nagy", method = "Pearson")  
  plot.heatmap(sigsSNME$CA$orig.celltype[core.genes,], main = "HCA", method = "Pearson")  
  dev.off()
  
  pdf(file = "../Results/Revisions/snme/Signature Heatmaps - All Celltypes (Core Genes).pdf", height = 6, width = 8)
  plot.heatmap(sigsSNME$VL$orig.celltype[core.genes,], main = "Vel")
  plot.heatmap(sigsSNME$NG$orig.celltype[core.genes,], main = "Nagy")
  plot.heatmap(sigsSNME$CA$orig.celltype[core.genes,], main = "HCA")
  dev.off()
  
# ## Produce correlation heatmaps across signatures
#   # first, find common genes
#   core.genes <- sapply(sigsSNME, function(x) rownames(x$orig.celltype))
#   core.genes <- do.call("c", core.genes)
#   core.genes <- (table((core.genes)) == 3)
#   core.genes <- names(core.genes)[which(core.genes)]
#   
#   ## Using all cell-types in a signature
#     pdf(file = "../Results/Revisions/snme/Signature Heatmaps - All Celltypes Across Signatures.pdf", height = 6, width = 8)
#     # vel vs. nagy
#     plot.heatmap(x = sigsSNME$VL$orig.celltype[core.genes,], 
#                  y = sigsSNME$NG$orig.celltype[core.genes,],
#                  main = "Vel (x) vs. Nagy (y)")    
#     
#     # vel vs. hca
#     plot.heatmap(x = sigsSNME$VL$orig.celltype[core.genes,], 
#                  y = sigsSNME$CA$orig.celltype[core.genes,],
#                  main = "Vel (x) vs. HCA (y)")    
#     
#     # nagy vs. hca
#     plot.heatmap(x = sigsSNME$NG$orig.celltype[core.genes,], 
#                  y = sigsSNME$CA$orig.celltype[core.genes,],
#                  main = "Nagy (x) vs. HCA (y)")    
#     dev.off()
#     
#   ## Of major cell classes
#     pdf(file = "../Results/Revisions/snme/Signature Heatmaps - merge1 Across Signatures.pdf", height = 6, width = 8)
#     # vel vs. nagy
#     plot.heatmap(x = sigsSNME$VL$merge1[core.genes,], 
#                  y = sigsSNME$NG$merge1[core.genes,],
#                  main = "Vel (x) vs. Nagy (y)")    
#     
#     # vel vs. hca
#     plot.heatmap(x = sigsSNME$VL$merge1[core.genes,], 
#                  y = sigsSNME$CA$merge1[core.genes,],
#                  main = "Vel (x) vs. HCA (y)")    
#     
#     # nagy vs. hca
#     plot.heatmap(x = sigsSNME$NG$merge1[core.genes,], 
#                  y = sigsSNME$CA$merge1[core.genes,],
#                  main = "Nagy (x) vs. HCA (y)")    
#     dev.off()
#      
#   ## Of neuronal subtypes
#     pdf(file = "../Results/Revisions/snme/Signature Heatmaps - Neuronal Subtypes Across Signatures.pdf", height = 6, width = 8)
#     # vel vs. nagy
#     plot.heatmap(x = sigsSNME$VL$orig.celltype[core.genes, merges$VL$Neurons], 
#                  y = sigsSNME$NG$orig.celltype[core.genes, merges$NG$Neurons,],
#                  main = "Vel (x) vs. Nagy (y)")    
#     
#     # vel vs. hca
#     plot.heatmap(x = sigsSNME$VL$orig.celltype[core.genes,merges$VL$Neurons], 
#                  y = sigsSNME$CA$orig.celltype[core.genes,merges$CA$Neurons],
#                  main = "Vel (x) vs. HCA (y)")    
#     
#     # nagy vs. hca
#     plot.heatmap(x = sigsSNME$NG$orig.celltype[core.genes,merges$NG$Neurons], 
#                  y = sigsSNME$CA$orig.celltype[core.genes,merges$CA$Neurons],
#                  main = "Nagy (x) vs. HCA (y)")    
#     dev.off()
#     
# ## Produce correlation heatmaps across signatures, but clustering across samples!
#   ## Of neuronal subtypes
#     pdf(file = "../Results/Revisions/snme/Signature Heatmaps - Neuronal Subtypes Across Signatures, Clustered.pdf", height = 6, width = 8)
#     # vel vs. nagy
#     plot.heatmap(x = sigsSNME$VL$orig.celltype[core.genes, merges$VL$Neurons], 
#                  y = sigsSNME$NG$orig.celltype[core.genes, merges$NG$Neurons,],
#                  main = "Vel vs. Nagy",
#                  combine = TRUE, xName = "V_", yName = "N_")    
#     
#     # vel vs. hca
#     plot.heatmap(x = sigsSNME$VL$orig.celltype[core.genes,merges$VL$Neurons], 
#                  y = sigsSNME$CA$orig.celltype[core.genes,merges$CA$Neurons],
#                  main = "Vel vs. HCA",
#                  combine = TRUE, xName = "V_", yName = "H_")    
#     
#     # nagy vs. hca
#     plot.heatmap(x = sigsSNME$NG$orig.celltype[core.genes,merges$NG$Neurons], 
#                  y = sigsSNME$CA$orig.celltype[core.genes,merges$CA$Neurons],
#                  main = "Nagy vs. HCA",
#                  combine = TRUE, xName = "N_", yName = "H_")    
#     dev.off()
    
    
    
################################################################################################################################ #
## Analyse properties of the simulated mixtures ---- 
    
# here, we want to explore how well the simulations approximate a bulk sample
    
## Collect data
  ## Prepare
    test.data <- list()
    test.data$counts <- list()
    test.data$CPM <- list()
    
  ## Get bulk brain data from Parikshak et al.
    load("Preprocessed/ASD.rda")
    pCPM <- apply(pSymbol, 2, function(x) {
      lib.size <- 10^6 / sum(x)
      x <- x * lib.size  
      return(x)  
    })
    
  ## Get SC data from Darmanis et al.
    # the scme mixtures
    load("Preprocessed/scme.rda") 
    scme <- addSymbol(scme)
    
    # the single cells
    DM.rpkm <- read.csv("Raw/Darmanis2015_scRNAseq.csv")  # load and process single cells
    rownames(DM.rpkm) <- DM.rpkm[,1]
    DM.rpkm <- DM.rpkm[,-1]
    
    DM.meta <- read.csv("Raw/Darmanis2015_Meta.txt", sep = "\t") # filtering cells to only those in the signature
    DM.meta <- DM.meta[-grep("hybrid", DM.meta$cell_type_s),]
    use <- which(DM.meta$cell_type_s %in% c("neurons", "astrocytes", "OPC", "microglia", "endothelia", "oligodendrocytes"))
    DM.rpkm <- DM.rpkm[,use]
    
    
    DM.rpkm <- apply(DM.rpkm, 2, function(x) {
      lib.size <- 10^6 / sum(x)
      x <- x * lib.size  
      return(x)  
    })
    
    
    
  ## Get CPM on nuclei
    VL.cpm <- apply(sigsSNME$VL$Full.counts, 2, function(x) {
      lib.size <- 10^6 / sum(x)
      x <- x * lib.size  
      return(x)  
    })
    CA.cpm <- apply(sigsSNME$CA$Full.counts, 2, function(x) {
      lib.size <- 10^6 / sum(x)
      x <- x * lib.size  
      return(x)  
    })
    
    n.sample <- 10
      
  # of 5 single-cells (nuclei, rather)
  use.sn <- sample(rownames(sigsSNME$VL$Meta), size = n.sample, replace = FALSE) # choose 5 samples
  test.data$counts$SN.Vel <- sigsSNME$VL$Full.counts[,use.sn]
  test.data$CPM$SN.Vel <- VL.cpm[,use.sn]
  
  use.sn <- sample(rownames(sigsSNME$CA$Meta), size = n.sample, replace = FALSE) # choose 5 samples
  test.data$counts$SN.HCA <- sigsSNME$CA$Full.counts[,use.sn]
  test.data$CPM$SN.HCA <- CA.cpm[,use.sn]
      
  # of 5 simulated Vel mixtures
  use.mix <- sample(colnames(snme$VL$mixtur$counts), size = n.sample, replace = FALSE)
  test.data$counts$Mix.Vel <- snme$VL$mixture$counts[,use.mix]
  test.data$CPM$Mix.Vel <- snme$VL$mixture$cpm[,use.mix]
  
  use.mix <- sample(colnames(snme$CA$mixture$counts), size = n.sample, replace = FALSE)
  test.data$counts$Mix.HCA <- snme$CA$mixture$counts[,use.mix]
  test.data$CPM$Mix.HCA <- snme$CA$mixture$cpm[,use.mix]
      
  # of 5 brain samples from Parikshak et al
  use.brain <- sample(colnames(pSymbol), size = n.sample, replace = FALSE)
  test.data$counts$Brain <- pSymbol[,use.brain]
  test.data$CPM$Brain <- pCPM[,use.brain]
  
  # of 5 mixtures and single cells from Darmanis et al.
  use.brain <- sample(colnames(DM.rpkm), size = n.sample, replace = FALSE)
  test.data$CPM$SC.DM <- DM.rpkm[,use.brain]
  
  use.brain <- sample(colnames(scme), size = n.sample, replace = FALSE)
  test.data$CPM$Mix.DM <- scme[,use.brain]
    
  # filter to common genes
  core.genes <- sapply(test.data$CPM, rownames)
  core.genes <- do.call("c", core.genes)
  core.genes <- (table((core.genes)) == 7)
  core.genes <- names(core.genes)[which(core.genes)]
  
  test.data <- lapply(test.data, function(x) {
    lapply(x, function(y) y[core.genes,])
  })
  
  # set column names
  for (j in names(test.data)) {
    for (k in names(test.data[[j]])) {
      colnames(test.data[[j]][[k]]) <- LETTERS[1:10]
    }
  }
      
## Properties on CPM
  plot.data <- melt(test.data$CPM)
  m <- which(is.na(plot.data$variable))
  plot.data$variable[m] <- plot.data$Var2[m]
  plot.data$L1 <- factor(plot.data$L1)
  levels(plot.data$L1) <- c("Bulk Brain", "DM Simulation", "CA Simulation", "VL Simulation", "DM Cells", "CA Nuclei", "VL Nuclei")
  
  pdf(file = "QC/SNME - Expression Distributions (CPM-level).pdf", height = 2.5, width = 8)
  ggplot(plot.data, aes(x = log10(value + 1), colour = variable)) +
    geom_density(position = "dodge") +
    facet_wrap(~L1, scale = "free_y", nrow = 1) +
    labs(y = "Density", x = "Log10(CPM + 1)") +
    # scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
    scale_colour_carto_d(palette = "Prism") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"), legend.position = "none", axis.ticks.y = element_blank())
  
  
  # proportion of zeroes
  plot.data$Zero <- plot.data$value == 0
  plot.data$Zero <- as.numeric((plot.data$Zero)) # sets TRUE to 1 and FALSE to 0
  plot.data$Zero <- plot.data$Zero / length(core.genes) # on the plot, this makes it a percentage!
  
  ggplot(plot.data, aes(x = variable, y = Zero, fill = variable)) +
    geom_col() +
    facet_wrap(~L1, nrow = 1) +
    theme_bw() +
    labs(y = "Percentage of Zeroes", x = "Randomly-chosen Sample") +
    scale_fill_carto_d(palette = "Prism") +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
    theme(legend.position = "None")
  dev.off()  

