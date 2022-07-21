## In this script, mixtures are made from single-nucleus data

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
  source("../Scripts/Fun_Composition.R")
  load(file = "Preprocessed/geneInfo.rda")
  load("Preprocessed/exonicLength.rda")

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
    
  # ## Function for brief UMAP visualisation. x must be the output of preprocess.fun(x)
  # UMAP.fun <- function(x, dims = 30) {
  #   x <- ScaleData(object = x, vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #   x <- RunPCA(x, npcs = dims)
  #   x <- FindNeighbors(object = x, dims = 1:dims) 
  #   x <- FindClusters(object = x, resolution = 1) # "We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets"
  #   x <- RunUMAP(object = x, dims = 1:dims)
  # }
  
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
  
# Download this and unzip: https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq  
 

## Read in
  dat <- read.csv("Raw/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_exon-matrix.csv") 
  
  # add gene symbol
  meta <- read.csv("Raw/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_genes-rows.csv")
  dat <- dat[,-1] # remove an annotation column
  rownames(dat) <- meta$gene
  
  # create Seurat object
  obj$CA <- CreateSeuratObject(counts = dat,
                                min.cells = round(ncol(dat) / 100),
                                min.features = min.features,
                                project = "HCA")
  
  # augment metadata
  meta <- read.csv("Raw/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_samples-columns.csv")
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
  
## Remove all (9) endothelial cells
  obj$CA <- subset(obj$CA, cells = grep("Endo", obj$CA$orig.celltype, invert = TRUE))
  
## The below section was commented out, and is now run prior to generating the snme mixtures!
  ## Quick filtering: restrict to cell-types with > min.celltype.n members
    # this includes many exc and inh subtypes, as well as ast and oli outgroups
    # keep <- table(obj$CA$orig.celltype)
    # keep <- names(keep[which(keep > min.celltype.n)])
    # keep <- keep[-grep("no class", keep)]
    # 
    # keep <- which(obj$CA$orig.celltype %in% keep)
    # obj$CA <- subset(obj$CA, cells = keep)
    

################################################################################################################################ #
## Velmeshev ----

## Load
  dat <- Read10X("Raw/Velmeshev2019/rawMatrix/")
  obj$VL <- CreateSeuratObject(counts = dat,
                            min.cells = round(ncol(dat) / 100),
                            min.features = min.features,
                            project = "Velmeshev")
  rm(dat)
  
## Annotate
  meta <- read.table("Raw/Velmeshev2019/meta.tsv", sep = "\t", header = TRUE)
  m <- match(colnames(obj$VL), meta$cell)
  
  obj$VL$orig.celltype <- meta$cluster[m]
  obj$VL$Region <- meta$region[m]
  obj$VL$Disorder <- meta$diagnosis[m]
  obj$VL$orig.ident <- paste0("Vel_", substr(meta$cell[m], start = 1, stop = 16))
  obj$VL$Individual <- meta$individual[m]

## Quick filtering: restrict to PFC control samples
  keep <- which(obj$VL$Region == "PFC" & obj$VL$Disorder == "Control") # from ~105K to  ~30K nuclei
  obj$VL <- subset(obj$VL, cells = keep)
  
## Preprocess and normalise
  max.depth <- get.max.depth(obj$VL)
  obj$VL <- preprocess.fun(obj$VL, max.depth = max.depth)
  


################################################################################################################################ #
## Nagy ----

## Download these three files from GEO:
  # GSE144136_GeneNames.csv.gz	
  # GSE144136_CellNames.csv.gz	
  # GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz	
  
  # Then put in the directory Data/Raw/Nagy2020, unzip, and rename as: features.tsv, barcodes.tsv, and matrix.mtx, respectively.
  # Note that this will require converting features.csv and barcodes.csv to tsvs
  
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
 
  
################################################################################################################################ #
## Tasic ----
  
# download the file 
# and access Supplementary Table 10 for metadata: https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0654-5/MediaObjects/41586_2018_654_MOESM3_ESM.zip
  
## Load
  dat <- read.csv("Raw/GSE115746_cells_exon_counts.csv")
  # load("Raw/Tasic2018_GSE115746_cells_exon_counts.rda") # instead of reading the csv
  
  # meta <- read.csv("Raw/Tasic2018_ST10_Full_Metadata.csv")
  meta <- readxl::read_xlsx("Raw/41586_2018_654_MOESM3_ESM/Supplementary_Table_10_Full_Metadata.xlsx", sheet = 1)
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$sample_name
  
## Convert to human EnsID
  # move gene names to rows
  rownames(dat) <- dat$X
  
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
  
  
################################################################################################################################ #
## Lake ----

## Load
  dat <- read.table("Raw/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt", sep = "\t")
  
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
  
  
## Preprocess and normalise
  max.depth <- get.max.depth(obj$LK)
  obj$LK <- preprocess.fun(obj$LK, max.depth = max.depth)
  
## Minor cell filtering: remove Pericytes and Endothelia
  keep <- which(!(obj$LK$orig.celltype %in% c("Per", "End")))
  obj$LK <- subset(obj$LK, cells = keep) 

################################################################################################################################ #
## Save ----
  
## Save
  obj <- lapply(obj, filter.features) # this removes mitochondrial and ribosomal reads
  save(obj, file = "Preprocessed/SeuratObjects.rda")
  
      


    
    
