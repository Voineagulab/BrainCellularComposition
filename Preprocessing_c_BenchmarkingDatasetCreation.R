## In this script, the simulated mixtures using CA, VL, and DM single-cells/nuclei are generated

################################################################################################################################ #
## Setup ----

## Start!
rm(list = ls()); gc()

## Directory
# root.dir <- "INSERT.YOUR.ROOT.DIRECTORY.PLEASE/"
root.dir <- "/Volumes/Data1/PROJECTS/BrainCellularComposition/BrainCellularComposition_GitHub/"
wd1 <- paste0(root.dir, "Data/") # working directory 1, for processing
setwd(wd1)

## Functions and packages
source("../Scripts/Fun_Composition.R")
source("../Scripts/Fun_Preprocessing.R")
source("../Scripts/Fun_Parameters.R")

## Annotation files
  load(file = "Preprocessed/geneInfo.rda")
  load("Preprocessed/exonicLength.rda")


################################################################################################################################ #
## Create a benchmarking dataset using the DM single-cells  ----    

## Preprocess data from Darmanis et al.
  # read in
  darmanis <- read.csv("Raw/brainTags_Darmanis2015.txt")
  rownames(darmanis) <- darmanis[,1]
  darmanis <- darmanis[,-1]
  
  # add ensID
  darmanis <- addENSID(darmanis)
  
  # annotate
  darmanis_meta <- read.table("Raw/SraRunTable_Darmanis2015.txt", sep = "\t", header = 1)
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
  
## Mixture creation
  nMix <- 100 # number of cells to mix per sample
  nReps <- 100 # number of samples
  
  bench_dm <- create.snmeRand(x = darmanis,
                              nMix = nMix,
                              nReps = nReps,
                              meta = darmanis_meta,
                              meta.columns = "cell_type_s",
                              rpkm = FALSE)
  
  # bench_dm is a list with two levels: 
    # $mixture is a list with dataframes for count- and cpm-level normalisations
    # $true is a list with a dataframe for true proportions, on the RNA content rather than cell content
  
## Save
  # rda
  save(bench_dm, file = "Preprocessed/BenchmarkingDataset_DM.rda")

  # cibersort-compatible scmes
  write.CIB(data = bench_dm$mixture$cpm, dir = "Preprocessed/CIB/Bench_DM.txt")
  
## Note that no gradient sampling approach was used for DM, due to its low number of cells
## Also due to the low number of cells, cells were not split in two for signature and mixture creation

  
################################################################################################################################ #
## Create a benchmarking dataset using the CA and VL single-nuclei  ----    


## Load
  load("Preprocessed/SeuratObjects.rda") # made in Preprocessing_a_SingleNucleusData.R
    
  
## First: prefilter to remove any lowly-abundant cell-types from being used in the simulation. 
  min.celltype.n <- 200 # minimum number of members in a celltype for it to be kept. applied to anything used for creating mixtures (at this stage, Vel and HCA, but the former passes this criterion for all celltypes anyway...)
  
  # for CA
  keep <- table(obj$CA$orig.celltype)
  keep <- names(keep[which(keep > min.celltype.n)])
  # keep <- keep[-grep("no class", keep)]
  
  keep <- which(obj$CA$orig.celltype %in% keep)
  obj$CA <- subset(obj$CA, cells = keep)
  
  # for VL, no cell-types are below this threshold
  
## Setup lists
  # snme <- list() # single-nucleus mixture experiments
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
  bench_CA <- create.snmeRand(x = dat,
                             nMix = nMix,
                             nReps = nReps,
                             meta = obj$CA@meta.data[cellsForX$CA.mix,],
                             meta.columns = "orig.celltype",
                             rpkm = FALSE)
  
  # velmeshev
  dat <- as.data.frame(obj$VL@assays$RNA@counts)
  dat <- dat[,cellsForX$VL.mix]
  bench_VL <- create.snmeRand(x = dat,
                             nMix = nMix,
                             nReps = nReps,
                             meta = obj$VL@meta.data[cellsForX$VL.mix,],
                             meta.columns = "orig.celltype",
                             rpkm = FALSE)
  
  
  
  ## Save
  write.CIB(data = bench_CA$mixture$cpm, dir = "Preprocessed/CIB/Bench_CA.txt")
  write.CIB(data = bench_VL$mixture$cpm, dir = "Preprocessed/CIB/Bench_VL.txt")
  
  save(bench_CA, bench_VL, file = "Preprocessed/bench_CA_VL.rda")  
  
  
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
    
  ## Set proportions of excitatory cells: balanced abundance of ct
    norm.abundance <- aggregate(meta$nCount_RNA, list(meta$grad1), mean)
    # norm.abundance$x <- max(norm.abundance$x) / norm.abundance$x
    norm.abundance$OverMax <- max(norm.abundance$x) / norm.abundance$x
    norm.abundance$OverMin <- norm.abundance$x / min(norm.abundance$x)
    
    base.n <- (nMix*2) / sum(norm.abundance$OverMax)
    
    norm.n <- floor(norm.abundance$OverMax * base.n)
    norm.n <- lapply(norm.n, function(y) { data.frame(c(0, y, 0, y)) })
    norm.n <- do.call("cbind", norm.n)
    colnames(norm.n) <- norm.abundance$Group.1
    
    # change cellnames for this run
    colnames(cells) <- meta$grad1
    
  ## Run
    grad_VL <- confound.proportion(dat = cells, # this balanced methods uses the same simulation approach as for generating datasets with confounded co
                                   ct = unique(colnames(cells)), 
                                   setProps = norm.n, 
                                   nPerGroup = 50, 
                                   nPerMixture = 500,
                                   keepPureExpression = FALSE)
    
    
    # cpm normalise
    
    grad_VL$confound.p <- apply(grad_VL$confound.p, 2, function(y) {
        lib.size <- 10^6 / sum(y)
        y <- y * lib.size    
    })
    
  ## Save
  save(grad_VL, file = "Preprocessed/snme_grad_VL.rda")  
  
  
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
    
  ## Set proportions: balanced abundance of ct
    norm.abundance <- aggregate(meta$nCount_RNA, list(meta$grad1), mean)
    # norm.abundance$x <- max(norm.abundance$x) / norm.abundance$x
    norm.abundance$OverMax <- max(norm.abundance$x) / norm.abundance$x
    norm.abundance$OverMin <- norm.abundance$x / min(norm.abundance$x)
    
    base.n <- (nMix*2) / sum(norm.abundance$OverMax)
    
    norm.n <- floor(norm.abundance$OverMax * base.n)
    norm.n <- lapply(norm.n, function(y) { data.frame(c(0, y, 0, y)) })
    norm.n <- do.call("cbind", norm.n)
    colnames(norm.n) <- norm.abundance$Group.1
    
    # change cellnames for this run
    colnames(cells) <- meta$grad1
  
  # run
  
    grad_CA <- confound.proportion(dat = cells, # this balanced methods uses the same simulation approach as for generating datasets with confounded co
                                   ct = unique(colnames(cells)), 
                                   setProps = norm.n, 
                                   nPerGroup = 50, 
                                   nPerMixture = 500,
                                   keepPureExpression = FALSE)
    
    
    # cpm normalise
    
    grad_CA$confound.p <- apply(grad_CA$confound.p, 2, function(y) {
      lib.size <- 10^6 / sum(y)
      y <- y * lib.size    
    })
    
    ## Save
    save(grad_CA, file = "Preprocessed/snme_grad_CA.rda")  
  
  
  
################################################################################################################################ #
## Create signatures from the above Seurat objects ----
  

merges <- list()
drops <- list()
  
  
## In this section, special versions of the CA and VL signatures are made, for deconvolving their respective mixtures
  
  
## First, get base signature with all cell-types
  for(j in c("VL", "CA")) {
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
    # x <- x[-grep("Full.cpm", names(x))]
    
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
  
  ## Merge celltypes in Velmeshev. each level denotes which cell-types are pooled, and renamed per the level's name
    merges$VL <- list()
  
    # main cell-types
    merges$VL$Neurons <- c("IN-PV", "IN-SST", "IN-SV2C", "IN-VIP", "L2/3", "L4", "L5/6", "L5/6-CC", "Neu-mat", "Neu-NRGN-I", "Neu-NRGN-II")
    merges$VL$Astrocytes <- c("AST-FB", "AST-PP")
    merges$VL$Endothelial <- "Endothelial"
    merges$VL$Microglia <- "Microglia"
    merges$VL$Oligodendrocytes <- "Oligodendrocytes"
    merges$VL$OPCs <- "OPC"
    
    # sub-types of neurons
    merges$VL$Neurons_Inh <- c("IN-PV", "IN-SST", "IN-SV2C", "IN-VIP")
    merges$VL$Neurons_Exc <- c("L2/3", "L4", "L5/6", "L5/6-CC", "Neu-mat")
    merges$VL$Neurons_NRGN <- c("Neu-NRGN-I", "Neu-NRGN-II")
    
    # apply
    sigsSNME <- merger(sig = "VL",
                       levels = c("Neurons",
                                  "Astrocytes",
                                  "Oligodendrocytes",
                                  "Microglia",
                                  "Endothelial",
                                  "OPCs"),
                       name = "merge1")
    sigsSNME <- merger(sig = "VL",
                       levels = c("Neurons_Exc",
                                  "Neurons_Inh",
                                  "Neurons_NRGN",
                                  "Astrocytes",
                                  "Oligodendrocytes",
                                  "Microglia",
                                  "Endothelial",
                                  "OPCs"),
                       name = "merge2")
   
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
  sigsSNME <- merger(sig = "CA",
                     levels = c("Neurons",
                                "Astrocytes",
                                "Oligodendrocytes",
                                "OPCs"),
                     name = "merge1")
  sigsSNME <- merger(sig = "CA",
                     levels = c("Neurons_Exc",
                                "Neurons_Inh",
                                "Astrocytes",
                                "Oligodendrocytes",
                                "OPCs"),
                     name = "merge2")
  

    
## Save
save(sigsSNME, file = "Preprocessed/Signatures - SNME.rda")
  
