################################################################################################################################ #
## Libraries ----

## For loading in data
  require(Seurat)

## For normalisation/preprocessing
  require(celda) # include decontX  

## For dataset integration
  require(harmony)
  require(SingleR)
  # library(Seurat) # already loaded!

## For plots
  require(ggplot2)
  require(reshape2)
  require(cowplot)

################################################################################################################################ #
## Generic functions ----

  
splitter <- function(
  x, # x is a character vector
  split, # split is the character at which the split shall occur (and a split occurs at EVERY match)
  side # side is which side of the split to be returned, any value from 1 - X where X is the number of splits + 1
) { 
  sapply(strsplit(x, split = split), "[", side)
}

################################################################################################################################ #
## Functions: loading in data and preliminary QC ----

## Function to load output from cellranger count
as.Seurat <- function(
  file, # a character vector pointing to EITHER the X_feature_bc_matrix directory, OR the .h5 file, which output from cellranger count
  sample.id, # a name to bestow upon the sample
  min.cells = 3,
  min.features = 200,
  h5 = TRUE
) {
  
  # read in raw data
  cat(paste0("\nReading in raw data for ", sample.id, " :)"))
  if (h5) {
    obj <- Read10X_h5(filename = file)
  } else {
    obj <- Read10X(data.dir = file)
  }
  
  # convert to Seurat object
  cat("\nDone! Now creating a Seurat object")
  obj <- CreateSeuratObject(
    counts = obj,
    min.cells = min.cells,
    min.features = min.features,
    project = sample.id
  )
  
  # return
  return(obj)
}

## Generic QC function
qc.features <- function(
  obj, # a seurat object, typically the output of as.Seurat
  cc = FALSE
) {
  
  # calculate mitochondrial percetage
  obj[["percent.mito"]] <- PercentageFeatureSet(object = obj, pattern = "^MT-") 
  
  # calculate ribosomal percentage
  obj[["percent.rps"]] <- PercentageFeatureSet(object = obj, pattern = "^RPS") 
  obj[["percent.rpl"]] <- PercentageFeatureSet(object = obj, pattern = "^RPL") 
  
  # calculate MALAT1 (nuclear marker) percentage
  obj[["percent.malat1"]] <- PercentageFeatureSet(object = obj, pattern = "MALAT1") 
  
  # a cell-cycle score
   if(cc) obj <- CellCycleScoring(obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
  
  # return
  return(obj)
}

## Plot QC metrics
qc.plots <- function(
  obj, # a seurat object that has passed through qc.features
  directory, # where the plots shall be saved!
  plot.width = 8, # set the pdf width
  colours,
  log10 = TRUE # whether or not to log transform data
) {
  
  ## Get plot data
    if (class(obj) == "list") {
      plot.data <- lapply(obj, function(x) { cbind(x@meta.data[,-1]) } )
      plot.data <- do.call("rbind", plot.data)
      plot.data$Library <- sapply(strsplit(rownames(plot.data), "\\."), "[", 1)
    } else {
      plot.data <- obj@meta.data
      plot.data$Library <- "obj"
    }
    
    colnames(plot.data)[1:2] <- c("nUMI", "nGene")
    if (log10) plot.data$nUMI <- log10(plot.data$nUMI + 1)
    if (class(obj) == "list") plot.data <- reshape2::melt(plot.data, id.vars = "Library")

  ## Plot
  plot.fun <- function(variable) {
    plot <- ggplot(plot.data[plot.data$variable == variable,], aes(x = Library, y = value, colour = Library)) +
      geom_violin(aes(fill = Library), scale = "width", width = 0.8, position = position_dodge(width = 0.8)) +
      geom_boxplot(fill = "white", outlier.colour = "black", width = 0.15, position = position_dodge(width = 0.8)) +
      theme_classic() +
      scale_colour_manual(values = colours) +
      scale_fill_manual(values = colours) +
      theme(panel.border = element_blank(), panel.background = element_rect(fill = "grey90"), 
            axis.line = element_line()) +
      labs(y = variable)
    
    if (class(obj) != "list") plot <- plot + theme(axis.title.x = element_blank())
    
    print(plot)
  }
  
  pdf(file = directory, height = 4, width = plot.width) 
  plot.fun("nUMI")
  plot.fun("nGene")
  plot.fun("percent.mito")
  plot.fun("percent.rps")
  plot.fun("percent.rpl")
  plot.fun("percent.malat1")
  dev.off()
  
  return("Done!")
}
  
################################################################################################################################ #
## Functions: normalisation ----

## Removal of ambient RNA by decontX
decontX.Seurat <- function(
  obj,
  label = "seurat_clusters", # the column name of the Seurat metdata to use as a grouping variable
  seed = 12345 # the default in decontX
) {
  # process counts
  mat <- as.matrix(obj@assays$RNA@counts)
  storage.mode(mat) <- "integer"
  
  # process clustering
  clust <- as.character(obj@meta.data[,label])
  
  # run
  res <- decontX(counts = dat, z = clust)
  
  # end
  return(res)
}

## Local perturbation
local.pert <- function(z, reduction, control.pool, n.neighbours) {
  # get normalised expression data
  test <- as.data.frame(z@assays$RNA@data[,-control.pool])
  con <- as.data.frame(z@assays$RNA@data[,control.pool])
  
  # get euclidian distance matrix from every test cell to every control cell
  tsne_test <- z@reductions$tsne@cell.embeddings[-control.pool,]
  tsne_con <- z@reductions$tsne@cell.embeddings[control.pool,]
  
  dist <- apply(tsne_test, 1, function(x) { 
    t1 <- (x[1] - tsne_con[,1]) ^ 2
    t2 <- (x[2] - tsne_con[,2]) ^ 2
    sqrt(t1 + t2)
  })
  
  # get the nearest negative control neighbours for each test cell
  nn <- apply(dist, 2, function(x) { # second argument must be 2
    y <- order(x, decreasing = FALSE)[1:n.neighbours]
    rownames(dist)[y]
  })
  
  nn <- as.data.frame(nn)
  
  # normalise expression
  reg <- test
  reg[,] <- 0
  for (j in 1:ncol(test)) {
    reg[,j] <- test[,j] - rowMeans(con[,nn[,j]])
  }
  
  return(reg)
}

## An as yet uncompleted normalisation function! TAKES ROUGHLY THIS FORM
norm.bySeurat <- function(
  obj, # a Seurat object
  min.nUMI = 500, # remove cells/nuclei with nUMI less than this
  max.percent.mito = 20, # remove cells/nuclei with percent.mito greater than this
  max.nUMI.percentile = 0.998 # remove cells/nuclei with nUMI greater than this percentile of all cells/nuclei
) {
  
  # subset cells
  topCount <- quantile(obj@meta.data$nCount_RNA, probs = max.nUMI.percentile)
  obj <- subset(x = obj, subset = (nCount_RNA > min.nUMI & nCount_RNA < topCount & percent.mito < max.percent.mito))
  
  # obj <- subset(x = obj, subset = nCount_RNA > 500 & percent.mito < 5)
  
  # normalise expression levels
  obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000) # standard parameters for Seurat

  # find variable genes (i.e. features)
  obj <- FindVariableFeatures(object = obj, selection.method = "vst", nfeatures = 2000)
  
  # scale data, and regress out covariates of nUMI and percent.mito. takes ~2 min. 
  # note that FindVariableFeatures returns identical results whether or not the data are scaled first
  obj <- ScaleData(object = obj, vars.to.regress = c("nCount_RNA", "percent.mito")) # this is using a Seurat V2 workflow; V3 recommends sctransform()
  
  # output
  return(obj)
}
   

## Clustering and visualisation, at ~10s per level
clust.bySeurat <- function(
  obj, # a Seurat object
  method = "tSNE", # choice 
  resolution = 1, # "We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets" 
  dim # gleaned from the Elbow plot heuristic
) {
  if(!(method %in% c("tSNE", "UMAP"))) {
    stop("Method must be one of tSNE or UMAP. Oh, but UMAP isn't working, so please consider tSNE for now :(")
  }
  
  cat("\n New library started!\n")
  
  if (method == "tSNE") obj <- RunTSNE(object = obj, dims = 1:dim, seed.use = 1)
  if (method == "UMAP") obj <- RunUMAP(object = obj, dims = 1:dim, seed.use = 1)
  
  obj <- FindNeighbors(object = obj, reduction = "pca", dims = 1:dim) 
  obj <- FindClusters(object = obj, resolution = 1) 
  
  return(obj)
}

## Run singleR
clust.annot <- function(
  obj, # a Seurat obj
  assay = "RNA", 
  new.column.name = "snglR.clust",
  sig # the signature of expression to which you wish to compare your data
) {
  
  print("Starting new library")
  
  # compare existing clusters to reference clusters. Correlation-based.
  annot <- SingleR(method = "cluster", 
                   sc_data = obj[[assay]]@data, 
                   ref_data = sig, 
                   types = colnames(sig),
                   clusters = obj$seurat_clusters, 
                   genes = "de")

  # reannotate the Seurat object
  new.ids <-  annot$labels1.thres
  names(new.ids) <- levels(obj)
  obj <- RenameIdents(object = obj, new.ids)
  
  # add a new metadata column for the new annotation
  obj@meta.data[,new.column.name] <- obj$seurat_clusters
  levels(obj@meta.data[,new.column.name]) <- annot$labels1.thres
  
  # return
  return(obj)
}

## Dimensionality reduction plots
clust.plots.list <- function(
  obj, # a seurat object that has passed through qc.features
  directory, # where the plots shall be saved!
  ncol,
  nrow,
  height, # set pdf height
  width # set the pdf width
) {
  
  # get plots
  orig.plots <- lapply(obj, function(x) DimPlot(object = x, group.by = "seurat_clusters"))
  annot.plots <- lapply(obj, function(x) DimPlot(object = x, group.by = "snglR.clust"))
  nUMI.plots <- lapply(obj, function(x) FeaturePlot(object = x, features = "nCount_RNA"))
  
  # rename
  for (j in names(orig.plots)) orig.plots[[j]] <- orig.plots[[j]] + labs(title = j)
  for (j in names(orig.plots)) annot.plots[[j]] <- annot.plots[[j]] + labs(title = j)
  for (j in names(orig.plots)) nUMI.plots[[j]] <- nUMI.plots[[j]] + labs(title = j)
  
  # output
  pdf(file = directory, height = height, width = width)
  print(plot_grid(plotlist = orig.plots, nrow = nrow, ncol = ncol))
  print(plot_grid(plotlist = annot.plots, nrow = nrow, ncol = ncol))
  print(plot_grid(plotlist = nUMI.plots, nrow = nrow, ncol = ncol))
  dev.off()
  return("Done!")
}


    
## Dataset integration using Seurat
int.bySeurat <- function(
  obj.list,
  dims = 1:30, # in the tutorial, dims = 1:20 was recommended, but here I'm using default parameters.
  scale = FALSE # "Only set to FALSE if you have previously scaled the features you want to use for each object in the object.list"
) {
  int <- FindIntegrationAnchors(object.list = obj.list, 
                                  dims = dims, 
                                  scale = FALSE) 

  int <- IntegrateData(anchorset = int, dims = dims)

  return(int)
}
  

## Dataset integration using Harmony


################################################################################################################################ #
## Functions: clustering ----

