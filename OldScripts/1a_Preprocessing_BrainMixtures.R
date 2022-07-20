## In this script, we process mixtures derived from human brain samples

################################################################################################################################ #
## Setup ----

## Generic
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)


## Set directory
# Not setting the working directory here for once... It's set individually in each section

## Functions and libraries
source("/Volumes/Data1/PROJECTS/BrainCellularComposition/Scripts/Fun_Preprocessing.R")
source("/Volumes/Data1/PROJECTS/BrainCellularComposition/Scripts/Fun_Composition.R")
source("/Volumes/Data1/PROJECTS/Lister_brainTimecourse/10X_Pipeline.R")
load("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/geneInfo.rda")
load("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/exonicLength.rda")


################################################################################################################################ #
## Read in and preprocess mixtures ----

  
setwd("/Volumes/Data1/DATA/2021/Lister_CrushedBrains/")

## Read in sample information
  meta.cr <- list.files("FASTQ/OriginalDownload/") # get the filenames of the original data Ryan sent
  meta.cr <- meta.cr[seq(1, length(meta.cr), by = 4)] # we only need every fourth entry, as the others are redundant
  meta.cr <- gsub("brain_", "", meta.cr) # removes the word brain from crushed samples, helps with downstream formatting 
  
  meta.cr <- strsplit(meta.cr, "_")
  meta.cr <- do.call("rbind", meta.cr)
  meta.cr <- as.data.frame(meta.cr)
  colnames(meta.cr) <- c("Boss", "Library", "Year", "Month", "Date", "BrainID", "Prep", "Sequencing1", "Sequencing2", "SampleID", "Lane", "Read", "Suffix")
  meta.cr <- meta.cr[,c("SampleID", "BrainID", "Library", "Prep")]
  # meta.cr$Million.Mapped.Reads <- colSums(crushed) / 10^6 # millions of reads mapped per cell
  meta.cr$Prep[grep("crushed", meta.cr$Prep)] <- "cell" # renaming 
  rownames(meta.cr) <- paste0(meta.cr$BrainID, "_", substr(meta.cr$Prep, 1, 1))

  # manually read RINe off the pdf Daniel sent
  meta.cr$RINe <- rep(c(2.6, 3.8, 4.3, 3.6, 2.9, 3), 2)
  
  
## Read in STAR summary data
  # nuclear samples
  # star.nuc <- read.table("STAR_ERCC_premRNA_Output/STARsummary.txt", sep = "\t", header = TRUE, row.names = 1)
  # colnames(star.nuc) <- substr(colnames(star.nuc), 6, 7)
  # 
  # whole cell samples
  star <- lapply(meta.cr$SampleID, function(x) {
    # read in
    dir <- paste0("STAR_ERCC_NoDups_Output/", x, "/", x, "Log.final.out")
    tab <- read.table(dir, sep ="\t", as.is = TRUE, fill = TRUE)
    
    # clean
    # tab$V1 <- gsub(" ", )
    tab$V1 <- trimws(tab$V1) # removes white space
    tab$V1 <- gsub("\\|", "", tab$V1)
    tab$V1 <- trimws(tab$V1) 
    
    tab <- data.frame(x = tab[,2], row.names = tab[,1])
    colnames(tab) <- x
    
    return(tab)
  })
  
  star <- do.call("cbind", star)
  
  write.csv(star, "STAR_ERCC_NoDups_Output/STAR Summary.csv")
  
  ## Append important STAR summary data to the metadara
    # library size
    meta.cr$LibSize <- as.numeric(star["Uniquely mapped reads number",]) / 10^6
    
    # uniquely mapped reads %
    meta.cr$pct.uniquely.mapped.reads <- as.numeric(substr(star["Uniquely mapped reads %",], 1, 4)) 
    
  
## Read in expression data
  counts <- lapply(meta.cr$SampleID, function(x) {
    y <- read.table(paste0("STAR_ERCC_NoDups_Output/", x, "/", x, "ReadsPerGene.out.tab"), sep = "\t") # note the lack of premrna
    y <- y[-c(1:4),] # removes unmapped reads
    y <- data.frame(Count = y$V4, row.names = y$V1) # convert to dataframe
    colnames(y) <- x # label the column by the sample number
    return(y)
  })
  
  counts <- do.call("cbind", counts) 
  

## Extract ERCC data
  # extract
  ercc.rows <- grep("ERCC", rownames(counts))
  ercc.counts <- counts[ercc.rows,]
  counts <- counts[-ercc.rows,]
  
  # annotate the metadata
  meta.cr$ercc.reads <- colSums(ercc.counts)
  meta.cr$ercc.percent <- meta.cr$ercc.reads / (meta.cr$LibSize * 10^6)
  
    
## Convert the ensembleID to your preferred form
  new.ids <- sapply(strsplit(rownames(counts), "\\."), "[", 1)
  
  # first, any gene that has a duplicated name (even original)
  dup.genes <- new.ids[which(duplicated(new.ids))]
  dup.genes <- which(new.ids %in% dup.genes)
  
  counts <- counts[-dup.genes,]
  new.ids <- new.ids[-dup.genes]
  
  # set
  rownames(counts) <- new.ids
  
## Rename samples per the metadata
  colnames(counts) <- colnames(ercc.counts) <- rownames(meta.cr)
  
## Now: we want to remove C55
  rem <- grep("55", rownames(meta.cr))
  
  counts <- counts[,-rem]
  ercc.counts <- ercc.counts[,-rem]
  meta.cr <- meta.cr[-rem,]

## Normalise
  ## First, load gene length information
    # get common genes
    common <- rownames(counts)[which(rownames(counts) %in% rownames(exonicLength))]
    counts <- counts[common,] 
    l <- exonicLength[common,]
    l <- do.call("c", l)
    
    # convert to kb
    m <- match(rownames(counts), names(l))
    l <- l[m]/1000
    
 
  ## Now normalise!
    norm <- list()
    
    # raw counts
    norm$counts <- counts
    
    # cpm
    norm$cpm <- apply(counts, 1, function(x) {
      x <- x / meta.cr$LibSize
      return(x)
    })
    norm$cpm <- as.data.frame(t(norm$cpm))
    
    # rpkm
    norm$rpkm <- apply(norm$cpm, 2, function(x) x / l)
    norm$rpkm <- as.data.frame(norm$rpkm)
    
    # tpm
    norm$tpm <- tpm(counts)
    
    # nuclear-count normalisation
      # first, threshold the expression data
      dat <- norm$cpm # get cpm data
      keep.genes <- rownames(norm$rpkm)[which(rowSums(norm$rpkm > 1) > 3)] # greater than 1rpkm in at least three samples (yes, rpkm)
      dat <- dat[keep.genes,]
      
      # run linear model
      dat <- t(log2(dat + 0.5))
      prep <- c(rep("nuclei", 5), rep("cell", 5))
      
      mod <- lm(dat ~ prep)
      sum <- summary(mod)
      
      pvals <- sapply(sum, function(y)  { y$coefficients["prepnuclei", "Pr(>|t|)"] } ) # fc is nuclear / wholecell
      log2fc <- sapply(sum, function(y)  { y$coefficients["prepnuclei", "Estimate"] } )
      res <- data.frame(pvals = pvals, log2fc = log2fc)
      rownames(res) <- sapply(strsplit(rownames(res), " "), "[", 2)
      
      res$padj <- p.adjust(res$pvals, method = "fdr")
      nuc.specific.genes <- rownames(res)[which(res$log2fc > log2(10) & res$padj < 0.05)]
      
      # collect the counts for these nuclear specific genes to use as a library size
      meta.cr$NucGeneCount <- colSums(counts[nuc.specific.genes,] / 10 ^ 6) # necessary to use counts here
      
      # generate normalised expression
      norm$nuc.count <- apply(counts, 1, function(x) {
        x <- x / meta.cr$NucGeneCount
        return(x)  
      }) 
      norm$nuc.count <- as.data.frame(t(norm$nuc.count))
      
## Final expression threshold
  min.exp <- 1
  n <- 10 # in all samples
  
  keep.genes <- which(rowSums(norm$rpkm > min.exp) >= n) # we are thresholding based on rpkm data...
  
  # ...and using this threshold for all normalisations
  norm <- lapply(norm, function(x) {
    x <- x[keep.genes,]  
    return(x)
  })
  
  ## ERCC-based expression thresholds
  ercc.info <- read.table("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/ERCC_Information.txt", sep = "\t", header = TRUE, row.names = 1)
  
  # molar ratios of the subgroups in Mix1:Mix2
  # A = 4
  # B = 1
  # C = 0.67
  # D = 0.5
  
  cor(ercc.counts)
  
## Add MALAT1 to metadata
  meta.cr$MALAT1 <- as.numeric(norm$rpkm["ENSG00000251562",])
  

## Save
  # expression data rda
  save(norm, file = "/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/Normalised RNAseq.rda")
  
  # metadata csv
  write.csv(meta.cr, file = "/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/Metadata.csv")
  
  x <- meta.cr[,-c(1,3,5)]
  write.csv(x, file = "/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/SuppTable8_Part1.csv")
  
  ## CIBERSORT-compatible data
    # all samples
    write.CIB(norm$rpkm, "/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/CIB/Crushed_all_rpkm.txt")
    write.CIB(norm$tpm, "/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/CIB/Crushed_all_tpm.txt")
    write.CIB(norm$cpm, "/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/CIB/Crushed_all_cpm.txt")
  
    # output with just a single individual, for individual-specific deconvolution
    for (j in meta.cr$BrainID) {
      write.CIB(norm$rpkm[,grep(j, colnames(norm$rpkm))], paste0("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/CIB/Crushed_", j, "_rpkm.txt"))
      write.CIB(norm$tpm[,grep(j, colnames(norm$tpm))], paste0("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/CIB/Crushed_", j, "_tpm.txt"))
      write.CIB(norm$cpm[,grep(j, colnames(norm$cpm))], paste0("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/CIB/Crushed_", j, "_cpm.txt"))
    }
  
  
## Basic QC: Was the nuclear extraction successful? Check the nuclear transcript MALAT1
  plot.data <- norm$rpkm["ENSG00000251562",]
  plot.data <- melt(plot.data)
  plot.data$Prep <- meta.cr$Prep
  
  pdf(file = "/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/QC - Nuclear Transcripts.pdf", height = 2.5, width = 8)
  ggplot(plot.data, aes(x = variable, colour = Prep, y = value)) +
    geom_point() +
    theme_bw() +
    scale_y_continuous(limits = c(0, 300), expand = c(0,0)) +
    labs(y = "MALAT1 Expression (RPKM)", y = "Preparation") +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none", axis.title.x = invis)
  dev.off()
  
 

## Basic QC: Histogram of expression
  plot.data <- log2(norm$rpkm + 0.5)
  # plot.data$Gene <- rownames(lister)
  plot.data <- melt(plot.data)
  
  pdf(file = "/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/QC - Expression Distribution.pdf", height = 4, width = 8)
  ggplot(plot.data, aes(x = variable, colour = variable, fill = variable, y = value)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.9) +
    geom_boxplot(fill = "white", outlier.colour = "black", width = 0.2, position = position_dodge(width = 0.8)) +
    theme_bw() +
    labs(y = "log2(RPKM + 0.5)") +
    geom_hline(yintercept = median(plot.data$value), linetype = 2) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none", axis.title.x = invis)
  dev.off()  
  
## Comparing Nuclei vs Cells: compare matched nuclear and cell expression, and overlay markers!
  load("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/SN/Lister Signature, Integrated, Nagy Annot.rda")
  listerXnagy <- addENSID(listerXnagy)
  
  ## An aside on marker genes!
  marks <- nTopFeatures(signature = listerXnagy, n = 100, alg = "diff")
  
  scatterplot <- function(ind, ct) { 
    dat <- data.frame(log2(norm$rpkm[,grep(ind, colnames(norm$rpkm))] + 0.5), Marker = FALSE)
    dat$Marker[which(rownames(norm$rpkm) %in% marks$EnsID[which(marks$forCelltype == ct)])] <- TRUE
    
    colnames(dat) <- c("N", "C", "Marker")
    
    dat <- dat[order(as.numeric(dat$Marker)),]
    
    dat$Marker[which(dat$Marker)] <- ct
    dat$Marker <- factor(dat$Marker, levels = c(FALSE, ct))
    
    ggplot(dat, aes(x = N, y = C, col = Marker)) +
      geom_point() +
      theme_bw() +
      geom_abline(intercept = c(-1,1), col = "red") + 
      scale_colour_manual(values = c("black", "dodgerblue")) +
      labs(x = "Nuclei", y = "Cells", title = paste0(ind, " ", ct)) +
      theme(legend.position = "none")
  }
  
  # inds <- meta.cr$BrainID[1:6]
  inds <- c("4899", "5144", "1541", "5387", "5006")
  
  pdf(file = "/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/N vs C - Marker Enrichment Across Preparations (rpkm).pdf", height = 6, width = 8)
  for (j in c("Neu", "Ast", "Oli", "OPC")) {
    plot.list <- lapply(inds, scatterplot, ct = j)  
    print(plot_grid(plotlist = plot.list, nrow = 2))
  }
  dev.off()
  
  
## Comparing nuclei and cells: fold change
  # distribution
  fc.hist <- function(ind) { 
    # dat <- data.frame(log2(norm$rpkm[,grep(ind, colnames(norm$rpkm))] + 0.5))
    dat <- data.frame(norm$rpkm[,grep(ind, colnames(norm$rpkm))])
    colnames(dat) <- c("Nuc", "Cel")
    dat$FC <- log10(dat$Nuc / dat$Cel)
    # dat <- dat[-which(is.na(dat$FC)),]
    
    ggplot(dat, aes(x = FC)) +
      geom_histogram(colour = "black", fill = "darkorange1", binwidth = 0.2 ) +
      theme_bw() +
      labs(x = "log10 (Nuclei / Cells) (rpkm)", y = "Count", title = paste0(ind)) +
      # labs(y = "log2 (Nuclei / Cells) (rpkm)", x = "log2(Nuclear RPKM)", title = paste0(ind)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "red") +
      # geom_hline(yintercept = 0, linetype = 2, colour = "red") +
      theme(legend.position = "none")
  }
  
  # relationship with expression
  fc.expression <- function(ind) { 
    # dat <- data.frame(log2(norm$rpkm[,grep(ind, colnames(norm$rpkm))] + 0.5))
    dat <- data.frame(norm$rpkm[,grep(ind, colnames(norm$rpkm))])
    colnames(dat) <- c("Nuc", "Cel")
    dat$FC <- log10(dat$Nuc / dat$Cel)
    # dat <- dat[-which(is.na(dat$FC)),]
    
    ggplot(dat, aes(x = log2(Nuc + 0.5), y = FC)) +
      geom_point() +
      # geom_histogram(colour = "black", fill = "darkorange1", binwidth = 0.2 ) +
      theme_bw() +
      # labs(x = "log10 (Nuclei / Cells) (rpkm)", y = "Count", title = paste0(ind)) +
      labs(y = "log2 (Nuclei / Cells) (rpkm)", x = "log2(Nuclear RPKM)", title = paste0(ind)) +
      # geom_vline(xintercept = 0, linetype = 2, colour = "red") +
      geom_hline(yintercept = 0, linetype = 2, colour = "red") +
      theme(legend.position = "none")
  }
  
  pdf(file = "/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/N vs C - Fold Change in Nuclei.pdf", height = 6, width = 8)
  plot.list <- lapply(inds, fc.hist)  
  print(plot_grid(plotlist = plot.list, nrow = 2))
  
  plot.list <- lapply(inds, fc.expression)  
  print(plot_grid(plotlist = plot.list, nrow = 2))
  dev.off()
  

################################################################################################################################ #
## Load signature data from snRNA-seq ----
  

setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/")


## Get samplenames and metadata
  # get a list of files
  sn.filelist <- c("/Volumes/Data1/DATA/2019/FragileX_snRNAseq/RL1390/filtered_gene_bc_matrices_h5.h5", # library RL1390, from brain F5006
                   "/Volumes/Data1/DATA/2020/Lister_ASD_snRNAseq/filtered_h5/RL2101_14yr_ASD_v3_FANS_filtered.h5", # library RL2101, from brain A4899
                   "/Volumes/Data1/DATA/2020/Lister_ASD_snRNAseq/filtered_h5/RL2123_20yr_Ctrl_v3_FANS_filtered.h5", # library RL2123, from brain C1541
                   "/Volumes/Data1/DATA/2020/Lister_ASD_snRNAseq/filtered_h5/RL2128_20yr_Ctrl_v3_FANS_filtered.h5", # library RL2128, from brain C55
                   "/Volumes/Data1/DATA/2020/Lister_ASD_snRNAseq/filtered_h5/RL2127_12yr_Ctrl_v3_FANS_filtered.h5", # library RL2127, from brain C5387 
                   "/Volumes/Data1/DATA/2020/Lister_ASD_snRNAseq/filtered_h5/RL2001_7yr_ASD_v3_sucrose_filtered.h5") # library RL2001, from brain A5144 
  
  names(sn.filelist) <- c("F5006", "A4899", "C1541", "C55", "C5387", "A5144")
  
  

  # get a list of sample names and metadata
  # meta <- data.frame(id = splitter(x = files, split = "_", side = 1),
  #                    age = as.numeric(gsub("yr", "", splitter(x = files, split = "_", side = 2))),
  #                    ASD = splitter(x = files, split = "_", side = 3),
  #                    chemistry = splitter(x = files, split = "_", side = 4),
  #                    isolation = splitter(x = files, split = "_", side = 5),
  #                    path = files,
  #                    row.names = splitter(x = files, split = "_", side = 1))
  # 
  # meta$col <- c(brewer_pal(palette = "Set3")(12), "black")
  
## Read in as Seurat objects
  obj <- list()
  for(j in names(sn.filelist)) {
    obj[[j]] <- as.Seurat(file = sn.filelist[j], 
                          h5 = TRUE, 
                          sample.id = j)
  }
  
  
## QC
  obj <- lapply(obj, qc.features)
  qc.plots(obj, "SN Preliminary QC Plots.pdf", colours = brewer_pal(palette = "Set3")(6))

## Normalisation and clustering 
  obj <- lapply(obj, norm.bySeurat)

## Linear dimensionality reduction (i.e. PCA)
  # run PCA
  obj <- lapply(obj, function(x) { RunPCA(object = x, npcs = 50) })
    
  # determine statistically-significant components using an instant heuristic (as our within-sample analyses are but a first pass) 
  plot.list <- lapply(obj, function(x) ElbowPlot(x, ndims = 50))
  plot_grid(plotlist = plot.list, ncol = 3)
  
  # for simplicity, I will use the same dimensionality (i.e. PC cut-off) for all samples
  # it is recommend to err on the high side, so:
  dim <- 35
  
## Clustering and visualisation, at ~10s per level
  obj <- lapply(obj, function(x) clust.bySeurat(x, dim = dim, resolution = 1.5)) # recommendation is 0.4-1.2 for small datasets (3k cells)

## Define cluster identity from marker enrichment
  load("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/snme_Signatures.rda")
  # sigs.snme <- sigs.snme[1]
  # sigs.snme$Vel <- sigs.snme$Vel[-c(1,2)] 
  
  obj <- lapply(obj, function(x) clust.annot(x, assay = "RNA", sigs.snme$Vel$orig.celltype))
  
## Save
  save(obj, file = "Preliminary Object.rda")
 
## Visualise 
clust.plots.list(obj, directory = "Preliminary tSNE Vel V2.pdf", ncol = 3, nrow = 2, height = 8, width = 12)
  
## More visualising
  ## Define markers to visualise
  markers <- list(Exc = c("NRG1", "RORB", "CUX2", "RBFOX3"),
                  Inh = c("GAD1", "SST", "VIP", "SV2C"),
                  Oli = c("OLIG2", "PLP1", "MBP", "MAG"),
                  OPC = c("PDGFRA"),
                  Ast = c("GFAP", "SLC2A2", "CD44"))
  
  for(i in names(obj)) {
    pdf(file = paste0("Marker Enrichment tSNE ", i, ".pdf"), height = 9, width = 8)
    plot.grid <- list()
    for(j in names(markers)) {
      plot.list <- list()
      for(k in markers[[j]]) {
        if (!(k %in% rownames(obj[[i]]@assays$RNA@counts))) next
        
        plot.list[[k]] <- FeaturePlot(obj[[i]], features = k) +
          theme(legend.position = "none", axis.title = invis, axis.text = invis, axis.ticks = invis) +
          labs(title = paste0(j, " (", k, ")"))
        
      }
      plot.grid[[j]] <- plot_grid(plotlist = plot.list, ncol = 4, nrow = 1)
    }  
    print(plot_grid(plotlist = plot.grid, ncol = 1, nrow = 5))
    dev.off()
  }
  
## Let's push on to integrating the data instead...
  anchors <- FindIntegrationAnchors(object.list = obj, 
                                    dims = 1:30, # in the tutorial, dims = 1:20 was recommended, but here I'm using default parameters.
                                    scale = FALSE) # scale is set to FALSE, as I've already scaled and regressed the data. per the help: "Only set to FALSE if you have previously scaled the features you want to use for each object in the object.list"
  
  # Warning message:
  # In CheckDuplicateCellNames(object.list = object.list) :
  # Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.
  
  ## Integrate using anchors
  int <- IntegrateData(anchorset = anchors, dims = 1:30)
  
  ## Set assay
  DefaultAssay(int) <- "integrated" # tells subsequent functions to look at int[["integrated"]]@...
  
  ## Scale data
  int <- ScaleData(int) # scale data again. i must debate whether to regress now, or previously. i chose the latter as the covariate of nUMI has a dataset-specific effect
  
  ## PCA, in ~30s
  int <- RunPCA(int, npcs = 50) # default is 50, though the tutorial showcases 30. calculating more doesn't change results, but grants more flexibility later
  
  ## tSNE and Clustering 
    ## First, optimise dimensionality and resolution parameters
  
    ElbowPlot(int, ndims = 50, reduction = "pca")
    dim <- 35 # based on the elbow plot
    
    # tSNE
    int <- RunTSNE(object = int, dims = 1:dim, seed.use = 1)
    
    # cluster
    int <- FindNeighbors(int, reduction = "pca", dims = 1:dim)
    int <- FindClusters(int, resolution = 1) 
    
  ## Label clusters
    int <- clust.annot(int, assay = "integrated", new.column.name = "VelAnnot", sigs.snme$Vel$orig.celltype)
    int <- clust.annot(int, assay = "integrated", new.column.name = "NagyAnnot", sigs.snme$Nagy$orig.celltype)
    int <- clust.annot(int, assay = "integrated", new.column.name = "HCAAnnot", sigs.snme$HCA$orig.celltype)

    # relabel to the level of merge2
    int$NagyAnnot2 <- as.character(int$NagyAnnot)
    int$NagyAnnot2[grep("Oligo", int$NagyAnnot2)] <- "Oli"
    int$NagyAnnot2[grep("OPC", int$NagyAnnot2)] <- "OPC"
    int$NagyAnnot2[grep("Ast", int$NagyAnnot2)] <- "Ast"
    int$NagyAnnot2[grep("Endo", int$NagyAnnot2)] <- "End"
    int$NagyAnnot2[grep("^Ex", int$NagyAnnot2)] <- "Exc"
    int$NagyAnnot2[grep("^Inh", int$NagyAnnot2)] <- "Inh"
    

  ## Plot
    pdf(file = "Integrated Libraries - tSNE.pdf", height = 6, width = 8)
    DimPlot(int, group.by = "seurat_clusters") + labs(title = "By Seurat Cluster")
    DimPlot(int, group.by = "orig.ident") + labs(title = "By Library")
    DimPlot(int, group.by = "VelAnnot") + labs(title = "By Vel")
    DimPlot(int, group.by = "NagyAnnot") + labs(title = "By Nagy")
    DimPlot(int, group.by = "HCAAnnot") + labs(title = "By HCA")
    dev.off()
    
    pdf(file = "Integrated Libraries - tSNE (Supp Fig).pdf", height = 6, width = 8)
    DimPlot(int, group.by = "NagyAnnot2")
    dev.off()
  
    pdf(file = "Integrated Libraries - Markers.pdf", height = 9, width = 8)
    
      for(j in names(markers)) {
        plot.list <- list()
        for(k in markers[[j]]) {
          if (!(k %in% rownames(int@assays$RNA@counts))) next
          
          plot.list[[k]] <- FeaturePlot(int, features = k) +
            theme(legend.position = "none", axis.title = invis, axis.text = invis, axis.ticks = invis) +
            labs(title = paste0(j, " (", k, ")"))
          
        }
        plot.grid[[j]] <- plot_grid(plotlist = plot.list, ncol = 4, nrow = 1)
      }  
      print(plot_grid(plotlist = plot.grid, ncol = 1, nrow = 5))
      dev.off()
      
    
  ## Save
    save(int, file = "Integrated snRNAseq Data.rda")
    
  ## Save signature
    int$NagyMerge <- int$NagyAnnot
    
    # levels(int$NagyMerge)[grep("Ex", levels(int$NagyMerge))] <- "Exc"
    # levels(int$NagyMerge)[grep("Inh", levels(int$NagyMerge))] <- "Inh"
    levels(int$NagyMerge)[grep("Ex|Inh", levels(int$NagyMerge))] <- "Neu"
    levels(int$NagyMerge)[grep("Oli", levels(int$NagyMerge))] <- "Oli"
    levels(int$NagyMerge)[grep("OPC", levels(int$NagyMerge))] <- "OPC"
    levels(int$NagyMerge)[grep("Ast", levels(int$NagyMerge))] <- "Ast"
    
    listerXnagy <- list() 
    for (j in levels(int$NagyMerge)) {
      print(j)
      listerXnagy[[j]] <- rowSums(as.matrix(int@assays$RNA@counts[,which(int$NagyMerge == j)]))
    }
    
    listerXnagy <- do.call("cbind", listerXnagy)
    listerXnagy <- as.data.frame(listerXnagy)
    
    listerXnagy <- apply(listerXnagy, 2, function(x) { # convert to cpm
      libsize <- sum(x) / 10^6
      x <- x / libsize
      return(x)
    })
  
    listerXnagy <- listerXnagy[which(apply(listerXnagy, 1, max) > 1),] # expression threshold
    listerXnagy <- as.data.frame(listerXnagy)
    listerXnagy <- listerXnagy[,-5] # remove End
    
    save(listerXnagy, file = "Lister Signature, Integrated, Nagy Annot.rda")
    
  ## Create signatures from each individual
    ind.sigs <- list()
    for (j in levels(as.factor(int$orig.ident))) {
      x <- list()
      
      for (k in levels(int$NagyMerge)) {
        print(k)
        x[[k]] <- rowSums(as.matrix(int@assays$RNA@counts[,which(int$NagyMerge == k & int$orig.ident == j)]))
      }
      
      x <- do.call("cbind", x)
      x <- as.data.frame(x)
      
      x <- apply(x, 2, function(x) { # convert to cpm
        libsize <- sum(x) / 10^6
        x <- x / libsize
        return(x)
      })
      
      x <- x[which(apply(x, 1, max) > 1),] # expression threshold
      x <- as.data.frame(x)
      x <- x[,-5] # remove End
      
      ind.sigs[[j]] <- x
    }
    
    save(ind.sigs, file = "Lister Signature, Library-specific, Nagy Annot.rda")
    
    
  ## Get true props
    true.props <- list()
    x <- table(int$orig.ident, int$NagyMerge)
    x <- as.matrix(x)
    x <- apply(x, 1, function(y) {y / sum(y)})
  
    write.csv(x, file = "/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/True Proportions, From Integrated Nagy.csv")
    
## Save metadata
  y <- lapply(obj, function(x) {
    a <- ncol(x)
    b <- median(x$nCount_RNA)
    c <- median(x$percent.mito)
    data.frame(nNuclei = a, median.nUMI = b, median.mitoPercent = c)
  })
  y <- do.call("rbind", y)
  
  write.csv(y, "/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/Lister_CrushedBrains/SuppTable8_Part2.csv")

    
################################################################################################################################ #
## Preprocess Zhang mixtures ----

## Downloaded as fpkm data in Table S4, Zhang et al. 2016, and modified to keep only healthy, human, pure cell-type transcriptomes
  zhang <- read.csv("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Raw/Zhang2016_RNASeq.csv")

  dup <- which(duplicated(zhang$GeneSymbol)) 
  zhang <- zhang[-dup,]

## reannotate gene names
  rownames(zhang) <- zhang$GeneSymbol
  zhang <- zhang[,-which(colnames(zhang) == "GeneSymbol")]
  zhang <- addENSID(zhang)

# filter out cell-types
  zhang <- zhang[,-grep("Foetal", colnames(zhang))]
  
# apply minimum threshold of 1 fpkm in any one cell-type
  exp_thresh <- 1
  zhang <- zhang[which(apply(zhang, 1, max) > exp_thresh),] # 12k genes is reasonable. only necessary in one sample as there's only one neuronal sample...

# write to disk
  write.CIB(zhang, "/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/CIB/Zhang_PureSamples.txt")

    
    
