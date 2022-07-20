################################################################################################################################ #
## Setup ----

## Generic
rm(list = ls())
gc()
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/")
options(stringsAsFactors = FALSE)

## Load some gene-level information to help with annotation and preprocessing
load(file = "Preprocessed/geneInfo.rda")
load("Preprocessed/exonicLength.rda")

## Functions for estimating composition
source("../Scripts/Fun_Composition.R")
source("../Scripts/Fun_Preprocessing.R")

## Functions for subsetting our large cross-tissue resources to tissues/cell-types of interest
  extract.gtex.tissue <- function(tissue, norm = "rpkm") {
    # filter to tissue samples
    keep <- which(gtexMeta$SMTS == tissue)
    exp <- gtex[,keep]
    meta <- gtexMeta[keep,]
    
    # process expression data
    for (j in 1:ncol(exp)) exp[,j] <- as.numeric(as.character(exp[,j])) # converts character to numeric
    
    if (norm == "rpkm") {
      exp <- rpkm(exp) # rpkm normalisation
    } else {
      exp <- tpm(exp)
    }
    
    n <- min(table(meta$SMTSD)) / 2 # n is set to be half the number of samples of the smallest subgroup of the tissue
    exp <- exp[rowSums(exp >= 1) >= n,] # keep genes expressed at > 1 rpkm/tpm in at least n samples
    
    # return
    return(list(exp = exp, meta = meta))
  }

  extract.encode.signature <- function(samples) {
    # get sample list
    fetch <- encode.meta$File.accession[which(encode.meta$Biosample.term.name %in% samples)]
    meta <- encode.meta[which(encode.meta$Biosample.term.name %in% samples),]
    
    # read in data
    dat <- list()
    for (j in fetch) {
      print(j)
      x <- read.table(paste0("Raw/ENCODE_PrimaryCells/", j, ".tsv"), header = 1)
      rownames(x) <- sapply(strsplit(x$gene_id, "\\."), "[", 1)
      x <- data.frame(x$FPKM, row.names = rownames(x))
      colnames(x) <- j
      
      dat[[j]] <- x
    }
    
    dat <- do.call("cbind", dat)
    
    # aggregate across celltypes of interest
    exp <- list()
    for(j in samples) {
      use <- meta$File.accession[which(meta$Biosample.term.name == j)]
      exp[[j]] <- rowMeans(dat[,use])
    }
    exp <- do.call("cbind", exp)
    exp <- as.data.frame(exp)
    
    # min threshold: where the expression must be above 1 fpkm in at least on cell-type
    exp <- exp[which(apply(exp, 1, max) > 1),]
    
    # return
    return(exp)
  }
  
  extract.F5.signature <- function(samples) {
    # get samples
    keep <- which(F5.meta$Tissue.category %in% samples)
    dat <- F5.exp[,keep]
    meta <- F5.meta[keep,]
    
    # get average expression across celltypes
    # exp <- aggregate(t(dat), by = list(meta$Tissue.category), mean)
    
    exp <- list()
    for(j in levels(as.factor(meta$Tissue.category))) {
      exp[[j]] <- rowMeans(dat[,which(meta$Tissue.category == j)])
    }
    exp <- do.call("cbind", exp)
    exp <- as.data.frame(exp)
    
    #  adding EnsIDs to the annotation!
    exp <- addENSID(exp)
  
    # min threshold: where the expression must be above 1 tpm in at least on cell-type
    exp <- exp[which(apply(exp, 1, max) > 1),]
    
    # return
    return(exp)
  }

################################################################################################################################ #
## Preprocess GTEx mixtures ----


## Raw data
  gtex <- read.table("Raw/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.txt")
  
  # adjust column names
  colnames(gtex) <- unlist(gtex[1,])
  gtex <- gtex[-1,]

  # move gene annotation to rownames
  rownames(gtex) <- sapply(strsplit(as.character(gtex$Name), split = "\\."), `[`, 1)
  gtex <- gtex[,-c(1,2)]
  

  
  ## Metadata
    gtexMeta <- read.csv("Raw/GTEx_v7_Annotations_SampleAttributesDS.csv")
    gtexMeta <- gtexMeta[which(gtexMeta$SAMPID%in%colnames(gtex)),]
    
    write.csv(table(gtexMeta$SMTS), file = "../Results/Revisions/CrossTissue/Samples GTEx Tissue.csv")
    write.csv(table(gtexMeta$SMTSD), file = "../Results/Revisions/CrossTissue/Samples GTEx Regions.csv")
    
    # reorder columns
    m <- match(colnames(gtex), gtexMeta$SAMPID)
    gtexMeta <- gtexMeta[m,]

## Load FANTOM5 data
  F5.exp <- read.csv("Raw/Fantom5Consortium_CAGE.csv", row.names = 1) 
  F5.meta <- read.csv("Raw/Fantom5Consortium_Meta.csv") 
  rownames(F5.meta) <- F5.meta$Library
  
  keep <- colnames(F5.exp)[which(colnames(F5.exp) %in% rownames(F5.meta))]
  F5.meta <- F5.meta[keep,]
  F5.exp <- F5.exp[,keep]
  
  # seems everything is primary cells...
  F5.meta$Tissue.category <- sapply(strsplit(F5.meta$Sample, ","), "[", 1)
  write.csv(table(F5.meta$Tissue.category), file = "../Results/Revisions/CrossTissue/Samples F5 Primary Cells.csv")


## Load ENCODE metadata
  encode.meta <- read.csv("Raw/ENCODE_PrimaryCells/metadata.csv", header = TRUE) # note: I converted the tsv to a csv in excel
  write.csv(table(encode.meta$Biosample.term.name), file = "../Revisions/Results/CrossTissue/Samples ENCODE.csv")
  # note: no need to load ENCODE data, that occurs during the function extract.encode.signature

################################################################################################################################ #
## Heart analyses ----

heart <- list()

## Prepare GTEx data
  heart$gtex <- extract.gtex.tissue("Heart")
  write.CIB(data = heart$gtex$exp, dir = "Preprocessed/CIB/GTEx_Heart.txt")

## Prepare F5 signature
  # locate primary cells of interest
  samples <- c("Cardiac Myocyte", 
               "Endothelial Cells - Aortic", # aorta is different to cardiac, but it's the best in the dataset
               "Fibroblast - Cardiac",
               "Smooth Muscle Cells - Coronary Artery")
  
  heart$F5 <- extract.F5.signature(samples)
  
## Prepare ENCODE signature
  samples <- c("regular cardiac myocyte",
               "endothelial cell of coronary artery",
               "cardiac atrium fibroblast",
               "smooth muscle cell of the coronary artery")
  
  heart$EN <- extract.encode.signature(samples)
  
  
## Prepare single-cell signature from Wang et al. 2020
  heart.sc <- read.csv("Raw/Wang2020_HeartSC_GSE109816_matrix.csv", row.names = 1)
  # heart.meta <- read.table("Raw/Wang2020_HeartSC_GSE109816_metadata.txt", sep = "\t", header = 1)
  # rownames(heart.meta) <- heart.meta$Sample.name..9994.single.cells.
  heart.clust <- read.table("Raw/Wang2020_HeartSC_GSE109816_clustering.txt", sep = "\t", header = 1)
  rownames(heart.clust) <- heart.clust$ID
  
  # filter to common samples
  keep <- heart.clust$ID[which(heart.clust$ID %in% colnames(heart.sc))]
  heart.sc <- heart.sc[,keep]
  heart.clust <- heart.clust[keep,]
  
  # filter to left atrial samples
  heart.clust <- heart.clust[which(heart.clust$Condition == "LA"),] 
  
  # the group column encodes 2 different sample preparations: use "CM" for any CM cells, and "NCM" for all others
  remove1 <- which(heart.clust$Group == "CM" & heart.clust$CellType != "CM")
  remove2 <- which(heart.clust$Group == "NCM" & heart.clust$CellType == "CM")
  
  remove3 <- which(heart.clust$CellType == "MP") # remove macrophage, too
  heart.clust <- heart.clust[-c(remove1, remove2, remove3),]
  
  # reannotate
  heart.sc <- addENSID(heart.sc)
  
  # rpkm
  heart.sc <- rpkm(heart.sc)
  
  exp <- list()
    for(j in levels(as.factor(heart.clust$CellType))) {
      use <- heart.clust$ID[which(heart.clust$CellType == j)]
      exp[[j]] <- rowMeans(heart.sc[,use])
    }
  exp <- do.call("cbind", exp)
  exp <- as.data.frame(exp)
    
  # expression threshold
  exp <- exp[which(apply(exp, 1, max) > 1),]
    
  # return
  heart$SC <- exp  

## Rename columns
  heart[2:4] <- lapply(heart[2:4], function(x) {
    colnames(x) <- c("CM", "EC", "FB", "SMC")
    return(x)
  })
  

## Save
  save(heart, file = "Preprocessed/heart.rda")
  
## Estimate composition
  est.heart <- list()
  est.heart$CIB <- list()
  est.heart$DRS <- list()
  est.heart$dtangle <- list()
  
  for(j in names(heart[-1])) {
    est.heart$CIB[[j]] <- run.CIB(from.file = FALSE,
                              sigObject = heart[[j]],
                              mixString = "GTEx_Heart.txt")
    
    est.heart$DRS[[j]] <- run.DRS(heart$gtex$exp, heart[[j]])
    
    est.heart$DTA[[j]] <- run.dtangle(heart$gtex$exp, heart[[j]], alg = "diff", 0.01)
  }
  
  save(est.heart, file = "../Results/Revisions/CrossTissue/Heart Estimates.rda")
  
## Evaluate composition
  ## Oh, do please split into Aortic samples!
  ventricle <- grep("Ventricle", heart$gtex$meta$SMTSD)
  atrialAppendage <- grep("Atrial Appendage", heart$gtex$meta$SMTSD)
  
  ## GoF
  gof.heart <- list()
  for(j in names(est.heart)) { # loop through algorithms
    gof.heart[[j]] <- list()
    for (k in names(est.heart[[j]])) { # loop through signatures
      print(paste0(j, "_", k))
      gof.heart[[j]][[k]] <- write.gof.v2(measuredExp = heart$gtex$exp, 
                                          estimatedComp = est.heart[[j]][[k]], 
                                          signatureUsed = heart[[k]])
    }
  }  
  
  save(gof.heart, file = "../Results/Revisions/CrossTissue/Heart Goodness of Fit.rda")  # load("../Results/Revisions/CrossTissue/Heart Goodness of Fit.rda")
 
  plot.data <- lapply(gof.heart, function(x) { sapply(x, function(y) y$r) })
  plot.data <- lapply(plot.data, function(x) { 
    x <- as.data.frame(x)
    x$Type <- "AA"
    x$Type[ventricle] <- "LV"
    return(x)
    })
  plot.data <- melt(plot.data)
  colnames(plot.data) <- c("Region", "Signature", "Pearson", "Algorithm")
  levels(plot.data$Signature) <- c("F5 (Cultured / Bulk)", "EN (Cultured / Bulk)", "SC (Fresh / scRNAseq)")
  
  pdf(file = "../Results/Revisions/CrossTissue/Heart Goodness of Fit Violins (Pg1=AA, Pg2=LV).pdf", height = 2, width = 8)
  ggplot(plot.data[which(plot.data$Region == "AA"),], aes(x = Algorithm, fill = Signature, y = Pearson, colour = Signature)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.8, colour = "black") +
    geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
    # facet_wrap(~Region, ncol = 1) +
    # geom_violin() +
    scale_y_continuous(limits = c(0.45,0.82), expand = c(0,0)) +
    theme_bw() +
    labs(y = "Goodness of Fit (r)") +
              theme(panel.border = element_blank(), axis.line = element_line())
  
  ggplot(plot.data[which(plot.data$Region == "LV"),], aes(x = Algorithm, fill = Signature, y = Pearson, colour = Signature)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.8, colour = "black") +
    geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
    # facet_wrap(~Region, ncol = 1) +
    # geom_violin() +
    scale_y_continuous(limits = c(0.45,0.82), expand = c(0,0)) +
    theme_bw() +
    labs(y = "Goodness of Fit (r)") +
              theme(panel.border = element_blank(), axis.line = element_line())
  dev.off()
  
  ## Correlation
  plot.data <- est.heart
  plot.data <- lapply(plot.data, function(x) {
    x$ENCODE <- x$ENCODE[,c(2,3,4,1)]
    return(x)
  })
  
  plot.data <- lapply(plot.data, function(x) {
    x <- lapply(x, function(y) {
      colnames(y) <- c("CM", "EC", "FB", "SMC")
      return(y)
    })
    return(x)
  })
  
  x <- lapply(plot.data, function(x) {
    do.call("cbind", x)
  })
  x <- do.call("cbind", x)
  x <- as.data.frame(x)
  
  cor.atrial <- cor(x[atrialAppendage,], method = "s")
  cor.ventricle <- cor(x[ventricle,], method = "s")
  
  plot.cor <- function(celltype = "CM", corMx = cor.ventricle) {
    g <- grep(celltype, colnames(corMx))
    dat <- as.data.frame(corMx[g,g])
    # dat[upper.tri(dat)] <- which(upper.tri(dat))
    labels <- strsplit(colnames(dat), "\\.")
    dat$Group1 <- gsub(paste0(".", celltype), "", colnames(dat))
    # 
    # 
    # dat$Alg <- sapply(labels, "[", 1)
    # dat$Sig <- sapply(labels, "[", 2)
    # 
    dat <- melt(dat)
    dat$Group2 <- gsub(paste0(".", celltype), "", dat$variable)
    dat$Rho <- round(dat$value, 2)
    
    ggplot(dat, aes(x = Group1, y = Group2, fill = Rho, label = Rho)) +
      geom_tile(colour = "black") +
      geom_text(size = 2) +
      theme_bw() +
      scale_fill_viridis(limits = c(-1, 1)) +
      labs(title = celltype) +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
            panel.border = element_blank())
    
  }
  
  pdf(file = "../Results/Revisions/CrossTissue/Heart Proportion Spearman (Ventricle).pdf", height = 7, width = 9)
  pA <- plot.cor("CM", corMx = cor.ventricle)
  pB <- plot.cor("EC")
  pC <- plot.cor("SMC")
  pD <- plot.cor("FB")
  plot_grid(pA, pB, pC, pD, nrow = 2)
  dev.off()
  
  pdf(file = "../Results/Revisions/CrossTissue/Heart Proportion Spearman (Atrial Appendage).pdf", height = 7, width = 9)
  pA <- plot.cor("CM", corMx = cor.atrial)
  pB <- plot.cor("EC", corMx = cor.atrial)
  pC <- plot.cor("SMC", corMx = cor.atrial)
  pD <- plot.cor("FB", corMx = cor.atrial)
  plot_grid(pA, pB, pC, pD, nrow = 2)
  dev.off()
  
  ## Distribution of proportions
  plot.data <- est.heart
  plot.data <- lapply(plot.data, function(x) {
    x$ENCODE <- x$ENCODE[,c(2,3,4,1)]
    return(x)
  })
  
  plot.data <- lapply(plot.data, function(x) {
    x <- lapply(x, function(y) {
      colnames(y) <- c("CM", "EC", "FB", "SMC")
      y$Region <- heart$gtex$meta$SMTSD
      return(y)
    })
    return(x)
  })
  
  plot.data <- melt(plot.data)
  
  pdf(file = "../Results/Revisions/CrossTissue/Heart Proportion Distributions (Atrial Appendage).pdf", height = 3, width = 10)
  ggplot(plot.data[which(plot.data$Region == "Heart - Atrial Appendage"),], aes(x = variable, y = value, fill = L2, colour = L2, shape = L2)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.35, size = 1) +
    # geom_violin() +
    stat_summary(fun = median, colour = "black", position = position_dodge(width = 0.7)) +
    scale_shape_manual(values = c(21, 23, 24)) +
    facet_wrap(~L1) +
    theme_bw() +
    scale_y_continuous(limits = c(-0.01,1.01), expand = c(0,0)) +
    labs(y = "Estimated Proportion", x = "Celltype") +
    theme(legend.position = "right", axis.line.y = element_line(), legend.title = element_blank())
  dev.off()
  
  pdf(file = "../Results/Revisions/CrossTissue/Heart Proportion Distributions (Left Ventricle).pdf", height = 3, width = 10)
  ggplot(plot.data[which(plot.data$Region == "Heart - Left Ventricle"),], aes(x = variable, y = value, fill = L2, colour = L2, shape = L2)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.35, size = 1) +
    # geom_violin() +
    stat_summary(fun = median, colour = "black", position = position_dodge(width = 0.7)) +
    scale_shape_manual(values = c(21, 23, 24)) +
    facet_wrap(~L1) +
    theme_bw() +
    scale_y_continuous(limits = c(-0.01,1.01), expand = c(0,0)) +
    labs(y = "Estimated Proportion", x = "Celltype") +
    theme(legend.position = "right", axis.line.y = element_line(), legend.title = element_blank())
  dev.off()
  
################################################################################################################################ #
## Pancreas analyses ----
  
## Setup
pancreas <- list()
  
## Get GTEx data
pancreas$gtex <- extract.gtex.tissue("Pancreas", norm = "tpm")
write.CIB(data = pancreas$gtex$exp, dir = "Preprocessed/CIB/GTEx_Pancreas.txt")  
  
# ## Load signatures (from pancreas scme, used in script 7)
#   load("Preprocessed/pancreas.rda")
#   rm(sim.de)
#   
#   # add signatures to the list...
#   pancreas <- c(pancreas, sigsPancreas)

## Signature 1: single-cell data from cells freshly sorted from the human pancreas (Enge et al., Cell, 2019)
  ## Read in count-level data
    # get filenames
    files <- list.files(path = "Raw/Enge2017_GSE81547_RAW/")
    file.list <- paste0("Raw/Enge2017_GSE81547_RAW/", files)
    file.list <- as.list(file.list)
    names(file.list) <- gsub(".csv.gz", "", files)
    
    # read in RNA-seq from each file
    EN <- lapply(file.list, fread) # fread takes ~1m, requiring 3.9GB
    
    # convert to dataframe
    gene.names <- EN$GSM2171880_1000010011.A07$V1 # the column vector of gene names
    EN <- sapply(EN, function(x) x$V2)
    rownames(EN) <- gene.names
    EN <- as.data.frame(EN)
    colnames(EN) <- sapply(strsplit(colnames(EN), "_"), "[", 1)
    
  ## Read in metadata
    EN.meta <- read.table("Raw/Enge2017_GSE81547_series_matrix_edited.txt", sep= "\t")
    EN.meta <- as.data.frame(t(EN.meta))
    EN.meta <- EN.meta[-1,-1]
    colnames(EN.meta) <- c("SampleID", "Age", "Gender", "Celltype")
    EN.meta$Celltype <- gsub("inferred_cell_type: ", "", EN.meta$Celltype)
    
  ## Confirm matching order in expression data and metadata
    plot(match(colnames(EN), EN.meta$SampleID))
  
  ## Convert to EnsID
    EN <- addENSID(EN)
  
  ## Convert to TPM
    EN <- tpm(EN)
  
  ## Create signature
    EN.sig <- data.frame(alpha = rowMeans(EN[,which(EN.meta$Celltype == "alpha")]),
                         beta = rowMeans(EN[,which(EN.meta$Celltype == "beta")]))
    
    EN.sig <- EN.sig[which(apply(EN.sig, 1, max) > exp_thresh),]
  
  ## Output
    pancreas$EN <- EN.sig

## Signature 2: FACS-based isolation of adult alpha and beta cells (Blodgett et al, 2015)
  ## Read in preprocessed to TPM
  BL <- read.csv("Raw/Blodgett2015_TPM.csv")
  rownames(BL) <- BL$GeneName
  BL <- BL[,-1]
  
  ## Convert to EnsID
  BL <- addENSID(BL)
  
  ## Collect signature
  BL.sig <- data.frame(alpha = rowMeans(BL[,grep("Adult.Alpha", colnames(BL))]),
                       beta = rowMeans(BL[,grep("Adult.Beta", colnames(BL))]))
  
  BL.sig <- BL.sig[which(apply(BL.sig, 1, max) > exp_thresh),]
  pancreas$BL <- BL.sig

## Signatures 3+4: from Furuyama et al 2019
  ## Read in data
    fur <- read.table("Raw/Furuyama2019_GSE117454_RNA_seq_matrix.txt", sep = "\t", header = TRUE, row.names = 1)
    colnames(fur) <- gsub("X", "", colnames(fur))
  
  ## Read in metadata
    fur.meta <- read.table("Raw/Furuyama2019_GSE117454_series_matrix_edited.txt", sep = "\t")
    fur.meta <- as.data.frame(t(fur.meta))
    fur.meta <- fur.meta[-1,c(1,3)]
    colnames(fur.meta) <- c("SampleID", "Celltype")
    fur.meta$SampleID <- sapply(strsplit(fur.meta$SampleID, "\\["), "[", 2)
    fur.meta$SampleID <- gsub("]", "", fur.meta$SampleID)
    
  ## Confirm identical order in expression and metadata
    plot(match(colnames(fur), fur.meta$SampleID))
  
  ## Convert to EnsID
    fur <- addENSID(fur)
  
  ## Convert to TPM
    fur <- tpm(fur)
  
  ## Create signatures
    ## Signature 3: Freshly sorted
    FS.sig <-  data.frame(alpha = rowMeans(fur[,which(fur.meta$Celltype == "sorted α-cells")]),
                          beta = rowMeans(fur[,which(fur.meta$Celltype == "sorted β-cells")])) 
    FS.sig <- FS.sig[which(apply(FS.sig, 1, max) > exp_thresh),]
    pancreas$FS <- FS.sig
    
    ## Signature 4: cultured for 1-week and transduced with a GFP vector
    FG.sig <-  data.frame(alpha = rowMeans(fur[,which(fur.meta$Celltype == "αGFP")]),
                          beta = rowMeans(fur[,which(fur.meta$Celltype == "βGFP")])) 
    FG.sig <- FG.sig[which(apply(FG.sig, 1, max) > exp_thresh),]
    pancreas$FG <- FG.sig
    
## Save
  save(pancreas, file = "Preprocessed/pancreas.rda") # load("Preprocessed/pancreas.rda")

## Deconvolve
  est.pancreas <- list()
  est.pancreas$CIB <- list()
  est.pancreas$DRS <- list()
  est.pancreas$DTA <- list()
  
  for(j in names(pancreas[-1])) {
    print(paste("Starting", j, "at", Sys.time()))
    
    est.pancreas$CIB[[j]] <- run.CIB(from.file = FALSE,
                                  sigObject = pancreas[[j]],
                                  mixString = "GTEx_Pancreas.txt")
    
    est.pancreas$DRS[[j]] <- run.DRS(pancreas$gtex$exp, pancreas[[j]])
    
    est.pancreas$DTA[[j]] <- run.DTA(pancreas$gtex$exp, pancreas[[j]], alg = "diff", 0.01)
    
    save(est.pancreas, file = "../Results/Revisions/CrossTissue/Pancreas Estimates (in progress).rda")
  }
  
  save(est.pancreas, file = "../Results/Revisions/CrossTissue/Pancreas Estimates (Final).rda")
  

## GoF
  gof.pancreas <- list()
  for(j in names(est.pancreas)) { # loop through algorithms
    gof.pancreas[[j]] <- list()
    for (k in names(est.pancreas[[j]])) { # loop through signatures
      print(paste0(j, "_", k))
      gof.pancreas[[j]][[k]] <- write.gof.v2(measuredExp = pancreas$gtex$exp, 
                                       estimatedComp = est.pancreas[[j]][[k]], 
                                       signatureUsed = pancreas[[k]])
    }
  }  
  
  save(gof.pancreas, file = "../Results/Revisions/CrossTissue/Pancreas Goodness of Fit.rda")  # load("../Results/Revisions/CrossTissue/Pancreas Goodness of Fit.rda")
  
  plot.data <- lapply(gof.pancreas, function(x) { sapply(x, function(y) y$r) })
  plot.data <- melt(plot.data)
  colnames(plot.data) <- c("X", "Signature", "Pearson", "Algorithm")
  levels(plot.data$Signature) <- c("EN (Fresh / scRNAseq)", "BL (Fresh / Bulk RNAseq)", "FS (Fresh / Bulk RNAseq)", "FG (Cultured / Bulk RNAseq)")
  
  pdf(file = "../Results/Revisions/CrossTissue/Pancreas Goodness of Fit Violins.pdf", height = 2, width = 8)
  ggplot(plot.data, aes(x = Algorithm, fill = Signature, y = Pearson, colour = Signature)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.8, colour = "black") +
    geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
    # geom_violin() +
    scale_y_continuous(limits = c(0.4,NA), expand = c(0,0)) +
    theme_bw() +
    labs(y = "Goodness of Fit (r)") +
    theme(panel.border = element_blank(), axis.line = element_line())
  dev.off()
  
## Correlation
  plot.data <- est.pancreas
  
  plot.data <- lapply(plot.data, function(x) {
    do.call("cbind", x)
  })
  plot.data <- do.call("cbind", plot.data)
  plot.data <- as.data.frame(plot.data)
  plot.data <- cor(plot.data, method = "s")

  plot.cor <- function(celltype = "CM", corMx = plot.data) {
    g <- grep(celltype, colnames(corMx))
    dat <- as.data.frame(corMx[g,g])
    colnames(dat) <- gsub(paste0(".", celltype), "", colnames(dat))
    colnames(dat) <- rownames(dat) <- gsub("\\.", "/", colnames(dat))
    
    # dat[upper.tri(dat)] <- which(upper.tri(dat))
    
    # labels <- strsplit(colnames(dat), "\\.")
    dat$Group1 <- rownames(dat)
    # 
    # 
    # dat$Alg <- sapply(labels, "[", 1)
    # dat$Sig <- sapply(labels, "[", 2)
    # 
    dat <- melt(dat)
    # dat$Group2 <- gsub(paste0(".", celltype), "", dat$variable)
    dat$Rho <- round(dat$value, 2)
    
    ggplot(dat, aes(x = Group1, y = variable, fill = Rho, label = Rho)) +
      geom_tile(colour = "black") +
      geom_text(size = 2) +
      theme_bw() +
      # scale_fill_viridis(limits = c(-1, 1)) +
      scale_fill_carto_c(palette = "ArmyRose", limits = c(NA, NA)) +
      labs(title = celltype) +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
            panel.border = element_blank(), axis.ticks = element_blank())
    
  }
  
  pdf(file = "../Results/Revisions/CrossTissue/Pancreas Proportion Spearman.pdf", height = 4, width = 5)
  plot.cor("alpha") + labs(title = "")
  plot.cor("beta") # corMx identical to alpha as beta = 1 - alpha
  dev.off()
  
## Distribution of proportions
  plot.data <- melt(est.pancreas)
  plot.data$L2 <- factor(plot.data$L2, levels = c("EN", "BL", "FS", "FG"))
  levels(plot.data$L2) <- c("EN (Fresh / scRNAseq)", "BL (Fresh / Bulk RNAseq)", "FS (Fresh / Bulk RNAseq)", "FG (Cultured / RNAseq)")

  pdf(file = "../Results/Revisions/CrossTissue/Pancreas Proportion Distributions.pdf", height = 3, width = 8)
  ggplot(plot.data[which(plot.data$variable == "alpha"),], aes(x = L1, y = value, fill = L2, colour = L2, shape = L2)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), alpha = 0.35, size = 1) +
    # geom_violin() +
    stat_summary(fun = median, colour = "black", position = position_dodge(width = 0.7)) +
    scale_shape_manual(values = c(21, 23, 24, 25)) +
    # facet_wrap(~L1) +
    theme_bw() +
    scale_y_continuous(limits = c(-0.01,1.01), expand = c(0,0)) +
    labs(y = "Estimated Proportion of Alpha Cells", x = "Algorithm") +
    theme(legend.position = "right", axis.line.y = element_line(), legend.title = element_blank(),
          panel.border = element_blank(), axis.line = element_line())
  dev.off()
  
 
