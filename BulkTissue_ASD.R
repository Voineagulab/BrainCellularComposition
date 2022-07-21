## Here, we apply composition-estimation algorithms to ASD and CTL transcriptomes from Parikshak et al. 2016

################################################################################################################################ #
## Setup ----

## Start!
  rm(list = ls()); gc()

## Directory
  # root.dir <- "INSERT.YOUR.ROOT.DIRECTORY.PLEASE/"
  root.dir <- "/Volumes/Data1/PROJECTS/BrainCellularComposition/BrainCellularComposition_GitHub/"
  wd1 <- paste0(root.dir, "Data/") # working directory 1, for processing
  wd2 <- paste0(root.dir, "Results/BulkTissue/Autism/") # and for analysis
  setwd(wd1)
  
  

## Functions and packages
  source("../Scripts/Fun_Composition.R")
  source("../Scripts/Fun_Preprocessing.R")
  source("../Scripts/Fun_ASD.R")
  source("../Scripts/Fun_Parameters.R")

  # please download these files from the GitHub, in ./Files
  load("Preprocessed/exonicLength.rda") # exonic lengths of genes, used for gene length normalisation
  load("Preprocessed/geneInfo.rda") # information about gene symbols, ids, and type
  
## Load signatures
  load("../../../../Data/Preprocessed/Multibrain.rda")
  

################################################################################################################################ #
## Process expression data ----
  
## Data from Parikshak et al, 2016
  load("Raw/ASD_RNAseq_ExpressionData.Rdata")
  
## Remove vermis (cerebellar) samples
  pCounts <- datExpr.HTSC.unionexon[,-grep("vermis", colnames(datExpr.HTSC.unionexon))] # removing samples from the cerebellar vermins
  pMeta <- datMeta[-grep("vermis", rownames(datMeta)),] # removing samples from the cerebellar vermis
  
## Clean environment
  rm(datExpr.Cufflinks, datExpr.HTSC.unionexon, datExpr.HTSC.wholegene, datLB.gene, datMeta, datOK.gene)
  
## RPKM normalisation
  pRPKM <- rpkm(pCounts)
  
## Expression thresholding
  exp_thresh <- 1 # minimum rpkm cutoff
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
  
  
## Filter to the subset analysed in Parikshak et al.
  # get filter
  subset <- which(pMeta$ASD.vs.CTL..CTX) # these 106 samples are those which Parikshak et al used in their DE analyses; it excludes young and/or Dup15q samples

  # apply filter
  pRPKM <- pRPKM[,subset]
  pMeta <- pMeta[subset,]
  pCounts <- pCounts[,subset]
  
## Save
  # rda
  save(pCounts, pRPKM, pMeta, file = "Preprocessed/ASD.rda")
  
  # cibersort-compatible files
  write.CIB(data = pRPKM, dir = "Preprocessed/CIB/ASD.txt") 
  
## Shorthand
  asd <- which(pMeta$ASD.CTL == "ASD")
  ctl <- which(pMeta$ASD.CTL == "CTL")
  
################################################################################################################################ #
## Estimate composition ----

setwd(wd2)
  
## Use CIBERSORT and the MultiBrain signature
  est.asd <- run.CIB(from.file = FALSE,
                            sigObject = multibrain,
                            mixString = "ASD.txt")
  
## Calculate goodness of fit
  gof.asd <- write.gof.v2(measuredExp = pRPKM, estimatedComp = est.asd, signatureUsed = multibrain)
 
  # use $r 
  summary(gof.asd$r) # median of 0.71

  

## Save
  save(est.asd, gof.asd, file = "Composition Estimates.rda")

  
## Evaluate trends in composition
  # calculate the difference in composition between ASD and CTL samples
  comp.diff <- data.frame(p = apply(est.asd, 2, function(y) { wilcox.test(y ~ pMeta$ASD.CTL)$p.value }), 
                          sign = apply(est.asd, 2, function(y) { sign(mean(y[asd]) - mean(y[ctl])) }), 
                          magnitude = apply(est.asd, 2, function(y) { (mean(y[asd]) - mean(y[ctl])) }))
  
  # violin plot
  plot.data <- est.asd
  plot.data$disease <- pMeta$ASD.CTL
  plot.data <- melt(plot.data)
  colnames(plot.data)[2] <- "Celltype"
  levels(plot.data$Celltype) <- substr(levels(plot.data$Celltype), 1, 3)

  pdf(file = "Cell-type Proportion Distributions.pdf", height = 3, width = 4)
  ggplot(plot.data, aes(x = Celltype, y = value, fill = disease, colour = disease)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.8, colour = "black") +
    geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
    theme_bw() +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_fill_manual(values = asd.colours) +
    scale_colour_manual(values = asd.colours) +
    labs(y = paste0("Estimated Proportion")) +
    theme(panel.grid = invis, legend.title = invis, axis.title.x = invis, axis.line.y = element_line(),
          legend.position = c(0.88, 0.84), panel.border = invis, legend.background = element_rect(colour = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
  dev.off()
      
  
## Plot correlations of cell-types within CIB.IP
  plot.data <- as.data.frame(cor(est.asd, method = "p"))
  plot.data$ct <- rownames(plot.data)
  plot.data <- melt(plot.data)
  colnames(plot.data)[3] <- "Pearson"

  plot.data$ct <- factor(plot.data$ct, levels = colnames(est.asd))
  plot.data$variable <- factor(plot.data$variable, levels = colnames(est.asd))
  
  levels(plot.data$ct) <- substr(levels(plot.data$ct), 1, 3)
  levels(plot.data$variable) <- substr(levels(plot.data$variable), 1, 3)
  
  pdf(file = "Inter cell-type correlations.pdf", height = 4, width = 4.2)
  ggplot(plot.data, aes(x = ct, y = variable, fill = Pearson)) +
    geom_tile(colour = "black", width = 0.95, height = 0.95) +
    geom_text(aes(label = round(Pearson, 2)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1,1)) +
    theme_bw() +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title = invis)
  dev.off()


################################################################################################################################ #
## DESeq2 ----
    
## First, bind proportions to the metadata tables
  pMeta$Astrocytes <- est.asd$Astrocytes
  pMeta$Oligodendrocytes <- est.asd$Oligodendrocytes
  pMeta$Microglia <- est.asd$Microglia
    

## DESeq2 Setup and parameters
  DE.output <- list() # will contain DE output for all genes
  degs <- list() # will contain DE output for significant genes...
  p <- 0.05 # at adjusted p < 0.05
  runNames <- c("CD", "CI") # composition-dependent and -independent models, respectively
  formulae <- vector("list", length = length(runNames))
  names(formulae) <- runNames
  formulae$CD <- "~ Age + Sex + RIN + RegionID + BrainBank + SeqBatch + seqStatPC1 + seqStatPC2 + ASD.CTL"
  formulae$CI <- "~ Age + Sex + RIN + RegionID + BrainBank + SeqBatch + seqStatPC1 + seqStatPC2 + Astrocytes + Oligodendrocytes + Microglia + ASD.CTL" # also models for Ast + Oli + Mic


## Analysis loop.  
  for (j in 1:length(runNames))  {  
    # counter
    cat(paste0("Starting DESeq2 loop ", j, " of ", length(runNames)))
    
    # run DESeq2
    dds <- DESeqDataSetFromMatrix(countData = pCounts, colData = pMeta, design = as.formula(formulae[[j]])) 
    dds <- DESeq(dds) 
    res <- results(dds, contrast = c("ASD.CTL", "ASD", "CTL"), alpha = 0.05)
    DE.output[[j]] <- as.data.frame(res@listData); rownames(DE.output[[j]]) <- rownames(res)
    names(DE.output)[j] <- runNames[j]
    
    # collect degs if > 0 genes passing the FDR threshold
    if (length(which(DE.output[[j]]$padj <= p)) > 0) { 
      
      # genes passing pvalue threshold
      degs[[j]] <- DE.output[[j]][which(DE.output[[j]]$padj <= p),]
      
      # annotate direction
      degs[[j]]$Direction <- "-"
      degs[[j]]$Direction[which(degs[[j]]$log2FoldChange > 0)] <- "Up"
      degs[[j]]$Direction[which(degs[[j]]$log2FoldChange < 0)] <- "Down"
      
    } else degs[[j]] <- "None"
    
    # rename
    names(degs)[j] <- runNames[j]
    print(Sys.time())
  }

## Save
  save(DE.output, degs, file = "DESeq2.rda") 

## Save DEG list
  up <- lapply(degs, function(x) rownames(x)[which(x$Direction == "Up")])
  dn <- lapply(degs, function(x) rownames(x)[which(x$Direction == "Down")])
  
  ci.only <- list(up = up$CI[-which(up$CI %in% up$CD)],
                  dn = dn$CI[-which(dn$CI %in% dn$CD)])
  
  cd.only <- list(up = up$CD[-which(up$CD %in% up$CI)],
                  dn = dn$CD[-which(dn$CD %in% dn$CI)])


################################################################################################################################ #
## Characterise ----


## Gene ontology and pathway analyses using gprofiler
  # setup
  go <- list()
  library(gprofiler2)
  
  run.go <- function(list) {
    gost(list, organism = "hsapiens", ordered_query = FALSE,
         multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
         measure_underrepresentation = FALSE, evcodes = FALSE,
         user_threshold = 0.05, correction_method = "fdr", custom_bg = rownames(pRPKM),
         numeric_ns = "")$result
  }
  
  go$CD.up <- run.go(list = up$CD)
  go$CD.dn <- run.go(list = dn$CD)
  go$CI.up <- run.go(list = up$CI)
  go$CI.dn <- run.go(list = dn$CI)
  go$CI.only.up <- run.go(list = ci.only$up)
  go$CI.only.dn <- run.go(list = ci.only$dn)
  go$CD.only.up <- run.go(list = cd.only$up)
  go$CD.only.dn <- run.go(list = cd.only$dn)
  
  
  # filter list to relevant categories (note: all categories were used for multiple-testing correction)
  go <- lapply(go, function(x) { 
    x[which(x$source %in% c("GO:CC", "GO:MF", "GO:BP", "KEGG", "REAC", "WP", "HP")), -c(1,2,13,14)]
  })
  
  # output
  folder <- "GO"
  if (file.exists(folder)) {
    cat("The folder already exists")
  } else {
    dir.create(folder)
  }
  
  for (j in names(go)) write.csv(x = go[[j]], file = paste0("GO/", j, ".csv"), row.names = FALSE, quote = FALSE)
  save(go, file = "GO/GO.rda")
  
  
## Overlap to cell-type markers
  marks <- nTopFeatures(multibrain, n = 100, alg = "diff")  

  # overlap gene lists to cell-type markers
  marker.overlap <- list()
    for (k in names(table(marks$forCelltype))) {
      marker.overlap[[paste0(k, "_CD_Up")]] <- overEnrich(list1 = marks$EnsID[which(marks$forCelltype == k)],
                                                        list2 = up$CD,
                                                        backgroundList = rownames(pRPKM))
      
      marker.overlap[[paste0(k, "_CD_Dn")]] <- overEnrich(list1 = marks$EnsID[which(marks$forCelltype == k)],
                                                        list2 = dn$CD,
                                                        backgroundList = rownames(pRPKM))
      
      marker.overlap[[paste0(k, "_CI_Up")]] <- overEnrich(list1 = marks$EnsID[which(marks$forCelltype == k)],
                                                        list2 = up$CI,
                                                        backgroundList = rownames(pRPKM))
      
      marker.overlap[[paste0(k, "_CI_Dn")]] <- overEnrich(list1 = marks$EnsID[which(marks$forCelltype == k)],
                                                        list2 = dn$CI,
                                                        backgroundList = rownames(pRPKM))
      
      print(k)
      
    }
  
  marker.overlap <- do.call("rbind", marker.overlap)
  write.csv(marker.overlap, file = "DEG Overlap with Cell-type Markers.csv")
  
  