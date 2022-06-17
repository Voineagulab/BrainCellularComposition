## Here, we apply composition-estimation algorithms to ASD and CTL transcriptomes from Parikshak et al. 2016

################################################################################################################################ #
## Setup ----

## Generic
rm(list=ls())
options(stringsAsFactors = FALSE)

## Parameters
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/ASD_revised/")

## Functions and packages
source("../../../../Scripts/Fun_Composition.R")
source("../../../../Scripts/Fun_Preprocessing.R")
source("../../../../Scripts/Fun_ASD.R")

## Load data
load("../../../../Data/Preprocessed/ASD.rda")
asd <- which(pMeta$ASD.CTL == "ASD")
ctl <- which(pMeta$ASD.CTL == "CTL")


## Load signatures
  # load
  load("../../../../Data/Preprocessed/Signatures - Brain.rda")
  load("../../../../Data/Preprocessed/Signatures - Brain (Full Ct and Subct).rda")
  load("../../../../Data/Preprocessed/Multibrain.rda")
  
  # reprocess
  names(sigsFull) <- paste0(names(sigsFull), "z")
  sigsBrain <- c(sigsBrain, sigsFull)
  
## Other
  invis <- element_blank()
  
################################################################################################################################ #
## Basic data exploration ----

## Principal components analysis
pca <- princomp(log2(pRPKM + 0.5), cor = TRUE)
pca$eigen <- pca$sdev^2
pca$propVar <- round(pca$eigen / sum(pca$eigen) * 100, digits = 3)

################################################################################################################################ #
## Estimate composition ----

## Load from previous analyses
load("../../Parikshak_revised/Composition Estimates (Final).rda")

## Filter
est.asd <- lapply(est[1:4], function(x) {
  x <- x[-grep("\\.", names(x))] # removes stable deconvolutions
  
  lapply(x, function(y) {
    y[rownames(pMeta),]
  })
})


################################################################################################################################ #
## Evaluate using Goodness of Fit ----

## Setup
gof <- list()

## Run logged
for (j in names(est.asd)) {
  
  a <- est.asd[[j]] # shorter name
  gof[[j]] <- list()
  
  for (k in names(a)) {
    
    print(paste0(j, ":", k))
    
    if (substr(k, 3, 3) == "_") { 
      s <- sigsBrain[[substr(k, 1, 2)]]  
    } else {
      s <- sigsBrain[[k]]  
    }
    
    gof[[j]][[k]] <- write.gof(measuredExp = pRPKM,
                               signatureUsed = s,
                               estimatedComp = a[[k]],
                               log = TRUE)
    
  }
}


## Save
save(gof, file = "Goodness of Fit.rda") # load("Goodness of Fit.rda")

## Plot
  plot.data <- lapply(gof, function(x) {
    sapply(x, function(y) y$r)
  })
  
  plot.data <- melt(plot.data)
  
## Plot
  plot.data$Var2 <- gsub("_", "\nSubneurons", plot.data$Var2)
  plot.data$Var2 <- gsub("z", "\nAll subtypes", plot.data$Var2)
  
  pdf(file = "GoF, Boxplots.pdf", height = 3, width = 8)
  ggplot(plot.data, aes(x = Var2, fill = L1, y = value)) +
    geom_boxplot()
  dev.off()
  
  pdf(file = "GoF, Mean Lollipop.pdf", height = 3, width = 8)
  ggplot(plot.data, aes(x = Var2, fill = L1, y = value)) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.8), shape = 21) +
    stat_summary(fun = mean, geom = "col", colour = "black", position = position_dodge(width = 0.8), width = 0.01) +
    theme_bw()
  dev.off()
  
## Plot means
  plot.data <- lapply(gof, function(x) {
    x <- sapply(x, function(y) mean(y$r))
    x <- data.frame(Mean = x, Sig = names(x))
  })
  
   plot.data <- melt(plot.data)
   
   plot.data$Type <- "Standard 5ct"
   plot.data$Type[grep("_", plot.data$Sig)] <- "Exc/Inh"
   plot.data$Type[grep("z", plot.data$Sig)] <- "All Subtypes"
   plot.data$Type <- factor(plot.data$Type, levels = rev(levels(as.factor(plot.data$Type))))
   plot.data$Sig <- substr(plot.data$Sig, 1, 2)
  
  pdf(file = "GoF, Mean Heatmap.pdf", height = 5, width = 5)
  ggplot(plot.data, aes(x = Sig, y = L1, fill = value, label = round(value, 3))) +
    geom_tile() +
    facet_wrap(~Type, ncol = 1) +
    geom_text(size = 2.5) +
    scale_fill_viridis_b() +
    # stat_summary(fun = mean, geom = "col", colour = "black", position = position_dodge(width = 0.8), width = 0.01) +
    theme_bw() +
    theme(axis.title = invis)
  dev.off()
  

################################################################################################################################ #
## Trends in cell-type proportion estimates ----
  

## Stats 
comp.diff <- list()
  for (j in names(est.asd)) {
    comp.diff[[j]] <- lapply(est.asd[[j]], function(x) {
      
      p <- apply(x, 2, function(y) { wilcox.test(y ~ pMeta$ASD.CTL)$p.value })
      sign <- apply(x, 2, function(y) { sign(mean(y[asd]) - mean(y[ctl])) })
      mag <- apply(x, 2, function(y) { (mean(y[asd]) - mean(y[ctl])) })
      
      output <- data.frame(p = p, sign = sign, magnitude = mag)
      output$Outcome <- "-"
      output$Outcome[which(p < 0.05 & sign == 1)] <- "Increased"
      output$Outcome[which(p < 0.05 & sign == -1)] <- "Decreased"
      
      return(output)
    })
  }

save(comp.diff, file = "Confound In Ct Proportions, Raw List.rda")

## Regroup comp.diff by the type of deconvolution
x <- lapply(comp.diff, function(x) {
  x <- do.call("rbind", x)
  
  # add extra information
  annot <- strsplit(rownames(x), "\\.")
  x$Sig <- sapply(annot, "[", 1)
  x$Ct <- sapply(annot, "[", 2)
  
  return(x)
  
})

x <- do.call("rbind", x)
x$Alg <- sapply(strsplit(rownames(x), "\\."), "[", 1)

## Get frequencies!
  ## When using base ct
    freq <- x[which(substr(x$Sig, 3, 3) == ""),]
    freq <- table(freq$Outcome, freq$Ct)
    freq <- t(freq)
    write.csv(freq, file = "Confound In Ct Proportions.csv")
    
  ## When using subneuronal divisions
    freq <- x[which(substr(x$Sig, 3, 3) == "_"),]
    freq <- table(freq$Outcome, freq$Ct)
    freq <- t(freq)
    write.csv(freq, file = "Confound In Ct Proportions (Subneurons).csv")
    

## Plotting global and algorithm-specific systematic composition changes
  ## Algorithm-specific: CIB.CA
    # Plotting dataframes
    plot.data <- est.asd$CIB$MB2
    # colnames(plot.data) <- c("Neu", "Ast", "Oli", "Mic")
    plot.data$disease <- pMeta$ASD.CTL
    plot.data <- melt(plot.data)
    colnames(plot.data)[2] <- "Celltype"
    levels(plot.data$Celltype) <- substr(levels(plot.data$Celltype), 1, 3)
    # plot.data$Celltype <- factor(plot.data$Celltype, levels = c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia"))
  

    # Plot 
    pdf(file = "ASD Compositional Changes, CIB.MB2.pdf", height = 2.5, width = 3.5)
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
    
    # Stats
    apply(est.asd$CIB$MB2, 2, function(x) wilcox.test(x~pMeta$ASD.CTL)$p.value)
    
    apply(est.asd$CIB$MB2, 2, function(x) mean(x[asd]) - mean(x[ctl]))
    # Neurons       Astrocytes Oligodendrocytes        Microglia       Endothelia 
    # 0.0726313242     0.0002051972     0.2520234828     0.0032741091     0.1915044004 
    
    cor(est.asd$CIB$MB2)
      
  
## Plot correlations of cell-types within CIB.IP
  plot.data <- as.data.frame(cor(est.asd$CIB$MB2, method = "p"))
  plot.data$ct <- rownames(plot.data)
  plot.data <- melt(plot.data)
  colnames(plot.data)[3] <- "Pearson"

  plot.data$ct <- factor(plot.data$ct, levels = colnames(est.asd$CIB$MB2))
  plot.data$variable <- factor(plot.data$variable, levels = colnames(est.asd$CIB$MB2))
  
  levels(plot.data$ct) <- substr(levels(plot.data$ct), 1, 3)
  levels(plot.data$variable) <- substr(levels(plot.data$variable), 1, 3)
  
  pdf(file = "Inter cell-type correlations, CIB.MB2.pdf", height = 4, width = 4.2)
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
    
## First, bind neuronal proportions to the metadata tables
pMeta$Astrocytes <- est.asd$CIB$MB2$Astrocytes
pMeta$Oligodendrocytes <- est.asd$CIB$MB2$Oligodendrocytes
pMeta$Microglia <- est.asd$CIB$MB2$Microglia
    

## DESeq2 Setup and parameters
DE.output <- list()
degs <- list()
p <- 0.05
runNames <- c("CD", "CI")
formulae <- vector("list", length = length(runNames))
names(formulae) <- runNames
formulae$CD <- "~ Age + Sex + RIN + RegionID + BrainBank + SeqBatch + seqStatPC1 + seqStatPC2 + ASD.CTL"
formulae$CI <- paste0("~ Age + Sex + RIN + RegionID + BrainBank + SeqBatch + seqStatPC1 + seqStatPC2 + Astrocytes + Oligodendrocytes + Microglia + ASD.CTL")


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
  save(DE.output, file = "Temp DE.rda")
}

## Save

save(DE.output, degs, file = "DESeq2.rda") # load this for a quick restart: load(file = "DESeq2.rda")


## Save as supplementary table
supp <- data.frame(EnsID = rownames(pRPKM),
                     Symbol = pFeatures$Gene.Symbol,
                     CD.FC = DE.output$CD$log2FoldChange,
                     CD.direction = sign(DE.output$CD$log2FoldChange),
                     CD.fdr = DE.output$CD$padj,
                     CD.sig = DE.output$CD$padj < 0.05,
                     CI.FC = DE.output$CI$log2FoldChange,
                     CI.direction = sign(DE.output$CI$log2FoldChange),
                     CI.fdr = DE.output$CI$padj,
                     CI.sig = DE.output$CI$padj < 0.05)

  # more details on directoin
  supp$CD.direction <- factor(supp$CD.direction)
  levels(supp$CD.direction) <- c("Down", "Up")
  
  supp$CI.direction <- factor(supp$CI.direction)
  levels(supp$CI.direction) <- c("Down", "Up")
  
   # categorise 
  supp$fn <- supp$CI.sig == TRUE & supp$CD.sig == FALSE 
  supp$fp <- supp$CI.sig == FALSE & supp$CD.sig == TRUE 
  supp$concordant <- supp$CI.sig == TRUE & supp$CD.sig == TRUE 
  
  # save
  write.csv(supp, file = "ST_DE.csv", row.names = FALSE)

## Save DEG list
up <- lapply(degs, function(x) rownames(x)[which(x$Direction == "Up")])
dn <- lapply(degs, function(x) rownames(x)[which(x$Direction == "Down")])

ci.only <- list(up = up$CI[-which(up$CI %in% up$CD)],
                dn = dn$CI[-which(dn$CI %in% dn$CD)])

cd.only <- list(up = up$CD[-which(up$CD %in% up$CI)],
                dn = dn$CD[-which(dn$CD %in% dn$CI)])


################################################################################################################################ #
## Characterise ----

## Overlap
  # enrichment of overlaps
  oe <- list()
  oe$uu <- overEnrich(up$CD, up$CI, backgroundList = rownames(pCounts))
  oe$dd <- overEnrich(dn$CD, dn$CI, backgroundList = rownames(pCounts))
  oe$ud <- overEnrich(up$CD, dn$CI, backgroundList = rownames(pCounts))
  oe$du <- overEnrich(dn$CD, up$CI, backgroundList = rownames(pCounts))
  oe <- do.call("rbind", oe)

  # venn diagrame
  ## Venn diagrams
  # upregulated genes
  venn.up <- list(CI = up$CI, CD = up$CD)
  
  grid.newpage()
  venn.diagram(venn.up, filename = "Venn of CD vs CI.CIBMB2, upregulated.png",
               scaled = TRUE, fill = comp.colours,
               height = 350, width = 350, resolution = 200, ext.text = FALSE, imagetype = "png")
  
  # downregulated genes
  venn.down <- list(CI = dn$CI, CD = dn$CD)
  
  grid.newpage()
  venn.diagram(venn.down, filename = "Venn of CD vs CI.CIBMB2, downregulated.png",
               scaled = TRUE, fill = comp.colours,
               height = 350, width = 350, resolution = 200, ext.text = FALSE, imagetype = "png")
  
## Gene ontology and pathway analyses using gprofiler
  # setup
  go <- list()
  
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
  for (j in names(go)) write.csv(x = go[[j]], file = paste0("GO/", j, ".csv"), row.names = FALSE, quote = FALSE)
  save(go, file = "GO/GO.rda")
  
  
## Overlap to cell-type markers
  marks <- nTopFeatures(multibrain2, n = 100, alg = "diff")  

  # overlap gene lists to cell-type markers
  marker.overlap <- list()
  for(j in names(genes)) {
    for (k in names(table(marks$forCelltype))) {
      marker.overlap[[paste0(j, ".", k)]] <- overEnrich(list1 = marks$EnsID[which(marks$forCelltype == k)],
                                                        list2 = up$CD,
                                                        backgroundList = rownames(pRPKM))
      
    }
  }
  marker.overlap <- do.call("rbind", marker.overlap)
  write.csv(marker.overlap, file = "Overlap to markers.csv")
  
  
## Length not in Parikshak DEGs
  load("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Raw/Parikshak2016_DEGs.rda")
  
  length(which(!(ci.only$up %in% Parikshak2016_DEGs$Upregulated)))
  length(which(!(ci.only$up %in% Parikshak2016_DEGs$Downregulated)))
  
  length(which(!(ci.only$dn %in% Parikshak2016_DEGs$Upregulated)))
  length(which(!(ci.only$dn %in% Parikshak2016_DEGs$Downregulated)))

