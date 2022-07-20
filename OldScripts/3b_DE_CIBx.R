##################################################################################################################################
## Setup

rm(list = ls())
gc()
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/DE_simulations/CIBx/")

## Load functions and packages
source("../../../../Scripts/Fun_Preprocessing.R")
library(ochRe)
library(rcartocolor)

##################################################################################################################################
## Write data for the webtool


## Load simulations
load("../Sims (Final).rda")

## Save simulations in CIBERSORT format
  write.CIB(sims$up0.rep1$confound.pe.1.5, dir = "CIBx/up0.rep1.confound.pe.1.5.txt")
  write.CIB(sims$up0.rep1$confound.pe.exc.1.5, dir = "CIBx/up0.rep1.confound.pe.exc.1.5.txt")
  
  write.CIB(sims$up50.rep1$confound.pe.1.5, dir = "CIBx/up50.rep1.confound.pe.1.5.txt")
  write.CIB(sims$up50.rep1$confound.pe.exc.1.5, dir = "CIBx/up50.rep1.confound.pe.exc.1.5.txt")
  
  reup <- list()
  reup$up0.rep1.confound.pe.1.5 <- read.table("up0.rep1.confound.pe.1.5.txt", row.names = 1, header = TRUE)
  reup$up0.rep1.confound.pe.1.5.exc <- read.table("up0.rep1.confound.pe.exc.1.5.txt", row.names = 1, header = TRUE)
  reup$up50.rep1.confound.pe.1.5 <- read.table("up50.rep1.confound.pe.1.5.txt", row.names = 1, header = TRUE)
  reup$up50.rep1.confound.pe.1.5.exc <-read.table("up50.rep1.confound.pe.exc.1.5.txt", row.names = 1, header = TRUE)
  reup <- lapply(reup, function(x) {
    x <- apply(x, 2, function(y) y / sum(y) * 10^6)
    return(x)
  })
  
  for (j in names(reup)) write.CIB(reup[[j]], dir = paste0(j, ".txt"))
  
## Plot the proportions of cell-type abundance within these
  abundances <- list(up00 = sims$up0.rep1$meta.prop,
                     up10 = sims$up50.rep1$meta.prop)

  abundances <- do.call("rbind", abundances)  
  abundances <- as.data.frame(abundances)
  abundances$Confound <- substr(rownames(abundances), 3, 3)
  abundances$Group <- substr(rownames(abundances), 6, 6)
  abundances <- melt(abundances, id.vars = c("Confound", "Group"))
  abundances$Confound <- factor(abundances$Confound)F
  levels(abundances$Confound) <- c("Simulation 1:\nNo Change in Exc Proportion",
                                   "Simulation 2:\n+10% Exc Proportion in Group B")
  abundances$variable <- substr(abundances$variable, 1, 3)

  pdf(file = "Composition Violinplots.pdf", height= 3, width = 5)  
  ggplot(abundances, aes(x = variable, fill = Group, y = value)) +
    geom_violin(scale = "width") +
    facet_wrap(~Confound) +
    theme_bw() +
    labs(y = "True Abundance Within Mixture") +
    scale_fill_manual(values = c("darkorange1", "black")) +
    theme(axis.title.x = element_blank())
    dev.off()

## Write a list of 1k genes, including 500 markers, all perturbed genes, and negative controls
  ## First, get 500 markers
    load("../../../../Data/Preprocessed/SeuratObjects.rda")
    CA.cells <- as.data.frame(obj$CA@assays$RNA@counts)
    CA.meta <- obj$CA@meta.data
    
  ## Remomve
    abundances <- table(CA.meta$orig.celltype)
    keep <- names(abundances)[which(abundances > 200)]
    keep <- which(CA.meta$orig.celltype %in% keep)
    CA.cells <- CA.cells[,keep]
    CA.meta <- CA.meta[keep,]
    
  ## Merge celltypes
    y <- x <- as.character(CA.meta$orig.celltype)
    y[grep("^Ex", x)] <- "Excitatory"
    y[grep("^In", x)] <- "Inhibitory"
    y[grep("^Astro", x)] <- "Astrocytes"
    y[grep("^Oli", x)] <- "Oligodendrocytes"
    y[grep("^OPC", x)] <- "OPC"
    
    CA.meta$merged <- y
  
    # set celltype label to the colnames of the expression dataframe
    colnames(CA.cells) <- CA.meta$merged
    
    # filter to highly expressed genes
    min.n <- min(table(colnames(CA.cells)))
    keep <- which(rowSums(CA.cells > 10) > min.n) # greater than 1 count in at least 230 samples (i.e. the number of samples in the smallest group)
    
    CA.cells <- CA.cells[keep,]
    
    ## Create new HCA signature for marker generation
      CA.sig <- list()
      
      for(k in levels(as.factor(y))) {
        use <- which(colnames(CA.cells) == k)
        
        temp <- apply(CA.cells[,use], 2, function(x) {
          lib.size <- 10^6 / sum(x)
          x <- x * lib.size
          return(x)
        })
        
        CA.sig[[k]] <- rowMeans(temp)
      }
      CA.sig <- as.data.frame(do.call("cbind", CA.sig))
      
    ## Clean workspace
      rm(obj)
      
    ## Finally, get markers  
      marks <- nTopFeatures(signature = CA.sig, n = 100, alg = "diff")
      
  ## Get perturbed genes
    x <- sims$dn195.rep1$meta.exp$confound.pe.1.1$genes
    
  ## Get 400 random other genes (because the 100 exc already overlap...)
    use <- unique(c(marks$EnsID, x$EnsID))
    
    use <- c(use, sample(rownames(sims$dn195.rep1$confound.p[-which(rownames(sims$dn195.rep1$confound.p) %in% use),]), size = 400, replace = FALSE))
  
  ## Write
    write.table(use, file = "Genelist.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  ## Save signature
    write.CIB(CA.sig, dir = "HCA_signature.txt")
    CA.sig <- read.table("HCA_signature.txt", header = TRUE, row.names = 1)
    
    
## Oh, and create inhibitory-specific perturbations 
  ct <- c("Astrocytes", "Excitatory", "Inhibitory", "Oligodendrocytes", "OPC") # which celltypes appear in the mixture
  # change.ct <- "Inhibitory" # one of the cts to perturb in both abundance and expression
  nPerGroup <- 50  # number of mixtures to make for each group
    
  ## Parameters for perturbing celltype proportions
    ## On the range of proportions to simulate
    # each level defines a different simulation. [1] and [2] define the min and max proportion for group "alpha", while [3] and [4] do so for group "beta"
    props <- list() 
    props$up0 <- data.frame(c(200, 300, 200, 300))
    props$up50 <- data.frame(c(200, 300, 250, 350))
    colnames(props$up0) <- colnames(props$up50) <- "Excitatory"
    
    
  ## Parameters for perturbing gene expression
    geneList <- list()
    
    # get marker genes
    marks.inh <- nTopFeatures(signature = CA.sig, n = 100, alg = "diff")
    marks.inh <- geneList$Markers[which(geneList$Markers$forCelltype %in% "Inhibitory"),]
    
    # half will be upregulated, and half downregulated
    which.up <- sample(x = 1:100, replace = FALSE, size = 50)
    
    genesUp <- rbind(geneList$Markers[which.up,],
                     data.frame(EnsID = names(rand.up), forCelltype = "Random"))
    
    genesDown <- rbind(geneList$Markers[-which.up,],
                       data.frame(EnsID = names(rand.dn), forCelltype = "Random"))
    
  ## Simulate
    # use an existsing  simulations as a base
    sims.inh$up0 <- sims$up0.rep2[1:3]
    sims.inh$up50 <- sims$up50.rep2[1:3]
    
    # layer on perturbed expression in the inhibitory component of expression
    sims.inh$up0 <- confound.expression(sims.inh$up0, genesUp = genesUp, replace.starting.level = FALSE, genesDown = genesDown, absoluteFC = 1.5, noise = 0, id = "confound.pe.inh.1.5", ct.specific = "Inhibitory")
    sims.inh$up50 <- confound.expression(sims.inh$up50, genesUp = genesUp, replace.starting.level = FALSE, genesDown = genesDown, absoluteFC = 1.5, noise = 0, id = "confound.pe.inh.1.5", ct.specific = "Inhibitory")
    
  ## Convert to cpm
    sims.inh <- lapply(sims.inh, function(x) {
      x$confound.p.cpm <- apply(x$confound.p, 2, function(x) x / sum(x) * 10^6)
      x$confound.pe.inh.1.5.cpm <- apply(x$confound.pe.inh.1.5, 2, function(x) x / sum(x) * 10^6)
      return(x)
    })
    
  ## Save
    write.CIB(sims.inh$up0$confound.pe.inh.1.5.cpm, dir = "up0.confound.pe.inh.1.5.txt")
    write.CIB(sims.inh$up50$confound.pe.inh.1.5.cpm, dir = "up50.confound.pe.inh.1.5.txt")
    
    save(sims.inh, file = "Inhibitory Simulations.rda")

    
## Abundance plot
    pdf(file = "Abundance Boxplots.pdf", height = 2.5, width = 3.5)
    
    # for the confounded simulation
    plot.data <- sims.inh$up50$meta.prop
    plot.data$Group <- rep(c("Group A", "Group B"), each = 50)
    plot.data <- melt(plot.data)
    colnames(plot.data)[2] <- "Celltype"
    levels(plot.data$Celltype) <- substr(levels(plot.data$Celltype), 1, 3)
    # plot.data$Celltype <- factor(plot.data$Celltype, levels = c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia"))
    
    
    ggplot(plot.data, aes(x = Celltype, y = value, fill = Group, colour = Group)) +
      geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.8, colour = "black") +
      geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
      theme_bw() +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      scale_fill_manual(values = c("darkorange1", "dodgerblue")) +
      scale_colour_manual(values = c("darkorange1", "dodgerblue")) +
      labs(y = paste0("Abundance In Mixture")) +
      theme(panel.grid = invis, legend.title = invis, axis.title.x = invis, axis.line.y = element_line(),
            legend.position = c(0.88, 0.84), panel.border = invis, legend.background = element_rect(colour = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
    
    # for the non-confounded simulation
    plot.data <- sims.inh$up0$meta.prop
    plot.data$Group <- rep(c("Group A", "Group B"), each = 50)
    plot.data <- melt(plot.data)
    colnames(plot.data)[2] <- "Celltype"
    levels(plot.data$Celltype) <- substr(levels(plot.data$Celltype), 1, 3)
    # plot.data$Celltype <- factor(plot.data$Celltype, levels = c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia"))
    
    
    ggplot(plot.data, aes(x = Celltype, y = value, fill = Group, colour = Group)) +
      geom_violin(scale = "width", position = position_dodge(width = 0.8), width = 0.8, colour = "black") +
      geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
      theme_bw() +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      scale_fill_manual(values = c("darkorange1", "dodgerblue")) +
      scale_colour_manual(values = c("darkorange1", "dodgerblue")) +
      labs(y = paste0("Abundance In Mixture")) +
      theme(panel.grid = invis, legend.title = invis, axis.title.x = invis, axis.line.y = element_line(),
            legend.position = c(0.88, 0.84), panel.border = invis, legend.background = element_rect(colour = "transparent"), legend.box.background = element_rect(fill = "transparent")) 
    dev.off()
    
##################################################################################################################################
## Analyse output
  
  
## First, here's the README on CIBx's output!!

# --------------------------------------------------------------------------------
# CIBERSORTx High Resolution Mode - OUTPUT FILES
# --------------------------------------------------------------------------------
# 
# MAIN OUTPUT FILES
# --------------------------------------------------------------------------------
# The main result of CIBERSORTx High Resolution is a set of .txt files and
# heatmaps for each individual cell type, showing the cell-type specific
# expression of individual genes at the sample level.
# 
# The name of these files contain the name of the cell type, followed by the
# window size used for the job (e.g. CIBERSORTxHiRes_job1_Mastcells_Window40.txt).
# 
# The "1" values in the expression matrix txt files are genes with insufficient
# evidence of expression (these genes are either not expressed or have inadequate
# statistical power to be imputed).
# 
# The NA values are genes that have inadequate statistical power to be imputed.
# 
# If you have provided a gene subset list, these file will have all the genes in
# the list that are found in the original mixture file given as input. If some
# genes are still missing, this could be due to different annotations or gene
# symbols between the gene subset list and the mixture file.
# 
# ADDITIONAL OUTPUT FILES
# --------------------------------------------------------------------------------
# CIBERSORTx High Resolution Mode runs both CIBERSORTx Fractions Mode and Group
# Mode, and a set of output files is generated from each analysis.


## First, load genesets
  # the 1k genes submitted to CIBx
  genes.1k <- read.table(file = "Genelist.txt")
  genes.1k <- as.character(genes.1k$V1)  
  
  # perturbed genes
  load("../Meta/Perturbed Genes (Final).rda")
  names(exc.dn) <- names(exc.up) <- names(rand.dn) <- names(rand.up) <- all.genes
  exc.dn <- which(exc.dn)
  exc.up <- which(exc.up)
  rand.dn <- which(rand.dn)
  rand.up <- which(rand.up)
  
## Next, state your celltypes
  ct <- c("Excitatory", "Inhibitory", "Astrocytes", "Oligodendrocytes", "OPC")
    
## Read in
  ## A function
    read.CIBx <- function(job.number) { # job number is supplied by CIBx's webtool
      x <- list()
      
      path <- paste0("CIBERSORTx_Job", job.number, "_output/CIBERSORTxHiRes_Job", job.number, "_")
      for (y in ct) {
        x[[y]] <- read.table(paste0(path, y, "_Window20.txt"), sep = "\t", header = TRUE, row.names = 1)
      }
       
      return(x)
    }
  
  ## Apply function
    res <- list()
    
    res$zeroProp.inExc <- read.CIBx(49)
    res$up10Prop.inExc <- read.CIBx(51)
    res$zeroProp.inAll <- read.CIBx(50)
    res$up10Prop.inAll <- read.CIBx(52)
    

## QC these
  qc.CIBx <- function(x) { # where x is the output of read.CIBx(n)
    x <- lapply(x, function(y) {
      y <- y[which(!(is.na(y[,1]))),] # remove any gene where expression is NA in the first column (which means NA in all columns)
    })
    
    return(x)
  }
  
  res <- lapply(res, qc.CIBx)
  
## Run DE
  ## Setup
    group.info <- data.frame(group = c(rep(1, 50), rep(0, 50)))
    
  ## Function
    lm.CIBx <- function(x) {
      x <- lapply(x, function(y) {
        # convert to log
        log <- t(log2(y + 0.5))
        
        # run lm
        mod <- lm(log ~ group, data = group.info)
        sum <- summary(mod) 
        pvals <- sapply(sum, function(y)  { y$coefficients["group", "Pr(>|t|)"] } )
        log2fc <- sapply(sum, function(y)  { y$coefficients["group", "Estimate"] } )
        z <- data.frame(pvals = pvals, log2fc = log2fc)
        rownames(z) <- sapply(strsplit(rownames(z), " "), "[", 2)
        z$padj <- p.adjust(z$pvals, n = length(all.genes)) # please note the n here, adjusting for the full transcriptome
        
        # categorise genes
        z$Cat <- "Neg/Marker"
        z$Cat[which(rownames(z) %in% names(exc.dn) )] <- "exc.dn"
        z$Cat[which(rownames(z) %in% names(exc.up) )] <- "exc.up"
        z$Cat[which(rownames(z) %in% names(rand.dn) )] <- "rand.dn"
        z$Cat[which(rownames(z) %in% names(rand.up) )] <- "rand.up"
      
        return(z)
      })
      return(x)
    }
  
  ## Apply function
    de <- lapply(res, lm.CIBx)
    
  ## Analyse trends in DE
    ## Function
    success.CIBx <- function(x) {
      x <- lapply(x, function(y) {
        # check for consistency in expected and observed fold-change
        y$fc.sign <- sign(y$log2fc)
        y$fc.dir <- sapply(strsplit(y$Cat, "\\."), "[", 2)
        y$fc.cons <- FALSE
        y$fc.cons[is.na(y$fc.dir)] <- TRUE # NAs occur for genes with simulated DE. these can be TRUE, because there is no expected fc direction
        y$fc.cons[y$fc.dir == "up" & y$fc.sign == 1] <- TRUE
        y$fc.cons[y$fc.dir == "dn" & y$fc.sign == -1] <- TRUE
        
        # convert to factor
        y$Cat <- factor(y$Cat, levels = c("exc.dn", "exc.up", "Neg/Marker", "rand.dn", "rand.up"))
        
        # extract out hits with padj < 0.05
        sig <- y[which(y$padj < 0.05 & y$fc.cons),] # significant p-value
       
        # hit rates!
        total.genes <- table(y$Cat)
        hit.genes <- table(sig$Cat)
      
        hit.rate <- hit.genes / total.genes
        
       # return
        return(hit.rate)
        # return(hit.genes)
      })
        return(x)
    }
    
    hit.rates <- lapply(de, success.CIBx)
    hit.rates <- lapply(hit.rates, function(x) do.call("rbind", x))

        
## A final summary of each dataset...
  ## Function
    summarise.expression <- lapply(res, function(x) {
      # get numbers of expressed genes
      n.exp <- lapply(x, function(y) {
        
        
        n.exp.genes <- 
        n.exp.upmarks <- 
        n.exp.uprand <- length(which(c(names(rand.dn)) %in% rownames(y)))
        n.exp.dnmarks <- length(which(c(names(exc.dn), names(exc.up)) %in% rownames(y)))
        n.exp.dnrand <- length(which(c(names(rand.dn), names(rand.up)) %in% rownames(y)))
        n.exp.marks <- length(which(c(names(exc.dn), names(exc.up)) %in% rownames(y)))
        n.exp.rand <- length(which(c(names(rand.dn), names(rand.up)) %in% rownames(y)))
        
        output <- list(n.total.genes = 1000,
                       n.exp.genes = nrow(y),
                       n.exp.marks = length(which(c(names(exc.dn), names(exc.up)) %in% rownames(y))),
                       n.exp.rand = length(which(c(names(rand.dn), names(rand.up)) %in% rownames(y))),
                       n.exp.upmarks = length(which(c(names(exc.up)) %in% rownames(y))),
                       n.exp.dnmarks = length(which(c(names(exc.dn)) %in% rownames(y))),
                       n.exp.uprand = length(which(c(names(rand.up)) %in% rownames(y))),
                       n.exp.dnrand = length(which(c(names(rand.dn)) %in% rownames(y))))
        
        return(output)
      })
      
      # run de
      de <- lm.CIBx(x)
      
      # analyse de
      hits <- lapply(de, function(y) {
        y$fc.sign <- sign(y$log2fc)
        y$fc.dir <- sapply(strsplit(y$Cat, "\\."), "[", 2)
        y$fc.cons <- FALSE
        y$fc.cons[is.na(y$fc.dir)] <- TRUE # NAs occur for genes with simulated DE. these can be TRUE, because there is no expected fc direction
        y$fc.cons[y$fc.dir == "up" & y$fc.sign == 1] <- TRUE
        y$fc.cons[y$fc.dir == "dn" & y$fc.sign == -1] <- TRUE
        
        # convert to factor
        y$Cat <- factor(y$Cat, levels = c("exc.dn", "exc.up", "Neg/Marker", "rand.dn", "rand.up"))
        
        # extract out hits with padj < 0.05
        sig <- y[which(y$padj < 0.05 & y$fc.cons),] # significant p-value
        
        # hit rates!
        n.de.marks <- 
        n.de.rand <- 
        
        output <- list(n.de.marks = length(grep("exc", sig$Cat)),
                       n.de.rand = length(grep("rand", sig$Cat)),
                       n.de.upmarks = length(grep("exc.up", sig$Cat)),
                       n.de.dnmarks = length(grep("exc.dn", sig$Cat)),
                       n.de.uprand = length(grep("rand.up", sig$Cat)),
                       n.de.dnrand = length(grep("rand.dn", sig$Cat)))
        # return
        return(output)
      })
      
      # combine all information
      output <- rbind(do.call("cbind", n.exp),
                      do.call("cbind", hits))
        
      return(output)
    })
    
  save(summarise.expression, file = "Number of Genes Expressed and DE.rda")
  for(j in names(summarise.expression)) write.csv(summarise.expression[[j]], file = paste0("Number of Genes Expressed and DE (", j, ").csv"))
   
  
    
##################################################################################################################################
## Plot output
  

## Reorganise dataframe
  # create
  plot.data <- melt(hit.rates)
    
  # convert NaN to zero
  plot.data$value[which(is.na(plot.data$value))] <- 0
    
  ## Move information to new columns
    # on proportion confounds
    plot.data$Proportion <- "No Composition Change" 
    plot.data$Proportion[grep("up10", plot.data$L1)] <- "+10% Exc"
    plot.data$Proportion <- factor(plot.data$Proportion, levels = c("No Composition Change", "+10% Exc"))
    
    # on cell-type in which expression is confounded
    plot.data$Celltypes <- "DE in all ct"
    plot.data$Celltypes[grep("inExc", plot.data$L1)] <- "DE in Exc Only"
    
  ## Rename
    plot.data$Var1 <- substr(plot.data$Var1, 1, 3)
    # plot.data$Var1 <- paste0("CIBx-isolated ", plot.data$Var1)
    
    ## Separate dataframe for negative controls
    plot.data2 <- plot.data[which(plot.data$Var2 == "Neg/Marker"),]
    ct.colours <- c("#009392", "#eeb479", "#e88471", "#e9e29c", "#fdfbe4")
    names(ct.colours) <- substr(ct, 1, 3)
    
    ## Final processing of position controls' plotting data
    plot.data <- plot.data[which(plot.data$Var2 != "Neg/Marker"),]
    plot.data$Var2 <- factor(plot.data$Var2)
    levels(plot.data$Var2) <- c("Exc Marker | Downregulated", "Exc Marker | Upregulated", "Non-marker | Downregulated", "Non-marker | Upregulated")
    marker.cols <- brewer_pal(palette = "Paired")(8)[c(1,2,7,8)]
    
    
## Plot non-specific DE
    pdf(file = "Plots - AllCt, TPR.pdf", height = 3.5, width = 5)
    ggplot(plot.data[which(plot.data$Celltypes == "DE in all ct"),], aes(x = Var1, colour = Var2, fill = Var2, y = value)) +
      stat_summary(fun = function(x) {x}, geom = "point", position = position_dodge(width = 0.7), size = 2) +
      geom_col(position = position_dodge(width = 0.7), width = 0.1, show.legend = FALSE) +
      theme_bw() +
      facet_wrap(~Proportion, ncol = 2) +
      scale_colour_manual(values = marker.cols) +
      scale_fill_manual(values = marker.cols) +
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(nrow = 2), colour = guide_legend(nrow = 2)) +
      labs(y = "% of Perturbed Genes\nReaching Significance", x = "CIBx-isolated Component of Mixture Expression") 
    dev.off()
    
    pdf(file = "Plots - AllCt, FPR.pdf", height = 3.5, width = 3)
    ggplot(plot.data2[which(plot.data2$Celltypes == "DE in all ct"),], aes(x = Var1, y = value)) +
      geom_col(position = position_dodge(width = 0.7), colour = "black", fill = "black", width = 0.1, show.legend = FALSE) +
      stat_summary(fun = function(x) {x}, geom = "point", colour = "black", position = position_dodge(width = 0.7),  size = 3) +
      theme_bw() +
      facet_wrap(~Proportion, ncol = 1) +
      # scale_y_continuous(expand = c(0,0), limits = c(-0.01, )) +
      # theme(panel.border = element_blank(), axis.line.y = element_line()) +
      guides(colour = guide_legend(title = "Ct-specific expression"), fill = guide_legend(title = "Ct-specific expression")) +
      labs(y = "% of Non-Perturbed Genes\nReaching Significance", x = "CIBx-isolated Component of\nMixture Expression") 
    dev.off()
    
## Plot Exc-specific DE
     pdf(file = "Plots - Exc, TPR.pdf", height = 3.5, width = 5)
    ggplot(plot.data[which(plot.data$Celltypes == "DE in Exc Only"),], aes(x = Var1, colour = Var2, fill = Var2, y = value)) +
      stat_summary(fun = function(x) {x}, geom = "point", position = position_dodge(width = 0.7), size = 2) +
      geom_col(position = position_dodge(width = 0.7), width = 0.1, show.legend = FALSE) +
      theme_bw() +
      facet_wrap(~Proportion, ncol = 2) +
      scale_colour_manual(values = marker.cols) +
      scale_fill_manual(values = marker.cols) +
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(nrow = 2), colour = guide_legend(nrow = 2)) +
      labs(y = "% of Perturbed Genes\nReaching Significance", x = "CIBx-isolated Component of Mixture Expression") 
    dev.off()
    
    pdf(file = "Plots - Exc, FPR.pdf", height = 3.5, width = 3)
    ggplot(plot.data2[which(plot.data2$Celltypes == "DE in Exc Only"),], aes(x = Var1, y = value)) +
      geom_col(position = position_dodge(width = 0.7), colour = "black", fill = "black", width = 0.1, show.legend = FALSE) +
      stat_summary(fun = function(x) {x}, geom = "point", colour = "black", position = position_dodge(width = 0.7),  size = 3) +
      theme_bw() +
      facet_wrap(~Proportion, ncol = 1) +
      # scale_y_continuous(expand = c(0,0), limits = c(-0.01, )) +
      # theme(panel.border = element_blank(), axis.line.y = element_line()) +
      guides(colour = guide_legend(title = "Ct-specific expression"), fill = guide_legend(title = "Ct-specific expression")) +
      labs(y = "% of Non-Perturbed Genes\nReaching Significance", x = "CIBx-isolated Component of\nMixture Expression") 
    dev.off()
    
#     
#     pdf(file = "Plots - True Positive Rate (V1 Exc-specific DE).pdf", height = 3, width = 8)
#     ggplot(plot.data[which(plot.data$Celltypes == "DE in Exc Only"),], aes(x = Var1, colour = Var2, fill = Var2, y = value)) +
#       stat_summary(fun = function(x) {x}, geom = "point", position = position_dodge(width = 0.7), size = 2) +
#       geom_col(position = position_dodge(width = 0.7), width = 0.1, show.legend = FALSE) +
#       theme_bw() +
#       facet_wrap(~Proportion) +
#       scale_colour_manual(values = ochre_palettes$nolan_ned[1:4]) +
#       scale_fill_manual(values = ochre_palettes$nolan_ned[1:4]) +
#       theme_bw() +
#       labs(y = "True Positive Rate", x = "Cell-type specific expression") +
#       theme(legend.title = element_blank())
#     dev.off()
#     
#     pdf(file = "Plots - True Positive Rate (V2 Focused on Exc).pdf", height = 3, width = 8)
#     ggplot(plot.data[which(plot.data$Var1 == "Exc"),], aes(x = Proportion, colour = Var2, fill = Var2, y = value)) +
#       stat_summary(fun = function(x) {x}, geom = "point", position = position_dodge(width = 0.7), size = 3) +
#       geom_col(position = position_dodge(width = 0.7), width = 0.1, show.legend = FALSE) +
#       theme_bw() +
#       facet_wrap(~Celltypes) +
#       scale_colour_manual(values = ochre_palettes$nolan_ned[1:4]) +
#       scale_fill_manual(values = ochre_palettes$nolan_ned[1:4]) +
#       theme_bw() +
#       labs(y = "True Positive Rate", x = "Cell-type specific expression") +
#       theme(legend.title = element_blank())
#     dev.off()
#     
# ## Separately plot false-positives (i.e., those in the Neg/Marker category)
#   ct.colours <- c("#009392", "#eeb479", "#e88471", "#e9e29c", "#fdfbe4")
#   names(ct.colours) <- substr(ct, 1, 3)
#   
#   pdf(file = "Plots - False Positive Rate.pdf", height = 3, width = 8)
#   ggplot(plot.data2, aes(x = Proportion, colour = Var1, fill = Var1, y = value)) +
#     geom_col(position = position_dodge(width = 0.7), width = 0.1, show.legend = FALSE) +
#     stat_summary(fun = function(x) {x}, geom = "point", position = position_dodge(width = 0.7), size = 3) +
#     theme_bw() +
#     facet_wrap(~Celltypes) +
#     scale_colour_manual(values = ochre_palettes$lorikeet[1:5]) +
#     scale_fill_manual(values = ochre_palettes$lorikeet[1:5]) +
#     # scale_y_continuous(expand = c(0,0), limits = c(-0.01, )) +
#     theme(panel.border = element_blank(), axis.line.y = element_line()) +
#     guides(colour = guide_legend(title = "Ct-specific expression"), fill = guide_legend(title = "Ct-specific expression")) +
#     labs(y = "False Positive Rate", x = "Simulation Settings")
#   dev.off()
  

##################################################################################################################################
## Repeat analyses in the inhibitory-specific output
  

## Load the simulations
  load("Inhibitory Simulations.rda")
  inh.genes <- sims.inh$up0$meta.exp$confound.pe.inh.1.5$genes
  
  inh.inh.dn <- inh.genes$EnsID[which(inh.genes$forCelltype == "Inhibitory" & inh.genes$Direction == "Down")]
  inh.inh.up <- inh.genes$EnsID[which(inh.genes$forCelltype == "Inhibitory" & inh.genes$Direction == "Up")]
  inh.rand.dn <- inh.genes$EnsID[which(inh.genes$forCelltype == "Random" & inh.genes$Direction == "Down")]
  inh.rand.up <- inh.genes$EnsID[which(inh.genes$forCelltype == "Random" & inh.genes$Direction == "Up")]
  
  
  
## Process data
  inh <- list()
  inh$zeroProp <- read.CIBx(55)
  inh$up10Prop <- read.CIBx(56)
  inh <- lapply(inh, qc.CIBx)
  
## Run DE
  ## Function
    lm.CIBx.2 <- function(x) { # this only changes the section on categorising genes
      x <- lapply(x, function(y) {
        # convert to log
        log <- t(log2(y + 0.5))
        
        # run lm
        mod <- lm(log ~ group, data = group.info)
        sum <- summary(mod) 
        pvals <- sapply(sum, function(y)  { y$coefficients["group", "Pr(>|t|)"] } )
        log2fc <- sapply(sum, function(y)  { y$coefficients["group", "Estimate"] } )
        z <- data.frame(pvals = pvals, log2fc = log2fc)
        rownames(z) <- sapply(strsplit(rownames(z), " "), "[", 2)
        z$padj <- p.adjust(z$pvals, n = length(all.genes)) # please note the n here, adjusting for the full transcriptome
        
        # categorise genes
        z$Cat <- "Neg/Marker"
        z$Cat[which(rownames(z) %in% inh.inh.dn )] <- "inh.dn"
        z$Cat[which(rownames(z) %in% inh.inh.up )] <- "inh.up"
        z$Cat[which(rownames(z) %in% inh.rand.dn )] <- "rand.dn"
        z$Cat[which(rownames(z) %in% inh.rand.up )] <- "rand.up"
      
        return(z)
      })
      return(x)
    }
  
    inh.de <- lapply(inh, lm.CIBx.2)
    
    success.CIBx.2 <- function(x) {
      x <- lapply(x, function(y) {
        # check for consistency in expected and observed fold-change
        y$fc.sign <- sign(y$log2fc)
        y$fc.dir <- sapply(strsplit(y$Cat, "\\."), "[", 2)
        y$fc.cons <- FALSE
        y$fc.cons[is.na(y$fc.dir)] <- TRUE # NAs occur for genes with simulated DE. these can be TRUE, because there is no expected fc direction
        y$fc.cons[y$fc.dir == "up" & y$fc.sign == 1] <- TRUE
        y$fc.cons[y$fc.dir == "dn" & y$fc.sign == -1] <- TRUE
        
        # convert to factor
        y$Cat <- factor(y$Cat, levels = c("inh.dn", "inh.up", "Neg/Marker", "rand.dn", "rand.up"))
        
        # extract out hits with padj < 0.05
        sig <- y[which(y$padj < 0.05 & y$fc.cons),] # significant p-value
       
        # hit rates!
        total.genes <- table(y$Cat)
        hit.genes <- table(sig$Cat)
      
        hit.rate <- hit.genes / total.genes
        
       # return
        return(hit.rate)
      })
        return(x)
    }
    
    inh.hit.rates <- lapply(inh.de, success.CIBx.2)
    inh.hit.rates <- lapply(inh.hit.rates, function(x) do.call("rbind", x))
  

## Reorganise dataframe
  # create
  plot.data <- melt(inh.hit.rates)
    
  # convert NaN to zero
  plot.data$value[which(is.na(plot.data$value))] <- 0
    
  ## Move information to new columns
    # on proportion confounds
    plot.data$Proportion <- "No Composition Change" 
    plot.data$Proportion[grep("up10", plot.data$L1)] <- "+10% Exc"
    plot.data$Proportion <- factor(plot.data$Proportion, levels = c("No Composition Change", "+10% Exc"))
  
  ## Rename
    plot.data$Var1 <- substr(plot.data$Var1, 1, 3)
    # plot.data$Var1 <- paste0("CIBx-isolated ", plot.data$Var1)
    
    ## Separate dataframe for negative controls
    plot.data2 <- plot.data[which(plot.data$Var2 == "Neg/Marker"),]
    ct.colours <- c("#009392", "#eeb479", "#e88471", "#e9e29c", "#fdfbe4")
    names(ct.colours) <- substr(ct, 1, 3)
    
    ## Final processing of position controls' plotting data
    plot.data <- plot.data[which(plot.data$Var2 != "Neg/Marker"),]
    plot.data$Var2 <- factor(plot.data$Var2)
    levels(plot.data$Var2) <- c("Inh Marker | Downregulated", "Inh Marker | Upregulated", "Non-marker | Downregulated", "Non-marker | Upregulated")
    marker.cols <- brewer_pal(palette = "Paired")(8)[c(3,4,7,8)]
    
    
## Plot non-specific DE
    pdf(file = "Plots - Inh, TPR.pdf", height = 3.5, width = 5)
    ggplot(plot.data, aes(x = Var1, colour = Var2, fill = Var2, y = value)) +
      stat_summary(fun = function(x) {x}, geom = "point", position = position_dodge(width = 0.7), size = 2) +
      geom_col(position = position_dodge(width = 0.7), width = 0.1, show.legend = FALSE) +
      theme_bw() +
      facet_wrap(~Proportion, ncol = 2) +
      scale_colour_manual(values = marker.cols) +
      scale_fill_manual(values = marker.cols) +
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(nrow = 2), colour = guide_legend(nrow = 2)) +
      labs(y = "% of Perturbed Genes\nReaching Significance", x = "CIBx-isolated Component of Mixture Expression") 
    dev.off()
    
    pdf(file = "Plots - Inh, FPR.pdf", height = 3.5, width = 3)
    ggplot(plot.data2, aes(x = Var1, y = value)) +
      geom_col(position = position_dodge(width = 0.7), colour = "black", fill = "black", width = 0.1, show.legend = FALSE) +
      stat_summary(fun = function(x) {x}, geom = "point", colour = "black", position = position_dodge(width = 0.7),  size = 3) +
      theme_bw() +
      facet_wrap(~Proportion, ncol = 1) +
      # scale_y_continuous(expand = c(0,0), limits = c(-0.01, )) +
      # theme(panel.border = element_blank(), axis.line.y = element_line()) +
      guides(colour = guide_legend(title = "Ct-specific expression"), fill = guide_legend(title = "Ct-specific expression")) +
      labs(y = "% of Non-Perturbed Genes\nReaching Significance", x = "CIBx-isolated Component of\nMixture Expression") 
    dev.off()
    
 
