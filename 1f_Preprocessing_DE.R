## In this (revised) script, I will simulate mixtures for testing differential expression

################################################################################################################################ #
## Generic setup  ----

## Clear environment
rm(list = ls())
gc()

## Global options
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/DE_simulations/")
options(stringsAsFactors = FALSE)

## Functions
source("../../../Scripts/Fun_Composition.R")
source("../../../Scripts/Fun_Preprocessing.R")


################################################################################################################################ #
## Load and organise the source data  ----

# using HCA data from the MTG

## Load
  load("../../../Data/Preprocessed/snme_SeuratObjects.rda")
  HCA.cells <- as.data.frame(obj$HCA@assays$RNA@counts)
  HCA.meta <- obj$HCA@meta.data
  
## Merge celltypes
  y <- x <- HCA.meta$orig.celltype
  y[grep("^Ex", x)] <- "Excitatory"
  y[grep("^In", x)] <- "Inhibitory"
  y[grep("^Astro", x)] <- "Astrocytes"
  y[grep("^Oli", x)] <- "Oligodendrocytes"
  y[grep("^OPC", x)] <- "OPC"
  
  HCA.meta$merged <- y
  
  # set celltype label to the colnames of the expression dataframe
  colnames(HCA.cells) <- y
  
  # cpm normalise
  # HCA.cells <- apply(HCA.cells, 2, function(x) {
  #   lib.size <- 10^6 / sum(x)
  #   x <- x * lib.size  
  #   return(x)
  # })
  
  # HCA.cells <- as.data.frame(HCA.cells)
  
  # filter to highly expressed genes
  min.n <- min(table(colnames(HCA.cells)))
  keep <- which(rowSums(HCA.cells > 10) > min.n) # greater than 1 count in at least 230 samples (i.e. the number of samples in the smallest group)
  
  HCA.cells <- HCA.cells[keep,]
  
## Create new HCA signature for marker generation
  
  HCA.sig <- list()
  
  for(k in levels(as.factor(y))) {
    use <- which(colnames(HCA.cells) == k)
    
    temp <- apply(HCA.cells[,use], 2, function(x) {
      lib.size <- 10^6 / sum(x)
      x <- x * lib.size
      return(x)
    })
    
    HCA.sig[[k]] <- rowMeans(temp)
  }
  HCA.sig <- as.data.frame(do.call("cbind", HCA.sig))
  
## Clean workspace
  rm(obj)
  gc() # as it's a large file...
  
  
################################################################################################################################ #
## The simulation  ----       

## Changes
  # nSims reduced from 20 to 10. this is the number of simulations within each, err, band
  # dat is now HCA
  # ct reflects the celltypes in HCA
  # props spans a bigger range


## General simulation parameters 
  ct <- c("Astrocytes", "Excitatory", "Inhibitory", "Oligodendrocytes", "OPC") # which celltypes appear in the mixture
  change.ct <- "Excitatory" # one of the cts to perturb in both abundance and expression
  nPerGroup <- 50  # number of mixtures to make for each group
  
## Parameters for perturbing celltype proportions
  ## On the range of proportions to simulate
    # each level defines a different simulation. [1] and [2] define the min and max proportion for group "alpha", while [3] and [4] do so for group "beta"
    props <- list() 
    base.prop <- data.frame(c(200, 300, 200, 300)) # this sets the "baseline" proportion range the ct specified in change.ct
    colnames(base.prop) <- change.ct
    abundance.change <- -39:39 * 5 # the amount of change in min/max for beta relative to base.prop
  
    for (j in abundance.change) { 
      x <- base.prop
      x[3:4,] <- x[3:4,] + j
      props[[paste0("up", j, ".rep1")]] <- x
      
      # if (abs(j) < 100) { # for perturbations of up to 20%, repeat the simulation thrice total!
        props[[paste0("up", j, ".rep2")]] <- x
        # props[[paste0("up", j, ".rep3")]] <- x
      # }
    }
    
    names(props) <- gsub("up-", "dn", names(props)) # prior to this, the name upX indicate an increase of X in change.ct, while up-X meant a decrease

  ## Parameters for perturbing gene expression
    geneList <- list()
    
    # get marker genes
    geneList$Markers <- nTopFeatures(signature = HCA.sig, n = 100, alg = "diff")
    
    # get genes which are markers of no cell-types
    geneList$Random <- nRandomFeatures(signature = HCA.sig[-which(rownames(HCA.sig) %in% geneList$Markers$EnsID),], n = 20, seed = 100) # seed = 100 was chosen as one possibility where no genes in this list are also "markers", using length(intersect(geneList$Markers$EnsID, geneList$Random))
    geneList$Random <- data.frame(EnsID = geneList$Random, forCelltype = "Random")
    
    # restrict markers to just excitatory markers
    geneList$Markers <- geneList$Markers[which(geneList$Markers$forCelltype %in% c("Excitatory")),]
    
    # combine
    geneList <- rbind(geneList$Markers, geneList$Random)
    
    # half will be upregulated, and half downregulated
    which.up <- c(sample(x = 1:100, replace = FALSE, size = 50),
                  sample(x = 101:200, replace = FALSE, size = 50))
    genesUp <- geneList[which.up,]
    genesDown <- geneList[-which.up,]

## Simulate
  sims <- lapply(props, function (x) { 
    # create simulations with altered composition
    x <- confound.proportion(dat = HCA.cells, ct = ct, setProps = x, nPerGroup = nPerGroup, keepPureExpression = TRUE)
    
    # layers altered expression onto the altered composition 
    x <- confound.expression(x, genesUp = genesUp, replace.starting.level = FALSE, genesDown = genesDown, absoluteFC = 1.1, noise = 0, id = "confound.pe.1.1")
    x <- confound.expression(x, genesUp = genesUp, replace.starting.level = FALSE, genesDown = genesDown, absoluteFC = 1.3, noise = 0, id = "confound.pe.1.3")
    x <- confound.expression(x, genesUp = genesUp, replace.starting.level = FALSE, genesDown = genesDown, absoluteFC = 1.5, noise = 0, id = "confound.pe.1.5")
    x <- confound.expression(x, genesUp = genesUp, replace.starting.level = FALSE, genesDown = genesDown, absoluteFC = 2.0, noise = 0, id = "confound.pe.2.0")

    # per above, but perturb expression only in the Excitatory component of expression
    x <- confound.expression(x, genesUp = genesUp, replace.starting.level = FALSE, genesDown = genesDown, absoluteFC = 1.1, noise = 0, id = "confound.pe.exc.1.1", ct.specific = "Excitatory")
    x <- confound.expression(x, genesUp = genesUp, replace.starting.level = FALSE, genesDown = genesDown, absoluteFC = 1.3, noise = 0, id = "confound.pe.exc.1.3", ct.specific = "Excitatory")
    x <- confound.expression(x, genesUp = genesUp, replace.starting.level = FALSE, genesDown = genesDown, absoluteFC = 1.5, noise = 0, id = "confound.pe.exc.1.5", ct.specific = "Excitatory")
    x <- confound.expression(x, genesUp = genesUp, replace.starting.level = FALSE, genesDown = genesDown, absoluteFC = 2.0, noise = 0, id = "confound.pe.exc.2.0", ct.specific = "Excitatory")
    
    return(x)
  })
  
  # save
  save(sims, file = "../Results/Revisions/DE_simulations/Sims (Final).rda")
  
  # save a single simulation level, for when you want to explore its structure quickly
  example.simulation.no.confound <- sims$up0.rep1
  example.simulation.up.50 <- sims$up50.rep1
  
  save(example.simulation.no.confound, example.simulation.up.50, file = "Sims (Examples) (Final).rda")
