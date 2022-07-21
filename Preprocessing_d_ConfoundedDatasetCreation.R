## In this script, I will simulate mixtures for testing differential expression

## I have set 

################################################################################################################################ #
## Generic setup  ----

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


################################################################################################################################ #
## Load and organise the source data  ----

# using CA data from the MTG

## Load
  load("Preprocessed/SeuratObjects.rda")
  CA.cells <- as.data.frame(obj$CA@assays$RNA@counts)
  CA.meta <- obj$CA@meta.data
  
## Remove microglia
  remove <- grep("Micro", CA.meta$orig.celltype)
  CA.cells <- CA.cells[,-remove]
  CA.meta <- CA.meta[-remove,]
  
## Merge celltypes
  y <- x <- CA.meta$orig.celltype
  y[grep("^Ex", x)] <- "Excitatory"
  y[grep("^In", x)] <- "Inhibitory"
  y[grep("^Astro", x)] <- "Astrocytes"
  y[grep("^Oli", x)] <- "Oligodendrocytes"
  y[grep("^OPC", x)] <- "OPC"
  
  CA.meta$merged <- y
  
  # set celltype label to the colnames of the expression dataframe
  colnames(CA.cells) <- y
  
 
  # filter to highly expressed genes
  min.n <- min(table(colnames(CA.cells)))
  keep <- which(rowSums(CA.cells > 10) > min.n) # greater than 1 count in at least 238 samples (i.e. the number of samples in the smallest group)
  
  CA.cells <- CA.cells[keep,] # no need for filtering CA.meta, no longer being used
  
## Create new CA signature for marker generation
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
  write.csv(CA.sig, "Preprocessed/ConfoundingComposition_Signature.csv")
  
## Clean workspace
  rm(obj)
  gc() # as it's a large file...
  
  
################################################################################################################################ #
## The simulation  ----       

## Changes
  # nSims reduced from 20 to 10. this is the number of simulations within each, err, band
  # dat is now CA
  # ct reflects the celltypes in CA
  # props spans a bigger range


## General simulation parameters 
  ct <- c("Astrocytes", "Excitatory", "Inhibitory", "Oligodendrocytes", "OPC") # which celltypes appear in the mixture
  change.ct <- "Excitatory" # one of the cts to perturb in both abundance and expression
  nPerGroup <- 50  # number of mixtures to make for each group
  
## Parameters for perturbing celltype proportions
  ## On the range of proportions to simulate
    # each level defines a different simulation. [1] and [2] define the min and max proportion for group "alpha", while [3] and [4] do so for group "beta"
    props <- list() 
    base.prop <- data.frame(c(200, 300, 200, 300)) # this sets the "baseline" proportion range the ct specified in change.ct. This is as a fraction of 500 cells, i.e. 40-60%.
    colnames(base.prop) <- change.ct
    # abundance.change <- -39:39 * 5 # the amount of change in min/max for beta relative to base.prop
    abundance.change <- -5:5 * 5 # a smaller number of samples, for replication purposes
  
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

## Which genes are you going to perturb?
  geneList <- list()
  
  # get marker genes
  geneList$Markers <- nTopFeatures(signature = CA.sig, n = 100, alg = "diff")
  
  # get genes which are markers of no cell-types
  geneList$Random <- nRandomFeatures(signature = CA.sig[-which(rownames(CA.sig) %in% geneList$Markers$EnsID),], n = 20, seed = 100) # seed = 100 was chosen as one possibility where no genes in this list are also "markers", using length(intersect(geneList$Markers$EnsID, geneList$Random))
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
  # each simulated dataset (i.e., level of "props") requires about a minute to run, and ~200MB of RAM
  sims <- lapply(props, function (x) { 
    # create simulations with altered composition
    x <- confound.proportion(dat = CA.cells, ct = ct, setProps = x, nPerGroup = nPerGroup, keepPureExpression = TRUE)
    
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
  save(sims, file = "Preprocessed/ConfoundingComposition_Simulations.rda")
  