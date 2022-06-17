##################################################################################################################################
## Packages

## Composition estimation
require(DeconRNASeq)
require(BrainInABlender)
require(xCell)
require(dtangle)
require(linseed)
require(MuSiC)
require(xbioc) # for music
require(WGCNA) # for co-expression-based composition
require(limma) # for co-expression-based composition
require(multtest) # for co-expression-based composition
source("/Volumes/Data1/PROJECTS/BrainCellularComposition/Scripts/Fun_CIBERSORT.R")

## Plotting
require(ggplot2)
require(cowplot)
require(reshape2)
require(viridis)
require(pals)
require(variancePartition)


##################################################################################################################################
## Run composition estimation algorithms

## Generic function to text output...
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


## Run DeconRNASeq. Returns dataframe
run.DRS <- function(mixture, signature) {
  res <- as.data.frame(DeconRNASeq(mixture, signature, use.scale = TRUE)$out.all)
  rownames(res) <- colnames(mixture)
  return(res)
}


## Run xCell. 1 note
run.xCell <- function(mixture, raw = FALSE) {
  # calculate enrichment
  res <- list()
  if (raw) res$Raw <- rawEnrichmentAnalysis(expr = mixture, signatures = xCell.data$signatures, genes = xCell.data$genes, parallel.sz = 2)
  res$Transformed <- xCellAnalysis(expr = mixture, save.raw = TRUE, rnaseq = TRUE, parallel.sz = 2)
  
  # reformat, and filter to neurons and astrocytes
  res <- lapply(res, function(x) { 
    as.data.frame(t(x[c("Neurons", "Astrocytes"),])) 
  })
  
  # clean extraneous variables from the environment
  env <- ls(pos = ".GlobalEnv")
  rm(list = env[c(grep("progressBar", env),
                  grep("nSamples", env),
                  grep("iSample", env))], pos = ".GlobalEnv")
  
  # return
  if (raw) {
    return(res) # Note1: output is a list of two dataframes
  } else {
    return(res$Transformed)
  }
}


## Run BrainInABlender. 1 note.
run.Blender <- function(mixture, OPCs = FALSE) {
  # setup
  n <- ncol(mixture)
  mixture$geneColumn <- rownames(mixture)
  
  # calculate enrichment
  blend <- quiet(Sir_UnMixALot(userInput = mixture, dataColumns = 1:n, geneColumn = n+1, species = "Human"))
  
  # reformat
  res <- list()
  
  
  ct <- c("Astrocyte", "Neuron_All", "Oligodendrocyte", "Endothelial", "Microglia")
  if(OPCs) ct <- c(ct, "Oligodendrocyte_Immature")
  
  res <- t(blend$AveragePrimary_CellTypeIndex[ct,])
  res <- as.data.frame(res)
  
  # rename columns
  if(OPCs) {
    colnames(res) <- c("Astrocytes", "Neurons", "Oligodendrocytes", "Endothelia", "Microglia", "OPC")
  } else {
    colnames(res) <- c("Astrocytes", "Neurons", "Oligodendrocytes", "Endothelia", "Microglia")
  }
    
  # return
  return(res)
  
}


## Run CIBERSORT. 1 note.
write.CIB <- function(data, dir) {
  data <- cbind(rownames(data), data)
  colnames(data)[1] <- "EnsID"
  write.table(data, file = dir, sep = "\t", quote = FALSE, row.names = FALSE)
}

run.CIB <- function(directoryCIB = "/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/CIB/", mixString, sigString = FALSE, sigObject = FALSE, from.file = TRUE) {
  # deconvolve
  if (from.file) {
    res <- CIBERSORT(sig_matrix = paste0(directoryCIB, sigString), # Note1: CIBERSORT reads in its files from a directory, not R's environment. Specify this using directoryCIB
                     mixture_file = paste0(directoryCIB, mixString))
  } else {
    # write out signature
    seed.sig <- sample(1:1000, 1)
    name.sig <- paste0(directoryCIB, "TemporarySig_",  seed.sig, ".txt")
    write.CIB(data = sigObject, dir = name.sig)
    
    # path of mixture
    # note: not writing it to disk, due to the potential filesize (several hundred MB)
    name.mix <- paste0(directoryCIB, mixString)
    
    # run CIB
    res <- CIBERSORT(sig_matrix = name.sig, mixture_file = name.mix)
    
    # clean
    file.remove(name.sig)
  }
  
  # reformat
  res <- as.data.frame(res)
  res <- res[,1:(ncol(res) - 3)] # Note1: This removes some of CIBERSORT's output that is extraneous (for our purposes)
  return(res)
}

## Run dtangle
run.DTA <- function(mixture, signature, alg, q) { # alg %in% c("diff", "ratio", "p.value", "regression")
  # bind dataframes
  common <- rownames(mixture)[which(rownames(mixture) %in% rownames(signature))]
  dat <- cbind(log2(mixture[common,] + 0.5), log2(signature[common,] + 0.5))
  dat <- as.data.frame(t(dat))
  
  # define where signature samples reside in dat
  ps <- list()
  for (j in colnames(signature)) { ps[[j]] <- grep(paste0("^", j, "$"), rownames(dat)) }

  # find markers
  markers <- find_markers(dat, marker_method = alg, data_type = "rna-seq", pure_samples = ps)
  
  # deconvolve
  quant <- lapply(markers$V, function(x) { quantile(x, 1-q) } )
  for (j in 1:length(quant)) {
    if (quant[j] == Inf) { quant[j] <- max(grep("Inf", markers$V[[j]]) + 1) } # a hack for when there are more "infs" than q markers
  }
  n <- sapply(1:length(markers$V), function(i) { max(which(markers$V[[i]] > quant[[i]])) } )
  res <- as.data.frame(dtangle(Y = dat, 
                               pure_samples = ps, 
                               markers = markers$L, 
                               n_markers = n, 
                               marker_method = alg, 
                               data_type = "rna-seq")$estimates)
  
  # return(res)
  return(res[1:(nrow(res) - ncol(signature)),])
}

## Run linseed
run.linseed <- function(mixture, nCelltypes, write.plots = TRUE, write.data = FALSE, iters = 100) {
  # setup object
  lin <- LinseedObject$new(mixture) # this auto-filters to the 10,000 most highly-expressed genes

  # run collinearity network
  lin$calculatePairwiseLinearity() # what I find thoroughly displeasing about this format is it prevents the use of help
  lin$calculateSpearmanCorrelation()
  lin$calculateSignificanceLevel(iters = iters) # takes ~20m
  if (write.plots) print(lin$significancePlot(threshold = 0.01))
  
  # filter to genes that have significant ML
  lin$filterDatasetByPval(0.01) # From 10,000 genes to 4364 - I think there can be a some small variation in this number (as the first run yielded 4377)
  
  # plot singular value decomposition
  if (write.plots) print(lin$svdPlot())
  
  # set number of cell-types based not on SVD (as the pipeline would suggest), but instead by our a priori knowledge
  lin$setCellTypeNumber(nCelltypes)
    
  # project (verb)
  lin$project("full")
  if (write.plots) print(lin$projectionPlot(color = "filtered")) 
  
  # deconvolve
  lin$smartSearchCorners(dataset = "filtered", error = "norm") # "Final vector is: 10 9"
  lin$deconvolveByEndpoints()
  if (write.plots) print(plotProportions(lin$proportions))
  
  # output results!
  res <- list()
  res$Raw <- as.data.frame(t(lin$proportions))
  res$Transformed <- as.data.frame(t(apply(res$Raw, 1, function(x) x / sum(x))))
  if(write.data) res$Data <- lin
  
  return(res)
}

## Run co-expression analyses for cell-type enrichment
run.coex <- function(mixture, signature, only.threshold = FALSE, sft = "auto", output.all = TRUE) {
  # setup
  dat <- list()
  log <- log2(mixture + 0.5) # normally-distribute data
  log <- log[which(apply(log, 1, function(x) sd(x) != 0)),]

  # if required, determine optimal soft-thresholding power
  if (sft == "auto") {
    # choose a set of powers to test
    powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
  
    # calculate network structure
    dat$sft <- pickSoftThreshold(t(log), powerVector = powers, verbose = 5, networkType = "signed", corFnc = "bicor")$fitIndices
    
    # auto-detect the optimal sft as the lowest with r2 > 0.8 & mediank < 100
    sft <- min(dat$sft$Power[which(dat$sft$SFT.R.sq > 0.8 & dat$sft$median.k. < 100)])
  }

  # construct network
  dat$net <- blockwiseModules(t(log), power = sft, deepSplit = 4, minModuleSize = 150, 
                             mergeCutHeight = 0.2, detectCutHeight = 0.9999,
                             corType = "bicor", networkType = "signed", pamStage = FALSE, pamRespectsDendro = TRUE,
                             verbose = 3, saveTOMs = FALSE, maxBlockSize = 30000, numericLabels = TRUE)
  
  ## Network relabelling
    # labelling the modules with a colour tag
    dat$modules <- as.data.frame(table(dat$net$colors)); colnames(dat$modules) <- c("Label", "N") 
    dat$modules$Label <- paste0("M", dat$modules$Label)    
    dat$modules$Color <- c("grey", labels2colors(dat$modules$Label[-1]))
    dat$moduleLabel <- paste0("M", dat$net$colors) # variable storing each genes' module label
    dat$moduleColor <- dat$modules$Color[match(dat$moduleLabel, dat$modules$Label)] # Similar to above but for colour
      
    # calculating kMEs
    dat$KMEs <- signedKME(t(log), dat$net$MEs, outputColumnName = "M") # signedKME calculates eigengene-based connectivity i.e. module membership. Also, outputColumnName gives us a prefix. Also, the corFnc defaults to Pearson
    dat$kme <- data.frame(rownames(log), dat$moduleColor, dat$moduleLabel, dat$KMEs)
    colnames(dat$kme)[1] <- "Symbol"
      
    # re-assign genes to modules by kME. Any not passing the criteria of correlation pvalue < 0.05, and kME > 0.5, will be moved to a junk module (M0)
    info <- c(1:3)
    dat$kmedata <- dat$kme[,-info] 
    dat$pvalBH <- dat$kmedata; dat$pvalBH[,] <- NA 
  
    for (k in c(1:ncol(dat$pvalBH))) {
      p <- mt.rawp2adjp(corPvalueStudent(dat$kmedata[,k], nSamples = ncol(log)), proc = "BH") 
      dat$pvalBH[,k] <- p$adjp[order(p$index),2]
    }
  
    dat$kme$newModule <- "NA" 
    
    for (k in c(1:nrow(dat$kmedata))) {
      if (k%%1000 == 0) print(k); 
      m <- which(dat$kmedata[k,] == max(dat$kmedata[k,]))
      if ((dat$pvalBH[k,m] < 0.05) & (dat$kmedata[k,m] > 0.5)) dat$kme$newModule[k] <- as.character(colnames(dat$kmedata)[m])
    }
  
    # assign genes not associated to any module to M0, though that's in a new column so we can compare!
    dat$kme$newModule[which(dat$kme$newModule%in%"NA")] <- "M0" 
        
    # replacing the old moduleLabel and moduleColor columns with the updated newModule and newColor columns respectively
    dat$kme$newColor <- dat$kme$dat.moduleColor[match(dat$kme$newModule, dat$kme$dat.moduleLabel)] 
    dat$kme$dat.moduleLabel <- dat$kme$newModule; dat$kme$dat.moduleColor <- dat$kme$newColor 
    dat$kme <- dat$kme[,-grep("newModule", colnames(dat$kme))]; dat$kme <- dat$kme[,-grep("newColor", colnames(dat$kme))] 

    # saving kMEs 
    dat$modules <- dat$modules$Label[-1] 
    dat$kmeTable <- dat$kme[,info]
    for(k in c(1:length(dat$modules))) {
      dat$kmeTable <- cbind(dat$kmeTable, dat$kmedata[,match(dat$modules[k],colnames(dat$kmedata))]) 
      colnames(dat$kmeTable)[ncol(dat$kmeTable)] <- paste("kME", dat$modules[k], sep = "_")    
      
      dat$kmeTable <- cbind(dat$kmeTable, dat$pvalBH[,match(dat$modules[k],colnames(dat$pvalBH))]) 
      colnames(dat$kmeTable)[ncol(dat$kmeTable)] <- paste("pvalBH", dat$modules[k], sep = "_")
    }
    
    table(dat$kmeTable$dat.moduleLabel)

    # saving Module Eigengenes 
    dat$me <- data.frame(colnames(log), dat$net$MEs) # Sample names bound to module eigengenes
    colnames(dat$me)[-1] <- paste("M", substr(colnames(dat$me)[-1], start = 3, stop = 100), sep = "") # ... the most convoluted method of renaming I've seen 
    rownames(dat$me) <- dat$me[,1]
    dat$me <- dat$me[,-1]
    
    # get list of genes in every module
    dat$membership <- list()
    for(k in dat$modules) { dat$membership[[k]] <- dat$kmeTable$Symbol[which(dat$kmeTable$dat.moduleLabel == k)] }

  ## Compare modules to cell-types in the signature
    # find markers
    marks <- nTopFeatures(signature, n = 50)
    
    # odds ratio of marker enrichment
    e <- list()
    for (k in names(table(marks$forCelltype))) {
      e[[k]] <- list()
      for(l in dat$modules) e[[k]][[l]] <- as.numeric(overEnrich(marks$EnsID[which(marks$forCelltype == k)], dat$membership[[l]], backgroundList = rownames(log))[2]) # indexing 2 takes the odds ratio
      e[[k]] <- do.call("rbind", e[[k]])
      colnames(e[[k]]) <- k
    }

    ids <- list()
    e <- do.call("cbind", e)
    for(j in colnames(e)) {
      # filter to modules which have the highest OR in this celltype
      m <- colnames(e)[apply(e, 1, which.max)]
      m <- which(m == j)
      
      # get the module with the highest OR for this celltype 
      ids[[j]] <- rownames(e)[which.max(e[,j])]
    }
    
  # collect module eigengenes as composition estimates
  res <- as.data.frame(matrix(nrow = ncol(log), ncol = length(ids)))
  rownames(res) <- colnames(log)
  colnames(res) <- names(ids)
  
  for(j in colnames(res)) {
    res[,j] <- dat$me[,ids[[j]]]
  }
    
  # output
  if (output.all) {
    dat$composition <- res
    return(dat)
  } else {
    return(res)
  }
}

run.music <- function(mixture, # mixture counts 
                      signature, # signature counts
                      signature.meta, # metadata for the signatures. needs $Indivisual, $Celltype, and $Celltype.family (latter if heirarchical)
                      family = NA, # input a list to follow this path
                      use.meta.column = "Celltype",
                      drop.ct = FALSE, # set to TRUE if you want to remove any cells labelled "Drop"
                      rpkm = FALSE, # set to true if mixture and signature are rpkm normalised
                      family.markers = "generate", # if generate, the script generates them, otherwise please provide your own!
                      lib.sizes) { # necessary if(rpkm)
  if(rpkm) {
    est <- NA
  } else { # standard approach is to run on count-level data
    print("Starting preprocessing!")
    
    # setup meta
    signature.meta$Celltype <- signature.meta[,use.meta.column]
    
    if(drop.ct) {
      drop <- which(signature.meta$Celltype == "Drop")
      
      signature <- signature[,-drop]
      signature.meta <- signature.meta[-drop,]
      
    }
    
    # setup mixture
    music.mix <- ExpressionSet(assayData = as.matrix(mixture))
    
    # setup signature
    music.sig <- ExpressionSet(assayData = as.matrix(signature),
                               phenoData = AnnotatedDataFrame(signature.meta))
    
    
    if (class(family) != "list") { # that is, if no clusters of cell-types are defined
      est <- music_prop(bulk.eset = music.mix, sc.eset = music.sig, clusters = "Celltype", samples = "Individual") # here, we set the clusters to be each celltype
      est <- as.data.frame(est$Est.prop.weighted)
      return(est)
      
    } else {
      ## First, define markers for the family
        # create a cpm matrix
        # dat <- apply(signature, 2, function(x) {
        #   lib.size <- 10^6 / sum(x)
        #   x <- x * lib.size
        #   return(x)
        # })
        
        if (family.markers == "generate") {
          # for loop instead of apply
          dat <- signature
          for(j in colnames(dat)) {
            lib.size <- 10^6 / sum(dat[,j])
            dat[,j] <- dat[,j] * lib.size
          }
          
          # run dtangle function "find_markers"
          dat <- as.data.frame(t(log2(dat + 0.5)))
          pure_samples <- list(); for(j in names(family)) pure_samples[[j]] <- grep(j, signature.meta$Celltype.family)
          family.markers <- find_markers(dat, marker_method = "diff", data_type = "rna-seq", pure_samples = pure_samples)
          # family.markers <- sapply(family.markers$L, names)
          
          # and then remove those genes, keeping all others!
          family.markers <- lapply(family.markers$L, function(x) {
            x <- colnames(dat)[-which(colnames(dat) %in% names(x))]
            x <- sample(x, 5000, replace = FALSE)
          })
          
          rm(dat)
        }      
        
      ## Finally deconvolve...
        print("Starting deconvolution :)")
        est <- music_prop.cluster(bulk.eset = music.mix, sc.eset = music.sig, clusters = "Celltype", samples = "Individual",
                                  group.markers = family.markers, groups = "Celltype.family", clusters.type = family)
        est <- as.data.frame(est$Est.prop.weighted.cluster)
    }
  }
  
  return(est)
}


## Function for assigning names to automatically-determined celltypes from Linseed or Coex
define.celltypes <- function(x, algorithm.type, celltype.signature, return.confidence = TRUE, nMarks.x = 10000, nMarks.sig = 100, bg) {
  # algorithm type must be one of linseed or coex
  stopifnot(algorithm.type %in% c("Linseed", "Coex"))
  
  # parameters
  thresh.or <- 5
  thresh.p <- 10^-5
  ct <- colnames(celltype.signature)
  
  # collect markers from within linseed data
  if(algorithm.type == "Linseed") {
    # get markers within linseed's celltypes
    x.sig <- x$signatures
    x.marks <- nTopFeatures(x.sig, n = nMarks.x) 
    x.marks <- x.marks[which(!(is.na(x.marks$EnsID))),]
    x.marks.list <- list()
    for(j in names(table(x.marks$forCelltype))) x.marks.list[[j]] <- x.marks$EnsID[which(x.marks$forCelltype == j)] 
  }
  
  # collect markers from within coex data
  if(algorithm.type == "Coex") {
    x.marks.list <- x$membership
  }
  
  # collect markers within signatures
  if (nMarks.sig == FALSE) { # that is, if the supplied signature is simply a list of Ensembl IDs
    sig.marks.list <- celltype.signature
  } else {
    sig.marks <- nTopFeatures(celltype.signature, n = nMarks.sig) 
    sig.marks.list <- list()
    for(j in names(table(sig.marks$forCelltype))) sig.marks.list[[j]] <- sig.marks$EnsID[which(sig.marks$forCelltype == j)] 
  }
  
  # overlap
  conf <- lapply(sig.marks.list, function(y) {
    y <- sapply(x.marks.list, function(z) {
      overEnrich(list1 = y, list2 = z, backgroundList = bg)
    }) 
    y <- as.data.frame(t(y))
    y$OddsRatio <- as.numeric(as.character(y$OddsRatio))
    y$p.value <- as.numeric(as.character(y$p.value))
    return(y)
  })
  
  ## Assign
    # collect data
    ranks <- do.call("rbind", conf)
    ids <- strsplit(rownames(ranks), "\\.")
    ranks$sig <- sapply(ids, "[", 1)
    ranks$x <- sapply(ids, "[", 2)
    
    # remove any associations for which p < thresh.p
    ranks <- ranks[which(ranks$p.value < thresh.p),]
    
    # rank by odds ratio (so as to not be sensitive to imbalances in n) (p.value is accounted for by thresh.p)
    ranks <- ranks[order(ranks$OddsRatio, decreasing = TRUE),] 
    
    # iteratively assign celltypes
    for(j in 1:length(ct)) {
      # get the overlap with the highest odds ratio
      top.hit <- ranks[1,]
      
      # assign a derived cell-type to a supervised celltype
      names(ct)[grep(top.hit$sig, ct)] <- top.hit$x
      
      # remove any results for the derived or supervised celltype
      remove <- unique(c(which(ranks$sig == top.hit$sig),
                         which(ranks$x == top.hit$x)))
      ranks <- ranks[-remove,]
      
      # break the loop if ranks is empty (i.e., when there is a celltype which cannot be assigned)
      if(nrow(ranks) == 0) break
    }
    
    # label celltypes which failed to meet any assignment criteria
    names(ct)[which(is.na(names(ct)))] <- "Unassigned"  

  ## Collect composition estimates
    if (algorithm.type == "Linseed") {
      comp <- (apply(x$proportions, 2, function(x) x / sum(x)))
      comp <- as.data.frame(t(comp))
    }
    
    if (algorithm.type == "Coex") {
      comp <- x$me[,names(ct)[which(names(ct) != "Unassigned")]]
    }
  
  ## Output
  output <- list()
  output$comp <- comp
  output$assignments <- ct
  if(return.confidence) output$conf <- conf
  return(output)
}

##################################################################################################################################
## Calculate statistics

## Calculate the markers of a column relative to other columns in its dataframe
nTopFeatures <- function(signature, n, alg = "diff") {
  # the signature is a dataframe of aggregated cell-types x genes
  # n is the number of markers per cell-type, giving a total of n*ncol(signature) markers
  
  # setup
  sig <- t(log2(signature + 1))
  
  ps <- list() # ps is required by dtangle's find_markers() function.  
  for (j in rownames(sig)) { ps[[j]] <- match(j, rownames(sig)) }
  
  res <- as.data.frame(matrix(ncol = 2, nrow = n*nrow(sig)))
  colnames(res) <- c("EnsID", "forCelltype")
  
  # find markers, using find_markers() from the dtangle package
  allMarkers <- find_markers(sig, marker_method = alg, data_type = "rna-seq", pure_samples = ps)
  nMarkers <- sapply(allMarkers$V, function(x) { names(x)[1:n] } ) 
  
  # add to res
  counter <- 1
  for (j in colnames(nMarkers)) {
    for(k in 1:n) {
      res$EnsID[counter] <- nMarkers[k,j]
      res$forCelltype[counter] <- j
      counter <- counter + 1
    }
  }
  
  # alternate situation when n = 1
  if (n == 1) {
    res$EnsID <- nMarkers
    res$forCelltype <- names(nMarkers)
  }
  
  # output
  return(res)
}

## Function for full report on list overlap
overEnrich <- function(list1, list2, backgroundList) {
  # setup dataframes
  overlap <- length(which(list1 %in% list2))
  x <- length(list2) - overlap
  y <- length(list1) - overlap
  z <- length(backgroundList) - x - y + overlap
  contingencyTable <- rbind(c(overlap, x), c(y, z))
  
  # run test
  test <- fisher.test(contingencyTable, alternative = "greater")
  
  # output
  results <- c(paste0(overlap, " / ", length(list1), "&", length(list2)), 
               round(test$estimate, 3), 
               round(test$conf.int[1], 3), 
               round(test$conf.int[2], 3), 
               test$p.value)
  names(results) <- c("Overlap", "OddsRatio", "conf95_lower", "conf95_upper", "p.value")
  return(results)
}

## Euclidean norm of a vector
eNorm <- function(x) {
  sqrt(sum(x^2))
}

## Calculate goodness of fit
write.gof <- function(measuredExp, estimatedComp, signatureUsed, returnPred = FALSE, log = FALSE) {
  # set to common order
    commonGenes <- rownames(measuredExp)[which(rownames(measuredExp) %in% rownames(signatureUsed))]
    measuredExp <- measuredExp[commonGenes,]; signatureUsed <- signatureUsed[commonGenes,]
    
  # predict expression (predExp) from the estimatedComp * signatureUsed
    predExp <- as.data.frame(matrix(nrow = length(commonGenes), ncol = ncol(measuredExp)))
    rownames(predExp) <- commonGenes
    
    for(j in 1:ncol(predExp)) {
      # storage
      a <- list()
      
      # the contribution of each cell-type to predicted expression
      for(k in colnames(signatureUsed)) { a[[k]] <- estimatedComp[j,k] * signatureUsed[,k] }
      
      # sum expression from all cell-types to a single predicted value
      predExp[,j] <- rowSums(do.call("cbind", a))
    }
  
  # statistics
      # (optional) log transform
      if(log) {
        measuredExp <- log2(measuredExp + 0.5)
        predExp <- log2(predExp + 0.5)
      }
    
    
    # setup
    stats <- as.data.frame(matrix(ncol = 6, nrow = ncol(measuredExp)))
    colnames(stats) <- c("rho", "r", "mae", "rmse", "recon", "cosine")
    rownames(stats) <- colnames(measuredExp)
    
    # spearman correlation between predicted and estimated expression
    stats$rho <- diag(cor(measuredExp, predExp, method = "s")) 
    
    # pearson correlation (comments per spearman)
    stats$r <- diag(cor(measuredExp, predExp, method = "p")) 
    
    # mean absolute error, root mean squared error, reconstruction accuracy, and cosine similarity
    for(j in 1:ncol(measuredExp)) { 
      a <- measuredExp[,j]
      b <- predExp[,j]
      
      stats$mae[j] <- mae(a, b) 
      stats$rmse[j] <- rmse(a, b) 
      stats$recon[j] <- 1 - ((eNorm(a - b) ^ 2) / (eNorm(a) ^ 2)) # note: this (I believe) is what PsychEncode calls the "reconstruction accuracy"
      stats$cosine[j] <- (sum(a * b) / (eNorm(a) * eNorm(b))) # note: this is the cosine similarity of two vectors per Wikipedia: dot product / product of euclidean norms
    }
  
  # return
    if(returnPred) {
      res <- list()
      res$predExp <- predExp
      res$stats <- stats  
    } else {
      res <- stats
    }
  
    return(res)
}

write.gof.v2 <- function(measuredExp, estimatedComp, signatureUsed, returnPred = FALSE) {
  # set to common row order
    commonGenes <- rownames(measuredExp)[which(rownames(measuredExp) %in% rownames(signatureUsed))]
    measuredExp <- measuredExp[commonGenes,]; signatureUsed <- signatureUsed[commonGenes,]
    
  # quantile normalise
    
    qn <- data.frame(signatureUsed, measuredExp)
    qn <- as.data.frame(normalize.quantiles(as.matrix(qn), copy = FALSE))
    
    signatureUsed <- qn[,1:ncol(signatureUsed)]
    measuredExp <- qn[,-c(1:ncol(signatureUsed))]

    
  # predict expression (predExp) from the estimatedComp * signatureUsed
    predExp <- as.data.frame(matrix(nrow = length(commonGenes), ncol = ncol(measuredExp)))
    rownames(predExp) <- commonGenes
    
    for(j in 1:ncol(predExp)) {
      # storage
      a <- list()
      
      # the contribution of each cell-type to predicted expression
      for(k in colnames(signatureUsed)) { a[[k]] <- estimatedComp[j,k] * signatureUsed[,k] }
      
      # sum expression from all cell-types to a single predicted value
      predExp[,j] <- rowSums(do.call("cbind", a))
    }

  ## Calculate statistics
    stats <- as.data.frame(matrix(ncol = 5, nrow = ncol(measuredExp)))
    colnames(stats) <- c("rho", "r", "mae", "rmse", "nmae")
    rownames(stats) <- colnames(measuredExp)
    
    for(j in 1:ncol(measuredExp)) { 
      a <- measuredExp[,j]
      b <- predExp[,j]
      
      stats$r[j] <- cor(log2(a+0.5), log2(b+0.5), method = "p")
      stats$rho[j] <- cor(a, b, method = "s")
      
      stats$mae[j] <- mae(a, b) 
      stats$rmse[j] <- rmse(a, b) 
      
      stats$nmae[j] <- nmae(a, b) 
      
    }

  # return
    if(returnPred) {
      res <- list()
      res$predExp <- predExp
      res$stats <- stats  
    } else {
      res <- stats
    }
  
    return(res)
}


## Mean absolute error (MAE) 
mae <- function(true, est) {
  mean(abs(true - est))
}


## Normalised mean absolute error (nmae)
nmae <- function(true, est) {
  mean(abs(true - est)) / mean(true)
}


## Mean absolute error percentage error (nmae)
mape <- function(true, est) {
  if(length(which(true == 0)) != 0) {
  zeroes <- which(true == 0)
  true <- true[-zeroes]
  est <- est[-zeroes]
  }
  mean(abs((true - est) / true))
}


## Root mean squared error (rmse)
rmse <- function(true, est) {
  sqrt(mean((true - est)^2))
}


## Normalised root mean squared error (nrmse)
nrmse <- function(true, est) {
  rmse(true, est) / mean(true)
}

## General stats function
write.stats <- function(t, e, alg, error) {
  # t for true
  # e for estimate
  # alg for algorithm
  # error is a logical indicating whether an error statistic should be calculated
  
  # setup
  e <- e[,which(colnames(e) %in% colnames(t))]
  stats <- matrix(NA, 6, length(e))
  rownames(stats) <- c("rho", "r", "rmse", "nrmse", "mae", "nmae")
  colnames(stats) <- colnames(e)
  # if(!(error)) stats <- stats[1:2,] # it is handy to see these rows as NA

  # loop
  for (ct in colnames(e)) {
      stats["rho",ct] <- cor(t[,ct], e[,ct], method = "s")
      stats["r",ct] <- cor(t[,ct], e[,ct], method = "p")  
    
    
    if (error) {
      stats["rmse",ct] <- rmse(t[,ct], e[,ct])
      stats["nrmse",ct] <- nrmse(t[,ct], e[,ct])
      stats["mae",ct] <- mae(t[,ct], e[,ct])
      stats["nmae",ct] <- nmae(t[,ct], e[,ct])
    }
  }
  
  # final formatting    
  rownames(stats) <- paste0(rownames(stats), "_", alg)
  stats <- signif(stats, 3)
  
  # output
  return(stats)
}


##################################################################################################################################
## Plotting


## Shorten the name of element_blank()
invis <- element_blank()

## Generic-but-flexible-and-oh-so-complicated scatterplot
plot.scatter <- function(t, e, ct, ylab = "Estimated Proportion", 
                         calcCor = FALSE, calcError = FALSE, colour, abline.colour = colour,
                         limits = c(0, NA), abline = TRUE, annot.pos = c(-Inf, Inf, -0.1, 1.2)) {
  # t for true (dataframe of all cell-types)
  # e for estimated (dataframe of all cell-types)
  # ct to specify the cell-type to be plotted
  # ylab is the label on the y-axis
  # cor is the correlation method (one of "s", "p", or NA)
  # error is the error method (one of "mae", "rme", or NA)
  # colour is the point colour
  # limits defines the x and y axes' limts
  
  # format
  e <- e[,ct]
  t <- t[,ct]
  
  # draw error annotation
    # collect correlation
    if (calcCor == "rho") cor <- signif(cor(t, e, method = "s"), 2)
    if (calcCor == "r") cor <- signif(cor(t, e, method = "p"), 2)
    
    # collect error
    if (calcError == "mae") err <- signif(mae(t, e), 2)
    if (calcError == "rmse") err <- signif(rmse(t, e), 2)
    if (calcError == "nmae") err <- signif(nmae(t, e), 2)
    if (calcError == "nrmse") err <- signif(nrmse(t, e), 2)
  
    # print a label
    logic <- as.numeric(c(calcCor != FALSE, calcError != FALSE))  # a logical vector for whether error and cor were calculated
    if (logic[1] > logic[2]) label <- paste0(calcCor, "=", cor)
    if (logic[1] < logic[2]) label <- paste0(calcError, "=", err)
    if (sum(logic) == 2) label <- paste0(calcError, "=", err, "\n", calcCor, "=", cor)

  # plot
  sp <- qplot(t, e, col = I(colour)) +
    theme_bw() +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    geom_smooth(method = "lm", se = FALSE, colour = abline.colour) +
    labs(x = "True", y = ylab) +
    scale_x_continuous(limits = limits) +
    scale_y_continuous(limits = limits)
  
  # (optional) annotate plot statistics
    if(sum(logic) >= 1) sp <- sp + annotate("text", x = annot.pos[1], y = annot.pos[2], hjust = annot.pos[3], vjust = annot.pos[4], label = label, colour = colour, size = 3)

  # (optional) trendline
    if(abline) sp <- sp + geom_abline(slope = 1, intercept = 0, colour = "red", linetype = 2) 

  # output
  sp
}

## Colour schemes
ct.colours <- c("black", "firebrick","dodgerblue3", "darkorchid1", "chartreuse4")
names(ct.colours) <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelia")

sig.colours <- rev(kelly()[2:8])
names(sig.colours) <- c("F5", "LK", "MM", "IP", "SC", "CA", "MultiBrain")

reg.colours <- c("darkmagenta", "firebrick2", "deeppink", "springgreen4")
names(reg.colours) <- c("CTX", "CB", "sCTX", "SP")

## Plot goodness of fit
plot.gof <- function(list, stat) {
    plot.data <- sapply(list, function(x) {x$stats[[stat]]})
    plot.data <- melt(plot.data)
    
    ggplot(plot.data, aes(x = Var2, y = value)) +
      geom_boxplot() +
      theme_bw() +
      labs(y = stat) +
      theme(axis.title.x = element_blank())
}

## Plot an empty space
plot.empty <- qplot(1,1, colour = I(0.01)) + theme_void()

