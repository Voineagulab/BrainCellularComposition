################################################################################################################################ #
## Setup ----

rm(list = ls())
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/rme/")
options(stringsAsFactors = FALSE)

## Load in data 
load("../../Data/Preprocessed/rme.rda")
load("../../Data/Preprocessed/signatures.rda"); rm(sigsBrain)
load("../../Data/Preprocessed/Signatures - Brain.rda")
sigsRME <- c(sigsRME["IH"], sigsBrain)
sigsRME <- lapply(sigsRME, function(x) {
  x <- x[,c("Neurons", "Astrocytes")]
  x[which(apply(x, 1, max) > 1),]
})

## Functions & Packages
source("../../Scripts/Fun_Composition.R")

################################################################################################################################ #
## Setup lists  ----

## Estimates' list
stats <- list(DRS = list(),
              CIB = list(),
              DTA = list(),
              xCell = list(),
              Blender = list(),
              Linseed = list())

## Stats list
stats <- est

################################################################################################################################ #
## Estimate composition  ----

## Run algorithms with custom signatures
for (j in names(sigsRME)) { 
  print(paste0(j, ":", Sys.time()))
  est$DRS[[j]] <- run.DRS(rme, sigsRME[[j]])
  est$CIB[[j]] <- run.CIB(from.file = FALSE,
                          sigObject = sigsRME[[j]],
                          mixString = "rme.txt")
  est$DTA[[j]] <- run.DTA(rme, sigsRME[[j]], alg = "diff", q = 0.01)
}

  # # CIB / LK fails, as the algorithm cannot find a "nu" above its threshold. below, I manually set it to a zero matrix for compatibility with scripts
  # est$CIB$LK <- est$CIB$F5; est$CIB$LK[,] <- 0

## Run algorithms with in-built signatures
est$xCell <- run.xCell(symbols)
est$Blender <- run.Blender(symbols)

## Run full deconvolution by Linseed
  est$Linseed <- linseed <- list()   

  # on all five samples
  pdf(file = "Linseed Plots.pdf", height = 2.5, width = 2.5)
  linseed$Full <- run.linseed(mixture = rme, nCelltypes = 2, write.plots = TRUE, write.data = TRUE)
  dev.off()
  
  est$Linseed$Full <- linseed$Full$Transformed
  colnames(est$Linseed$Full) <- c("Neurons", "Astrocytes") # manual relabelling of cell-types based on correlation
  
  pdf(file = "Linseed Plots (SVD).pdf", height = 2.5, width = 4)
  linseed$Full$Data$svdPlot() + labs(y = "Cumulative Variance Explained", x = "Number of Dimensions")
  dev.off()
  
  # on the middle three samples
  pdf(file = "Linseed Plots (Mixed Samples).pdf", height = 4, width = 4)
  linseed$Mixed <- run.linseed(mixture = rme[,2:4], nCelltypes = 2, write.plots = TRUE, write.data = TRUE)
  est$Linseed$Mixed <- linseed$Mixed$Transformed
  dev.off()
  
  colnames(est$Linseed$Mixed) <- c("Astrocytes", "Neurons") # manual relabelling of cell-types based on correlation
  est$Linseed$Mixed <- est$Linseed$Mixed[,c("Neurons", "Astrocytes")]
  
  pdf(file = "Linseed Plots (Mixed Samples' SVD).pdf", height = 2.5, width = 4)
  linseed$Mixed$Data$svdPlot() + labs(y = "Cumulative Variance Explained", x = "Number of Dimensions")
  dev.off()
  
  
## Save
save(est, file = "RME Composition Estimates (Revised).rda") # load("RME Composition Estimates.rda")
save(linseed, file = "Raw Linseed Data.rda")

# Filter to relevant cell-types in enrichment algorithms
est$xCell <- lapply(est$xCell, function(x) x[,c("Neurons", "Astrocytes")])
est$Blender <- lapply(est$Blender, function(x) x[,c("Neurons", "Astrocytes")])


################################################################################################################################ #
## Statistics  ----

## Stats
  # compute for deconvolution algorithms
  for (j in c("DRS", "DTA", "CIB")) {
    stats[[j]] <- lapply(est[[j]], function (x) { write.stats(true_RME, x, alg = j, error = TRUE) } )
  }
  
  # compute for xCell
  stats$xCell <- list()
  stats$xCell$Raw <- write.stats(true_RME, est$xCell$Raw, alg = "xCell.Raw", error = FALSE) 
  stats$xCell$Trans <- write.stats(true_RME, est$xCell$Transformed, alg = "xCell.Trans", error = FALSE) 
  
  # compute for Blender
  stats$Blender <- list()
  stats$Blender$AverageIndex <- write.stats(true_RME, est$Blender$AverageIndex, alg = "Blender.Average", error = FALSE) 
  stats$Blender$DarmanisIndex <- write.stats(true_RME, est$Blender$DarmanisIndex, alg = "Blender.Darmanis", error = FALSE) 
  
  # compute for Linseed  (full version)
  stats$Linseed <- list()
  stats$Linseed$Raw <- write.stats(true_RME, est$Linseed$Full$Raw, alg = "Linseed.Raw", error = TRUE) 
  stats$Linseed$Transformed <- write.stats(true_RME, est$Linseed$Full$Transformed, alg = "Linseed.Trans", error = TRUE) 
  
  # save
  save(stats, file = "Statistics (Revised).rda")
  
  # condense lists to single matrices (for easy export)
  for (j in names(stats)) {
    stats[[j]] <- do.call("cbind", data.frame(stats[[j]]))
    rownames(stats[[j]]) <- c("rho", "r", "rmse", "nrmse", "mae", "nmae")
  }
  
  for(j in names(stats)) write.csv(stats[[j]], file = paste0("Statistics - ", j, ".csv"))


################################################################################################################################ #
## Plot using matched / default signature ----  
  
## Partial deconvolution algorithms
  plot.list <- list()
  plot.list$CIB <- plot.scatter(t = true_RME, e = est$CIB$IH, ct = "Neurons", calcCor = "r", calcError = "nmae", colour = ct.colours[["Neurons"]], abline.colour = "black") + 
    labs(y = "Estimated NP", x = "True NP") +
    annotate("text", x = 0.88, y = 0.1, label = "CIB")
  plot.list$DRS <- plot.scatter(t = true_RME, e = est$DRS$IH, ct = "Neurons", calcCor = "r", calcError = "nmae", colour = ct.colours[["Neurons"]], abline.colour = "black") + 
    labs(y = "Estimated NP", x = "True NP") + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
    annotate("text", x = 0.88, y = 0.1, label = "DRS")
  plot.list$dtangle <- plot.scatter(t = true_RME, e = est$dtangle$IH, ct = "Neurons", calcCor = "r", calcError = "nmae", colour = ct.colours[["Neurons"]], abline.colour = "black") + 
    labs(y = "Estimated NP", x = "True NP") + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
    annotate("text", x = 0.82, y = 0.1, label = "dtangle")
  # plot.list$blank <- plot.empty
  pdf(file = "Signature-optimised Scatterplots, Deconvolution.pdf", height = 2, width = 4.5)
  plot_grid(plotlist = plot.list, ncol = 3, rel_widths = c(1, 0.8, 0.8))
  dev.off()
  
## Blender and xCell
  plot.list <- list()
  plot.list$bn <- plot.scatter(t = true_RME, e = est$Blender$AverageIndex, ct = "Neurons", calcCor = "r", calcError = FALSE, colour = ct.colours[["Neurons"]], abline = FALSE) + 
    scale_y_continuous(limits = c(NA, NA)) + 
    labs(x = "True NP", y = "Neuronal Enrichment") +
    annotate("text", x = 0.82, y = -0.5, label = "Blender")
  plot.list$xn <- plot.scatter(t = true_RME, e = est$xCell$Transformed, ct = "Neurons", calcCor = "r", calcError = FALSE, colour = ct.colours[["Neurons"]], abline = FALSE) + 
    labs(x = "True NP") +
    annotate("text", x = 0.88, y = 0.005, label = "xCell") +
    theme(axis.title.y = element_blank())
  
  plot.list$ba <- plot.scatter(t = true_RME, e = est$Blender$AverageIndex, ct = "Astrocytes", calcCor = "r", calcError = FALSE, colour = ct.colours[["Astrocytes"]], abline = FALSE) + 
    scale_y_continuous(limits = c(NA, NA)) + 
    scale_y_continuous(limits = c(NA, NA)) + 
    labs(x = "True AP", y = "Astrocytic Enrichment") +
    annotate("text", x = 0.25, y = -0.8, label = "Blender")
  plot.list$xa <- plot.scatter(t = true_RME, e = est$xCell$Transformed, ct = "Astrocytes", calcCor = "r", calcError = FALSE, colour = ct.colours[["Astrocytes"]], abline = FALSE) + 
    labs(x = "True AP") +
    annotate("text", x = 0.88, y = 1e-18, label = "xCell") +
    theme(axis.title.y = element_blank()) +
    scale_y_continuous(breaks = c(0, 5e-18, 1e-17, 1.5e-17, 2e-17), limits = c(0,2.1e-17))
  
  pdf(file = "Signature-optimised Scatterplots, Enrichment.pdf", height = 2, width = 7.5)
  plot_grid(plotlist = plot.list, ncol = 4, rel_widths = c(1, 0.9, 1, 0.9))
  dev.off()
  
## Linseed
  ## Full run
pdf(file = "Linseed Scatterplots.pdf", height = 3.5, width = 3.5)
plot.scatter(t = true_RME, e = est$Linseed$Full$Transformed, ct = "Neurons", calcCor = "r", calcError = "nmae", colour = ct.colours[["Neurons"]]) +
  labs(x = "True NP", y = "Linseed Cell-type 1")
dev.off()

  ## The run on the mixed samples 
  pdf(file = "Linseed Scatterplots (Mixed Samples).pdf", height = 3.5, width = 3.5)
  plot.scatter(t = true_RME[2:4,], e = est$Linseed$Mixed$Transformed, ct = "Neurons", calcCor = "r", calcError = "nmae", colour = ct.colours[["Neurons"]]) +
    labs(x = "True NP", y = "Linseed Cell-type 1") +
    scale_y_continuous(limits = c(0.3, 0.55)) +
    scale_x_continuous(limits = c(0.3, 0.55)) 
  dev.off()

  
################################################################################################################################ #
## Plot deconvolution using mismatched signatures  ----  

## Scatterplot   
  plot.list <- list()
  for(k in c("CIB", "DRS", "DTA")) {
    for(j in names(sigsRME)) {
      
      annot.pos <- c(Inf, -Inf, 1, -0.5) 
      e <- est[[k]][[j]]
      
      plot.list[[j]] <- plot.scatter(t = true_RME, e = e, ct = "Neurons", ylab = "Estimated Proportion", calcCor = "r", 
                                     calcError = "nmae", colour = "black", abline.colour = "black", abline = TRUE, annot.pos = annot.pos)  +
        labs(title = j, x = "True Proportion") +
        scale_y_continuous(limits = c(0,1)) +
        scale_x_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) +
        theme(plot.title = element_text(hjust = 0.5), panel.grid.minor = element_blank())
      if (!(j %in% c("IH", "VL"))) plot.list[[j]] <- plot.list[[j]] + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
      if (!(j %in% names(sigsRME)[6:10])) plot.list[[j]] <- plot.list[[j]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
      if (j == "IH") plot.list[[j]] <- plot.list[[j]] + labs(title = "Matched (IH)")
    }
    
    pdf(file = paste0("Origin Test - ", k, ", Scatterplot (Revised).pdf"), height = 4.5, width = 8)
    print(plot_grid(plotlist = plot.list, ncol = 5, rel_widths = c(1.2,1,1,1,1), rel_heights = c(1,1.2)))
    dev.off()
    
  }
 

############################################################# FIN ################################################################