################################################################################################################################ #
## Setup ----

rm(list = ls())
setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Results/Revisions/scme_revised/")
options(stringsAsFactors = FALSE)

## Functions and packages
  source("../../../Scripts/Fun_Composition.R")
  source("../../../Scripts/Fun_Preprocessing.R")

## Load in scme data
  load("../../../Data/Preprocessed/scme.rda")
  
  ## Reprocess
    scme <- rpkm(counts = scme$mixture$counts)

  ## Find truth
    true <- scme.v2$true$cell_type_s
    
## Load in the gradient
  load("../../../Data/Preprocessed/scme.gradients.rda")

## Load in additional signatures
  load("../../../Data/Preprocessed/Signatures - Brain.rda")
  load("../../../Data/Preprocessed/Signatures - Brain (MuSiC).rda")
  
  ## Filter to just key cell-types
    sigsBrain <- lapply(sigsBrain, function(x) {
      x <- x[,which(colnames(x) %in% colnames(true))]
      return(x)
    })
  
  ## Remove endothelia from LK and CA (for consistency with other runs)
    sigsBrain$CA <- sigsBrain$CA[,-grep("End", colnames(sigsBrain$CA))]
    sigsBrain$LK <- sigsBrain$LK[,-grep("End", colnames(sigsBrain$LK))]

  ## Reprocess the Music signatures
    # rename a column in DM
    sigsMuSiC$DM$meta$brain.ct2 <- as.character(sigsMuSiC$DM$meta$orig.celltype)
    
    # remove OPCs
    sigsMuSiC <- lapply(sigsMuSiC, function(x) {
      remove <- grep("OPC", x$meta$brain.ct2)
      
      x$counts <- x$counts[,-remove]
      x$meta <- x$meta[-remove,]
      
      return(x)
    })
    
    # remove end from LK and CA
    remove <- grep("End", sigsMuSiC$CA$meta$brain.ct2)
    sigsMuSiC$CA$counts <- sigsMuSiC$CA$counts[,-remove]
    sigsMuSiC$CA$meta <- sigsMuSiC$CA$meta[-remove,]    
    
    remove <- grep("End", sigsMuSiC$LK$meta$brain.ct2)
    sigsMuSiC$LK$counts <- sigsMuSiC$LK$counts[,-remove]
    sigsMuSiC$LK$meta <- sigsMuSiC$LK$meta[-remove,]    
    
## Other
    invis <- element_blank()
    deconv <- c("CIB", "DRS", "DTA", "MUS")

################################################################################################################################ #
## Setup lists  ----
  
est <- list(DRS = list(),
            CIB = list(),
            DTA = list(),
            MUS = list(),
            xCell = list(),
            Blender = list(),
            Linseed = list(),
            Coex = list())

est.gradient <- list(Linseed = list(),
                     Coex = list())

## Stats list
stats <- est

################################################################################################################################ #
## Estimate composition 

## Estimate composition by partial deconvolution methods
  for (j in names(sigsBrain)) {
    # timestamp
    print(paste0(j, ": ", Sys.time()))
    
    # run deconrnaseq
    est$DRS[[j]] <- run.DRS(mixture = scme, 
                            signature = sigsBrain[[j]])
    
    # run cibersort
    est$CIB[[j]] <- run.CIB(from.file = FALSE,
                            sigObject = sigsBrain[[j]],
                            mixString = "scme.v2.txt")
        
    # run dtangle
    est$DTA[[j]] <- run.DTA(mixture = scme, 
                                    signature = sigsBrain[[j]],
                                    alg = "diff",
                                    q = 0.01)
    
    # run music
    if (j %in% names(sigsMuSiC)) {
      est$MUS[[j]] <- run.music(mixture = scme.v2$mixture$counts, 
                                    use.meta.column = "brain.ct2", 
                                    signature = sigsMuSiC[[j]]$counts, 
                                    signature.meta = sigsMuSiC[[j]]$meta) 
      
    }
    
  }

  save(est, file = "Composition Estimates.rda")
  
  

## Estimate composition by enrichment-based algorithms
  est$xCell <- run.xCell(addSymbol(scme))
  est$Blender <- run.Blender(addSymbol(scme))
  
  
################################################################################################################################ #
## Complete deconvolution ----
 
complete <- list()
  

## Linseed function
  apply.linseed <- function(mixture, title) {
    n <- 5
    
    pdf(file = paste0("Linseed Plots (", title, ").pdf"), height = 2.5, width = 5)
    output <- run.linseed(mixture = mixture,
                          nCelltypes = n,
                          write.plots = TRUE,
                          write.data = TRUE) 
    output$Data$svdPlot() + 
      labs(y = "Cumulative Variance Explained", x = "Number of Dimensions")
    dev.off()
    
    output <- output$Transformed  
    return(output)
  }
  
## A coex function
  analyse.coex <- function(coex) {
    # enrichment tests for ct assignment
    coex <- define.celltypes(x = coex, algorithm.type = "Coex", celltype.signature = sigsSCME, 
                             return.confidence = TRUE, nMarks.x = 10000, nMarks.sig = 100, bg = rownames(scme))
    
    # reannotate     
    m <- match(colnames(coex$comp), names(coex$assignments))
    colnames(coex$comp) <- coex$assignments[m]
    coex <- coex$comp
    
    # return
    return(coex)
  }
    
## Linseed: on the standard mixtures
  complete$linseed <- apply.linseed(scme, "")

## Linseed: on the gradient mixtures
  complete$linseed.gB <- apply.linseed(scme.gradients$confound.p, "gB")
  
    

## Coex: on the standard mixtures
    coex.standardOutput <- run.coex(mixture = scme,
                                  signature = sigsSCME,
                                  only.threshold = FALSE,
                                  sft = "auto", 
                                  output.all = TRUE)
    
    complete$coex <- analyse.coex(coex.standardOutput)

    
## Coex: on the gradient mixtures
    coex.gradientOutput.B <- run.coex(mixture = scme.gradients$confound.p,
                                      signature = sigsSCME,
                                      only.threshold = FALSE,
                                      sft = "auto", 
                                      output.all = TRUE)
    ## Analyse
    complete$coex.gB <- analyse.coex(coex.gradientOutput.B)
    
  
## Save
  save(coex.gradientOutput.B, coex.standardOutput, file = "Coexpression Data.rda")
  save(complete, file = "Complete Deconvolution Estimates.rda") 
  
  
################################################################################################################################ #
## Analyse randomly sampled mixtures ----

## Statistics
  stats <- list()
  
  for (j in deconv) { stats[[j]] <- write.stats(true, est[[j]]$DM, alg = j, error = TRUE) }
  stats$xCell <- write.stats(true, est$xCell, alg = "xCell", error = FALSE) 
  stats$Blender <- write.stats(true, est$Blender, alg = "Blender", error = FALSE) 


## Scatterplot
  # deconvolution algorithms using the DM signature
    # setup list of plots
    plot.list <- list()
    for (j in c(deconv)) { for (k in colnames(true)[c(4,1,5,3,2)]) {
        run <- paste0(j, k)
        plot.list[[run]] <- plot.scatter(t = true, e = est[[j]]$DM, abline.colour = "grey50", ct = k, calcCor = "r", calcError = "nmae", colour = ct.colours[[k]])  +
          theme(axis.title.x = invis) 
        if (j == "CIB") plot.list[[run]] <- plot.list[[run]] + labs(title = substr(k, 1, 3)) + theme(plot.title = element_text(hjust = 0.5, colour = ct.colours[[k]])) 
        if (k == "Neurons") plot.list[[run]] <- plot.list[[run]] + labs(y = paste0(j, "\nProportion")) else plot.list[[run]] <- plot.list[[run]] + theme(axis.title.y = invis)
        if (k == "Oligodendrocytes") plot.list[[run]] <- plot.list[[run]] + scale_x_continuous(breaks = c(0, 0.1, 0.2), limits = c(0, NA)) 
        if (k == "Endothelia") plot.list[[run]] <- plot.list[[run]] + scale_x_continuous(breaks = c(0, 0.05, 0.1), limits = c(0, NA)) 
        if (k == "Microglia") plot.list[[run]] <- plot.list[[run]] + scale_x_continuous(breaks = c(0, 0.04, 0.08), limits = c(0, NA))
        
    } }
    
    # plot
    pdf(file = "Random Sampling - Scatterplots Deconvolution.pdf", height = 6, width = 8)
    plot_grid(plotlist = plot.list, ncol = 5, nrow = 4, rel_widths = c(1.3, 1, 1, 1, 1), rel_heights = c(1, 0.8, 0.8, 0.8))
    dev.off()
  
  # blender
    # setup list of plots, but only for the "Average Index
    plot.list2 <- list()
    for (k in colnames(est$Blender)[c(2,1,3,5,4)]) {
        plot.list2[[k]] <- plot.scatter(t = true, e = est$Blender, ct = k, abline.colour = "grey50",
                 calcCor = "r", calcError = FALSE, colour = ct.colours[[k]], ylab = "Celltype Index", limits = c(NA, NA), abline = FALSE) +
          labs(x = "True Proportion")
        if (k == "Neurons") {
          plot.list2[[k]] <- plot.list2[[k]] + labs(y = "Blender\nScore") 
            } else {
          plot.list2[[k]] <- plot.list2[[k]] + theme(axis.title.y = invis) 
        }
    } 
    
    plot.list2$Endothelia <- plot.list2$Endothelia + scale_x_continuous(breaks = c(0, 0.05, 0.1), limits = c(0, 0.14))
    plot.list2$Microglia <- plot.list2$Microglia + scale_x_continuous(breaks = c(0, 0.04, 0.08), limits = c(0, NA))
    plot.list2$Oligodendrocytes <- plot.list2$Oligodendrocytes + scale_x_continuous(breaks = c(0, 0.1, 0.2), limits = c(0, NA)) 
    plot.list2$Neurons <- plot.list2$Neurons + scale_x_continuous(breaks = c(0.4, 0.48, 0.56))
    plot.list2$Astrocytes <- plot.list2$Astrocytes + scale_x_continuous(breaks = c(0.1, 0.2, 0.3))
    
    
    
    # plots
    pdf(file = "Random Sampling - Scatterplots Blender.pdf", height = 1.6, width = 8)
    plot_grid(plotlist = plot.list2, ncol = 5, rel_widths = c(1.3,1,1,1,1))
    dev.off()
  
  
  # xCell
    # setup list of plots
    plot.list3 <- list()
    for (k in c("Neurons", "Astrocytes")) {
      plot.list3[[k]] <- plot.scatter(t = true, e = signif(est$xCell,2), ct = k, calcCor = "r", 
                                   calcError = FALSE, colour = ct.colours[[k]], limits = c(NA, NA), abline = FALSE, abline.colour = "grey50") + 
          theme(plot.title = element_text(hjust = 0.5, colour = ct.colours[[k]]), axis.text.y = element_text(size = 6.7, angle = 45)) +
          scale_y_continuous(limits = c(0, NA)) + labs(x = "True Proportion")
        if (k == "Neurons") plot.list3[[k]] <- plot.list3[[k]] + labs(y = "xCell\nScore") else plot.list3[[k]] <- plot.list3[[k]] + theme(axis.title.y = invis)
    }
    plot.list3$Astrocytes <- plot.list3$Astrocytes + scale_y_continuous(breaks = c(0, 1e-17, 2e-17, 3e-17), limits = c(0, 3.2e-17))
    plot.list3$blank3 <- plot.list3$blank2 <- plot.list3$blank1 <- plot.empty
    
    # plot
    pdf(file = "Random Sampling - Scatterplots xCell.pdf", height = 1.5, width = 8)
    plot_grid(plotlist = plot.list3, ncol = 5, nrow = 1, rel_widths = c(1.5,1.2,1,1,1))
    dev.off()
    
## Plot error
  # setup
  plot.data <- melt(stats)
  colnames(plot.data) <- c("Metric", "Celltype", "Value", "Algorithm")

  
  # reannotate
  plot.data$Metric <- sapply(strsplit(as.character(plot.data$Metric), "_"), "[", 1)
  plot.data$Algorithm <- factor(plot.data$Algorithm, levels = c("CIB", "DRS", "DTA", "MUS", "Blender", "xCell"))

  # plot r
  new.data <- data.frame(Metric = "r", Celltype = c("Oligodendrocytes", "Microglia", "Endothelia"), Value = NA, Algorithm = "xCell")
  plot.data <- rbind(plot.data, new.data)

  pA <- ggplot(plot.data[which(plot.data$Metric == "r"),], aes(x = Algorithm, y = Value, fill = Celltype, colour = Celltype)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.1) +
    geom_point(position = position_dodge(width = 0.8), size = 2) +
    theme_bw() +
    scale_fill_manual(values = ct.colours, labels = c("Neu", "Ast", "Oli", "Mic", "End")) +
    scale_colour_manual(values = ct.colours, labels = c("Neu", "Ast", "Oli", "Mic", "End")) +
    theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis, axis.title.x = invis) +
    geom_hline(yintercept = 0) +
    labs(y = "r")

  # plot normalised error?
  keep <- which(plot.data$Algorithm %in% c("CIB", "DRS", "DTA", "MUS") & plot.data$Metric == "nmae")

  pB <- ggplot(plot.data[keep,], aes(x = Algorithm, y = Value, fill = Celltype, colour = Celltype)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.1) +
    geom_point(position = position_dodge(width = 0.8), size = 2) +
    theme_bw() +
    scale_fill_manual(values = ct.colours) +
    scale_colour_manual(values = ct.colours, labels = c("Neu", "Ast", "Oli", "Mic", "End")) +
    theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis, legend.position = "none", axis.title.x = invis) +
    geom_hline(yintercept = 0) +
    labs(y = "nmae")


  # combined plot
  pdf(file = "Random Sampling - Barplots.pdf", height = 1.5, width = 8)
  plot_grid(pB, pA, rel_widths = c(1, 1.7))
  dev.off()

  
################################################################################################################################ #
## The effect of signature choice ----

## Append data
  ot <- est[deconv]

  ## Alternate plotting scheme - all cts on a single plot
  new.plot <- function(x, alg) {
    plot.data <- list()
    for(j in colnames(ot[[alg]][[x]])) plot.data[[j]] <- data.frame(Est = ot[[alg]][[x]][,j], True = true[,j], ct = j)
    plot.data <- do.call("rbind", plot.data)
    
    ggplot(plot.data, aes(x = True*100, y = Est*100, col = ct)) +
      geom_point(alpha = 0.5) +
      scale_colour_manual(values = ct.colours) +
      theme_bw() +
      annotate("text", x = 15, y = 95, label = x) +
      # geom_smooth(method = "lm", se = FALSE) +
      theme(panel.border = invis, axis.line = element_line(), plot.title = element_text(hjust = 0.5)) +
      geom_abline(slope = 1, intercept = 0, col = "red", linetype = 2) +
      scale_y_continuous(expand = c(0,0), limits = c(-0.005, 101), breaks = c(0, 25, 50, 75)) +
      scale_x_continuous(expand = c(0,0), limits = c(-0.005, 101), breaks = c(0, 25, 50, 75)) +
      labs(y = "Estimated Proportion", x = "True Proportion") 
  }
  
  for (j in names(ot)) {
    # original publication plot
    plot.list <- list(DM = new.plot("DM", j) + theme(legend.position = "none", axis.title.x = invis) + labs(y = ""),
                      NG = new.plot("NG", j) + theme(legend.position = "none", axis.title = invis, axis.text.y = invis),
                      CA = new.plot("CA", j) + theme(legend.position = "none", axis.title = invis, axis.text.y = invis),
                      LK = new.plot("LK", j) + theme(legend.position = "none", axis.title.x = invis),
                      TS = new.plot("TS", j) + theme(legend.position = "none", axis.title = invis, axis.text.y = invis), 
                      F5 = new.plot("F5", j) + theme(legend.position = "none", axis.title = invis, axis.text.y = invis),
                      IP = new.plot("IP", j) + theme(legend.position = "none")  + labs(y = "", x = ""),
                      MM = new.plot("MM", j) + theme(legend.position = "none", axis.title.y = invis, axis.text.y = invis),
                      VL = new.plot("VL", j) + theme(legend.position = "none", axis.title.y = invis, axis.text.y = invis)  + labs(x = ""))
  
    if (j == "MUS") { 
      plot.list <- plot.list[which(names(plot.list) %in% names(est$MUS))]
      plot.list[4:6] <- lapply(plot.list[4:6], function(x) {
        x + theme(axis.title.x = element_text()) 
      })
      plot.list$DM <- plot.list$DM+ labs(y = "Estimated Proportion")
      
      pdf(file = paste0("Mismatched Signatures - Scatterplots, ", j, ".pdf"), height = 3, width = 3.8)
      print(plot_grid(plotlist = plot.list, ncol = 3, nrow = 2, rel_widths = c(1.2,1,1), rel_heights = c(1,1.2)))
      dev.off()    
      
    } else {
      pdf(file = paste0("Mismatched Signatures - Scatterplots, ", j, ".pdf"), height = 4.5, width = 3.8)
      print(plot_grid(plotlist = plot.list, ncol = 3, nrow = 3, rel_widths = c(1.2,1,1), rel_heights = c(1,1,1.2)))
      dev.off()    
    }
  }
  
## Statistics
  # calculate
  stats_ot <- lapply(ot, function(x) {
    y <- list()
    for(j in names(x)) {
      y[[j]] <- write.stats(t = true, e = x[[j]], alg = "", error = TRUE)
    }
    return(y)
  })

  # barplots
    # setup
    plot.data <- melt(stats_ot)
    
    # rename columns
    colnames(plot.data) <- c("Metric", "Celltype", "Value", "Signature", "Algorithm")
    
    
    # relabel factors
    plot.data$Signature <- factor(plot.data$Signature, levels = c("DM", "NG", "CA", "LK", "TS", "F5", "IP", "MM", "VL"))
    # levels(plot.data$Signature)[7] <- "Merged"
    plot.data$Metric <- gsub("_", "", plot.data$Metric)
    
    # filter
    plot.data <- plot.data[which(plot.data$Metric %in% c("nmae", "r")),]
    
    
    pdf(file = "Mismatched Signatures - Barplots.pdf", height = 4, width = 8) 
    ggplot(plot.data[which(plot.data$Metric == "r"),], aes(x = Signature, y = Value, fill = Celltype, colour = Celltype)) +
      # geom_col(position = position_dodge2(preserve = "single"), colour = "black", width = 0.8) +
      geom_col(position = position_dodge(width = 0.8), width = 0.1) +
      geom_point(position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = ct.colours, labels = c("Neu", "Ast", "Oli", "Mic", "End")) +
      scale_colour_manual(values = ct.colours, labels = c("Neu", "Ast", "Oli", "Mic", "End")) +
      # facet_grid(Metric~Algorithm, scales = "free_y") +
      facet_wrap(~Algorithm, nrow = 2) +
      theme_bw() +
      labs(y = "Correlation (r)") +
      scale_y_continuous(expand = c(0,0), limits = c(0, 1.01)) +
      theme(axis.title.x = invis, legend.title = invis, axis.line = element_line())
    
    ggplot(plot.data[which(plot.data$Metric == "nmae"),], aes(x = Signature, y = Value, fill = Celltype, colour = Celltype)) +
      # geom_col(position = position_dodge2(preserve = "single"), colour = "black", width = 0.8) +
      geom_col(position = position_dodge(width = 0.8), width = 0.1) +
      geom_point(position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = ct.colours, labels = c("Neu", "Ast", "Oli", "Mic", "End")) +
      scale_colour_manual(values = ct.colours, labels = c("Neu", "Ast", "Oli", "Mic", "End")) +
      # facet_grid(Metric~Algorithm, scales = "free_y") +
      facet_wrap(~Algorithm, nrow = 2) +
      theme_bw() +
      labs(y = "nmae") +
      scale_y_continuous(expand = c(0,0), limits = c(0, 8), trans = pseudo_log_trans(sigma = 0.1), breaks = c(0.1, 0.5, 1, 2, 3, 4, 6)) +
      theme(axis.title.x = invis, legend.title = invis, axis.line = element_line())
    dev.off()

## Goodness of fit across signatures!
  # here, I am calculating gof for exc/inh subtypes (where possible), for CIB-based deconvolutions
  gof <- list()
    
  # function
  get.gof <- function(e, sig) {
    sig <- sig[,colnames(e)]
    output <- write.gof.v2(estimatedComp = e, measuredExp = scme, signatureUsed = sig)
    output <- output$r
  }
  
  # calculate
  gof$VL <- get.gof(est$CIB$VL, sigsBrain$VL)
  gof$CA <- get.gof(est$CIB$CA, sigsBrain$CA)
  gof$TS <- get.gof(est$CIB$TS, sigsBrain$TS)
  gof$LK <- get.gof(est$CIB$LK, sigsBrain$LK)
  gof$NG <- get.gof(est$CIB$NG, sigsBrain$NG)
  gof$F5 <- get.gof(est$CIB$F5, sigsBrain$F5)
  gof$IP <- get.gof(est$CIB$IP, sigsBrain$IP)
  gof$MM <- get.gof(est$CIB$MM, sigsBrain$MM)
  gof$DM <- get.gof(est$CIB$DM, sigsBrain$DM)
  
  # generate plot data
  plot.data <- do.call("cbind", gof)
  plot.data <- melt(plot.data)
  order <- names(gof)[rev(order(sapply(gof, mean)))]
  plot.data$Var2 <- factor(plot.data$Var2, levels = order)
  
  
  
  # and plot! 
  pdf(file = "GoF.pdf", height = 3, width = 4)
  ggplot(plot.data, aes(x = Var2, y = value)) +
    # geom_violin(scale = "width", width = 0.9, fill = "white", colour = "black") +
    # geom_boxplot(fill = "grey50", colour = "black") +
    geom_point(position = position_jitter(), alpha = 0.2, size = 0.5) +
    # stat_summary(geom = "point", fun = function(x) mean(x), shape = "-", size = 8, show.legend = FALSE, colour = "firebrick1") +
    theme_bw() +
    theme(panel.grid = invis, panel.border = invis, axis.line = element_line()) +
    geom_hline(yintercept = c(0.5, 0.7), linetype = 2) +
    labs(y = "Goodness of Fit (r)", x = "Signature")
  dev.off()

  
################################################################################################################################ #
## Formal plots of complete deconvolution ----  
  
  
## On Linseed
  ## Heatmap function
    linseed.heatmap <- function(e, t, filename) {
      plot.data <- cor(e, t, method = "p")
      
      plot.data <- as.data.frame(plot.data)
      plot.data$Linseed <- LETTERS[1:nrow(plot.data)]
      plot.data <- melt(plot.data)
      colnames(plot.data)[3] <- "Correlation"
      
      pdf(file = filename, height = 2.7, width = 3)
      print(ggplot(plot.data, aes(x = variable, y = Linseed, fill = Correlation, label = round(Correlation, 2))) +
              geom_tile(colour = "black") +
              geom_text(size = 2.5) +
              theme_bw() +
              scale_fill_carto_c(palette = "ArmyRose", limits = c(-1,1)) +
              scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
              labs(x = "True Cell-type") +
              theme(panel.border = invis, axis.ticks = invis, panel.grid = invis, axis.title.y = invis, axis.text.x = element_text(angle = 90, hjust = 1)))
      dev.off()
    }  
    
  ## Apply
    linseed.heatmap(complete$linseed, true, "Complete Linseed Heatmap.pdf")
    linseed.heatmap(complete$linseed.gB, scme.gradients$meta.prop, "Complete Linseed Heatmap, gB.pdf")

    
  ## Scatterplot function
    complete.scatter <- function(t, l, ct, ylab, axis.min = 0, abline = TRUE) {
      r <- round(cor(t, l), 2)
      label <- paste0("r=", r)
      p <- qplot(t, l, colour = I(ct.colours[ct])) +
        theme_bw() +
        # geom_smooth(method = "lm", se = FALSE) +
        scale_y_continuous(limits = c(axis.min, NA), expand = c(0,0)) +
        labs(y = ylab, x = paste0("True ", ct, "\n", label)) +
        scale_x_continuous(limits = c(axis.min, NA), expand = c(0,0)) +
        theme(panel.border = invis, axis.line = element_line()) +
        NoLegend()
      
      if (abline) {
        return(p + geom_abline(slope = 1, intercept = 0, linetype = 2))
      } else {
        return(p)
      }
    }
    
  ## Apply
    pdf(file = "Complete Linseed Scatter.pdf", height = 2, width = 2)
    complete.scatter(true$Endothelia, complete$linseed$`Cell type 3`, ct = "Endothelia", ylab = "Linseed Celltype C")
    complete.scatter(true$Microglia, complete$linseed$`Cell type 5`, ct = "Microglia", ylab = "Linseed Celltype E") 
    dev.off()
    
    pdf(file = "Complete Linseed Scatter, gB.pdf", height = 2, width = 2)
    complete.scatter(scme.gradients$meta.prop$Endothelia, complete$linseed.gB$`Cell type 2`, ct = "Endothelia", ylab = "Linseed Celltype B")
    complete.scatter(scme.gradients$meta.prop$Microglia, complete$linseed.gB$`Cell type 5`, ct = "Microglia", ylab = "Linseed Celltype E") 
    dev.off()
    
    
    
## Coex
  ## Barplot
    coex.barplot <- function(e, u, filename) {
      dat <- diag(cor(e, u[,colnames(e)]))
      
      absent <- colnames(u)[-which(colnames(u) %in% colnames(e))]
      zeros <- rep(0, length(absent))
      names(zeros) <- absent
      dat <- c(dat, zeros)
      
      dat <- data.frame(Cor = dat, Celltype = names(dat))
      
      pdf(file = filename, height = 2.5, width = 2.5)
      print(ggplot(dat, aes(x = Celltype, y = Cor, fill = Celltype, colour = Celltype)) +
        geom_col() +
        theme_bw() +
        scale_fill_manual(values = ct.colours) +
        scale_colour_manual(values = ct.colours) +
        scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
        scale_y_continuous(limits = c(NA, 1)) +
        labs(y = "Pearson Correlation") +
          geom_hline(yintercept = 0.8, linetype = 2) +
        theme(panel.border = invis, axis.line = element_line(), legend.position = "none"))
      dev.off()
    }
    
    coex.barplot(e = complete$coex, u = true, filename = "Complete Coex Barplot.pdf")
    coex.barplot(e = complete$coex.gB, u = scme.gradients$meta.prop, filename = "Complete Coex Barplot, gradientB.pdf")
 
  
  ## Heatmap
    coex.heatmap <- function(e, u, filename) {
      dat <- diag(cor(e, u[,colnames(e)]))
      
      absent <- colnames(u)[-which(colnames(u) %in% colnames(e))]
      zeros <- rep(0, length(absent))
      names(zeros) <- absent
      dat <- c(dat, zeros)
      
      dat <- data.frame(Cor = dat, Celltype = names(dat), y = "")
      
      pdf(file = filename, height = 1, width = 2.2)
      print(ggplot(dat, aes(x = Celltype, y = y, fill = Cor, label = round(Cor, 2))) +
              geom_tile(colour = "black") +
              geom_text(size = 2.5) +
              theme_bw() +
              scale_fill_carto_c(palette = "ArmyRose", limits = c(-1,1)) +
              scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
              labs(x = "True Cell-type") +
              NoLegend() +
              theme(panel.border = invis, axis.ticks = invis, panel.grid = invis, axis.title.y = invis, axis.text.x = element_text(angle = 90, hjust = 1)))
      dev.off()
    }
    
    coex.heatmap(e = complete$coex, u = true, filename = "Complete Coex Heatmap.pdf")
    coex.heatmap(e = complete$coex.gB, u = scme.gradients$meta.prop, filename = "Complete Coex Heatmap, gradientB.pdf")
    
    
  ## Scatterplot
      pdf(file = "Complete Coex Scatterplot.pdf", height = 2, width = 8)
      plot.list <- lapply(colnames(true), function(x) { complete.scatter(true[,x], complete$coex[,x], ct = x, ylab = "Coex Enrichment", axis.min = NA, abline = FALSE)  })
      plot.list[2:5] <- lapply(plot.list[2:5], function(x) x + theme(axis.title.y = invis))
      plot_grid(plotlist = plot.list, nrow = 1)
      dev.off()
      
      pdf(file = "Complete Coex Scatterplot, gB.pdf", height = 2, width = 8)
      plot.list <- lapply(colnames(true), function(x) { complete.scatter(scme.gradients$meta.prop[,x], complete$coex.gB[,x], ct = x, ylab = "Coex Enrichment", axis.min = NA, abline = FALSE)  })
      plot.list[2:5] <- lapply(plot.list[2:5], function(x) x + theme(axis.title.y = invis))
      plot_grid(plotlist = plot.list, nrow = 1)
      dev.off()
    
  
## Aside: properties of mixtures
  ## Abundance

    pdf(file = "Complete QC Abundance.pdf", height = 1.8, width = 2.5)
    dat <- melt(true)
    ggplot(dat, aes(x = Var2, y = value)) +
      geom_violin(scale = "width", width = 0.8, draw_quantiles = 0.5, colour = "black", fill = "darkorange1") +
      theme_bw() +
      scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
      scale_y_continuous(limits = c(0, 0.85), expand = c(0,0)) +
      theme(panel.grid = invis, axis.title.x = invis, panel.border = invis, axis.line = element_line(), legend.position = "none") +
      labs(y = "Simulated Abundance")

    dev.off()
    
    
    
  ## CorMx
    # random
    qc.corMx <- function(t, filename, flip = FALSE) {
      plot.data <- cor(t, method = "p")
      plot.data <- as.data.frame(plot.data)
      plot.data[lower.tri(plot.data)] <- NA
      plot.data$Ct <- colnames(plot.data) <- substr(rownames(plot.data), 1, 3)
      plot.data <- melt(plot.data, na.rm = TRUE)
      colnames(plot.data)[3] <- "Correlation"
      
      if (flip) {
        plot.data$variable <- factor(plot.data$variable, levels = rev(levels(as.factor(plot.data$Ct)))) 
        
        }
      
      a <- ggplot(plot.data, aes(x = variable, y = Ct, fill = Correlation, label = round(Correlation, 2))) +
              geom_tile(colour = "black") +
              geom_text(size = 2.5) +
              theme_bw() +
              # coord_flip() +
              scale_fill_carto_c(palette = "Earth", limits = c(-1,1)) +
              scale_x_discrete(labels = function(x) substr(x, 1, 3)) +
              scale_y_discrete(position = "right") +
              NoLegend() +
              theme(panel.border = invis, axis.ticks = invis, panel.grid = invis, axis.title = invis, axis.text.x = element_text(angle = 90, hjust = 1))
      
      if (flip) {
        pdf(file = filename, height = 3, width = 2.7)
        print(a + theme(axis.text.y = invis))
        dev.off()
      } else {
        pdf(file = filename, height = 3, width = 2.7)
        print(a)
        dev.off()
      }
    }
    
   qc.corMx(true, "Complete QC Random Cors.pdf")
   qc.corMx(scme.gradients$meta.prop, "Complete QC Non-Random Cors.pdf", flip = TRUE)
   