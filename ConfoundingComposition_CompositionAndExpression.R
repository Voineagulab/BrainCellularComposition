## In this (revised) script, I analyse those simulate mixtures for testing differential expression

## Run ConfoundingComposition_CompositionOnly.R prior to this script. You'll need its models.rda

################################################################################################################################ #
## Generic setup  ----

## Start!
  rm(list = ls()); gc()

## Directory
  root.dir <- "INSERT.YOUR.ROOT.DIRECTORY.PLEASE/"
  setwd(paste0(root.dir, "Results/ConfoundingComposition/CompositionAndExpression/"))

## Functions
source("../../../Scripts/Fun_Composition.R")
source("../../../Scripts/Fun_Preprocessing.R")

## Packages
library(DESeq2)

## Load simulations
load("../../../Data/Preprocessed/ConfoundingComposition_Simulations.rda")

## Load DE results from the previous ConfoundingComposition_CompositionOnly.R script
load("../CompositionOnly/Models.rda")
load("../CompositionOnly/Percent Confounds.rda")
load("../CompositionOnly/Perturbed Gene List.rda")
alpha <- 1:50
beta <- 51:100


  



################################################################################################################################ #
## Visualise false negatives ----


  
## Plotting functions
  ## Function to extract plotting data
    get.tp.data <- function(fc) {
      
      plot.data <- lapply(tests, function(x) {
        # get genes which have an adjusted p-value of < 0.05
        sig <- sapply(x[[fc]], function(y) { p.adjust(y$pvals, method = "BH")  * sign(y$log2fc)  })
        rownames(sig) <- rownames(x$confound.pe.1.5$lm.uncorrected)
        
        
        ## Get information
        # number of false positives
        nFP <- apply(sig, 2, function(y) length(which(abs(y) < 0.05 & !(all.pert))))
        
        # number of true positives
        nTP.rand.dn <- apply(sig, 2, function(y) {
          y <- length(which(abs(y) < 0.05 & rand.dn & y < 0)) / length(which(rand.dn))
        })
        
        nTP.exc.dn <- apply(sig, 2, function(y) {
          y <- length(which(abs(y) < 0.05 & exc.dn & y < 0)) / length(which(exc.dn))
        })
        
        nTP.rand.up <- apply(sig, 2, function(y) {
          y <- length(which(abs(y) < 0.05 & rand.up & y > 0)) / length(which(rand.up))
        })
        
        nTP.exc.up <- apply(sig, 2, function(y) {
          y <- length(which(abs(y) < 0.05 & exc.up & y > 0)) / length(which(exc.up)) 
        })
        
        # number of true positives in the top n most significant genes
        n <- length(which(all.pert)) # as we have perturbed n genes, let's see what percent of the top n hits are truly perturbed!
        
        nTop <- apply(sig, 2, function(y) {
          lowest <- names(y[order(abs(y))])[1:n]
          length(which(lowest %in% pert.genes$EnsID)) / n
        })
        
        rbind(nFP, nTP.exc.up, nTP.exc.dn, nTP.rand.up, nTP.rand.dn, nTop)
      })
      
      # final processing
      plot.data <- melt(plot.data)
      colnames(plot.data) <- c("Metric", "Model", "Value", "Simulation")
      x <- sapply(cf, function(x) x["Excitatory"])
      names(x) <- gsub("\\.Excitatory", "", names(x))
      m <- match(plot.data$Simulation, names(x))
      plot.data$Confound <- x[m]
      
      # return
      return(plot.data)
    }
    

  ## Another reworking... 
    plot.tp.3 <- function(fc) {
      # setup
      filename <- paste0("TP Plot, Style 3 (", fc, ").pdf")
      print(paste0("Saving plot at ", filename))
      x <- plot.data[[fc]]
      x <- x[sample(rownames(x), nrow(x)),] # this makes the overplotting of points less of an issue
      
      x$Model <- factor(x$Model)
      levels(x$Model) <- c("Uncorrected LM", "Corrected LM", "Corrected LM Quadratic", "Corrected Spline", "Uncorrected DESeq2", "Corrected DESeq2", "Corrected DESeq2 Quadratic")
      x$Model <- factor(x$Model, levels = levels(x$Model)[c(1, 5, 2, 6, 3, 7, 4)])
      
      line <- 0.95 * x$Value[which(x$Simulation == "up0.rep1" & x$Model == "Uncorrected LM" & x$Metric == "nTop")]
      
      pdf(file = filename, height = 2.5, width = 3)
      print(ggplot(x[which(x$Metric == "nTop"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
        geom_point(size = 0.5, alpha = 0.5) +
        geom_smooth(se = FALSE, span = 0.3) +
        geom_hline(yintercept =  line, linetype = 2) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_colour_manual(values = model.colours) +
        labs(y = "Fraction of Perturbed Genes\n in 200 Most Significant", x = "Difference in Excitatory Proportion") +
        theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
              legend.background = element_rect(fill = alpha("white", 0)),
              panel.grid.minor = invis, legend.position = "none"))
      dev.off()
    }
  
  
    
## Plot!
  model.colours <- c("#5F4690","#1D6996","#38A6A5","#0F8554","#73AF48","#EDAD08","#E17C05")
    
  use <- names(tests[[1]])
  plot.data <- list()
  
  for(j in use) { 
    # print(j)
    x <- get.tp.data(j)
    x$Model <- factor(x$Model)
    plot.data[[j]] <- x

    # plot.tp.1(j)
    # plot.tp.2(j)
    plot.tp.3(j)
  }
  
  
## Effective correction plot!!
  ## Think of this as a summary plot
  
  ## First, get the data!
    # apply existing functions
    use <- names(tests[[1]])
    summary.data <- plot.data
    
    # restrict to nTop
    summary.ntop <- lapply(summary.data, function(x) { x[which(x$Metric == "nTop"),]})
    
    # define the baseline performance, i.e. for each fold-change, what's the discriminatory ability when there's no change in composition, using a basic model?
    baseline <- sapply(summary.ntop, function(x) x$Value[which(x$Simulation == "up0.rep1" & x$Model == "lm.uncorrected")])
    
    baseline.095 <- baseline * 0.95
  
  ## Next, for each model, determine the + and - confounds where the new baseline is reached!
    unacceptable <- list()
    models <- levels(as.factor(summary.ntop$confound.p$Model))
    
    for(j in use[-1]) { # [-1] removes the conditions with no confound in expression
      dat <- summary.ntop[[j]]
      # baseline[j] <-
      dat$below.baseline.095 <- dat$Value < baseline.095[j] 
      
      unacceptable[[j]] <- list()
      
      for (k in models) {
        
        negcf.limit <- max(dat$Confound[which(dat$Model == k & dat$below.baseline & dat$Confound < 0)])
        poscf.limit <- min(dat$Confound[which(dat$Model == k & dat$below.baseline & dat$Confound > 0)])
        
        unacceptable[[j]][[k]] <- data.frame(Neg = negcf.limit, Pos = poscf.limit)
      }
      
    }
    
    
    
  ## Now plot!
    dat <- melt(unacceptable)
    colnames(dat) <- c("ConfoundDir", "Unacceptable", "Model", "Condition")
    dat$FC <- (substr(dat$Condition, nchar(dat$Condition) - 2, nchar(dat$Condition)))
    dat$Ct <- "All Ct"; dat$Ct[grep("exc", dat$Condition)] <- "Excitatory-specific"
    
    dat$Model <- factor(dat$Model, levels = levels(as.factor(dat$Model))[c(6,3,4,1,5,2,7)])
    levels(dat$Model) <- names(model.colours)

    dat$Unacceptable <- abs(dat$Unacceptable)
    
    discrimination.plot <- function(dir = "Neg", in.ct = "All Ct") {
      ggplot(dat[which(dat$ConfoundDir == dir  & dat$Ct == in.ct),], aes(y = Unacceptable, x = FC, fill = Model)) +
        theme_bw() +
        theme(panel.border = invis, axis.line = element_line(), legend.text = element_text(size = 7)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.2) +
        stat_summary(aes(colour = Model), fun = function(x) x, geom = "point", position = position_dodge(width = 0.7), size = 1.5) +
        # geom_segment(aes(colour = Model, x = FC, xend = FC, y = 0, yend = Unacceptable), position = position_dodge(width = 0.7)) +
        scale_fill_manual(values = model.colours) +
        scale_colour_manual(values = model.colours) +
        scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0, 0.31), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)) +
        labs(y = "Confound in Composition\nTolerated by Model", x = "Simulated Fold-change") 
      
    }
    
    pdf(file = "Discrimination < 0.95, Negative Confound.pdf", height = 2.5, width = 5)
    discrimination.plot("Neg", "All Ct")
    dev.off()
    
    pdf(file = "Discrimination < 0.95, Positive Confound.pdf", height = 2.5, width = 5)
    discrimination.plot("Pos", "All Ct")
    dev.off()
    
    