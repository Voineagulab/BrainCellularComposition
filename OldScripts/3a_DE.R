## In this (revised) script, I analyse those simulate mixtures for testing differential expression

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

## Packages
library(DESeq2)

## Load simulations
load("Sims (Final).rda")


################################################################################################################################ #
## Characterisation of these simulated composition(s) ----

## Group definition
  alpha <- 1:50
  beta <- 51:100

## Boxplot true proportions
  pdf(file = "Meta/True Proportion Boxplots (Final).pdf", height = 2, width = 3)
  for(j in names(sims)) {
    print(j)
    dat <- cbind(sims[[j]]$meta.prop, rep(c("alpha", "beta"), each = 50))
    colnames(dat)[6] <- "Group"
    dat <- melt(dat, id.vars = "Group")

    print(ggplot(dat, aes(x = variable, y = value, colour = Group, fill = Group)) +
            geom_violin(position = position_dodge(width = 1), scale = "width") +
            geom_boxplot(width = 0.5, position = position_dodge(width = 1), outlier.shape =  NA, colour = "black") +
            theme_bw() +
            theme(axis.title.x = invis, panel.border = invis, axis.line = element_line(), legend.position = c(0.9, 0.8),
                  panel.background = element_rect(fill = "grey90")) +
            # scale_x_discrete(labels = c("Neu", "Ast", "Oli", "Mic", "End")) +
            labs(y = "Simulated Proportion", title = j))
  }
  dev.off()


## Estimate composition using DeconRNASeq
  DRS <- list()
  DRS$confound.p <- lapply(sims, function(x) run.DRS(x$confound.p, HCA.sig))
  DRS$confound.pe.1.1 <- lapply(sims, function(x) run.DRS(x$confound.pe.1.1, HCA.sig))
  DRS$confound.pe.1.3 <- lapply(sims, function(x) run.DRS(x$confound.pe.1.3, HCA.sig))
  DRS$confound.pe.1.5 <- lapply(sims, function(x) run.DRS(x$confound.pe.1.5, HCA.sig))
  DRS$confound.pe.2.0 <- lapply(sims, function(x) run.DRS(x$confound.pe.2.0, HCA.sig))
  
  DRS$confound.pe.exc.1.1 <- lapply(sims, function(x) run.DRS(x$confound.pe.exc.1.1, HCA.sig))
  DRS$confound.pe.exc.1.3 <- lapply(sims, function(x) run.DRS(x$confound.pe.exc.1.3, HCA.sig))
  DRS$confound.pe.exc.1.5 <- lapply(sims, function(x) run.DRS(x$confound.pe.exc.1.5, HCA.sig))
  DRS$confound.pe.exc.2.0 <- lapply(sims, function(x) run.DRS(x$confound.pe.exc.2.0, HCA.sig))

  # save and load
  save(DRS, file = "Meta/DeconRNASeq (Final).rda") 


## Extract the true difference in cell-type proportions: cf, standing for "confound"
  # percent confound
  cf <- lapply(sims, function(x) {
    # setup a vector to hold results
    j <- vector("numeric", length = 5)
    names(j) <- colnames(x$meta.prop)
    
    # calculate difference in mean proportion between alpha and beta, for each cell-type
    for (k in names(j)) {
      j[k] <- mean(x$meta.prop[alpha,k]) - mean(x$meta.prop[beta,k])
    }
    return(j)
  })
  
  # p-values of these changes
  cf.p <- lapply(sims, function(x) {
    # setup a vector to hold results
    j <- vector("numeric", length = 5)
    names(j) <- colnames(x$meta.prop)

    # calculate difference in mean proportion between alpha and beta, for each cell-type
    for (k in names(j)) {
      j[k] <- wilcox.test(x$meta.prop[alpha,k], x$meta.prop[beta,k], alternative = "t")$p.value
    }
    return(j)
  })
  
  # save
  save(cf, cf.p, file = "Meta/Percent Confounds (Final).rda") # load("Percent Confounds.rda")
  
  
## Extract the list of perturbed genes
  pert.genes <- sims$dn195.rep1$meta.exp$confound.pe.1.1$genes
  all.genes <- rownames(sims$dn195.rep1$confound.pe.1.1)
  
  all.pert <- all.genes %in% pert.genes$EnsID
  rand.up <- all.genes %in% pert.genes$EnsID[which(pert.genes$forCelltype == "Random" & pert.genes$Direction == "Up")]
  rand.dn <- all.genes %in% pert.genes$EnsID[which(pert.genes$forCelltype == "Random" & pert.genes$Direction == "Down")]
  exc.up <- all.genes %in% pert.genes$EnsID[which(pert.genes$forCelltype == "Excitatory" & pert.genes$Direction == "Up")]
  exc.dn <- all.genes %in% pert.genes$EnsID[which(pert.genes$forCelltype == "Excitatory" & pert.genes$Direction == "Down")]
  
  save(pert.genes, all.genes, all.pert, rand.up, rand.dn, exc.up, exc.dn, file = "Meta/Perturbed Genes (Final).rda")


################################################################################################################################ #
## Run differential expression ----

  
## Setup
  tests <- list()
  runs <- names(DRS)
  
  glm.formulae <- c("~ group", "~ true + group", "~ true + quad + group")
  names(glm.formulae) <- c("glm.uncorrected", "glm.corrected", "glm.q")
  
  
## Loop
  sims.to.run <- names(sims)[grep("rep2", names(sims))]
  
  for (i in sims.to.run) { # for each simulation
    
    tests[[i]] <- list() # new level for the simulation

    for(j in runs) { # test a range of fold-changes

      ## Setup
        print(paste0("Run ", grep(i, names(sims)), " of ", length(sims), ": Models in ", i, "-", j, " starting ", Sys.time()))
        mods <- list()

        # create alternative versions of expression data
        counts <- sims[[i]][[j]]

        logcpm <- apply(counts, 2, function(x) {
          lib.size <- 10^6 / sum(x)
          x <- x * lib.size
          return(x)
        })
        logcpm <- t(as.data.frame(log2(logcpm + 0.5)))

        # create covariates
        covar <- data.frame(true = sims[[i]]$meta.prop$Excitatory,
                             quad = sims[[i]]$meta.prop$Excitatory ^ 2,
                             est = DRS[[j]][[i]]$Excitatory,
                             group = as.factor(rep(c(1, 0), each = 50)))

        # create spline knots
        knots <- quantile(covar$true, p = c(0.25, 0.5, 0.75))


    # ## Run cpm-level models
      mods$lm.uncorrected <- lm(logcpm ~ group, data = covar) # lm without correction
      mods$lm.corrected <- lm(logcpm ~ group + true, data = covar)
      mods$lm.corrected.q <- lm(logcpm ~ group + true + quad, data = covar)
      mods$lm.corrected.e <- lm(logcpm ~ group + est, data = covar)
      mods$spline <- lm(logcpm ~ group + bs(covar$true, knots = knots, degree = 3), data = covar)

      # extract pvalues and effect size
      tests[[i]][[j]] <- lapply(mods, function(y) {
        sum <- summary(y)
        pvals <- sapply(sum, function(y)  { y$coefficients["group1", "Pr(>|t|)"] } )
        log2fc <- sapply(sum, function(y)  { y$coefficients["group1", "Estimate"] } )
        # t <- sapply(sum, function(y)  { y$coefficients["group1", "t value"] } )
        z <- data.frame(pvals = pvals, log2fc = log2fc)
        rownames(z) <- sapply(strsplit(rownames(z), " "), "[", 2)
        return(z)
      })
      

    ## Run GLMs on counts
      for (k in names(glm.formulae)) {
        dds <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = covar,
                                      design = as.formula(glm.formulae[k]))
        dds <- DESeq(dds)
        res <- results(dds, contrast = c("group", "1", "0"), alpha = 0.05)

        res <- as.data.frame(res@listData)
        rownames(res) <- rownames(counts)

        tests[[i]][[j]][[k]] <- data.frame(pvals = res$pvalue, log2fc = res$log2FoldChange)
      }
    }
  }

  ## Save
    save(tests, file = "Models (Final).rda") # load("Linear Models.rda")

  
## Run cell-type specific DE
  ct.specific <- list()  
  
  for (i in sims.to.run) { # for each simulation
    
    ct.specific[[i]] <- list() # new level for the simulation

    for(j in runs) { # test a range of fold-changes
      
      ## Setup
        print(paste0("Run ", grep(i, names(sims)), " of ", length(sims), ": Models in ", i, "-", j, " starting ", Sys.time()))
        ct.specific[[i]][[j]] <- list()

        # create alternative versions of expression data
        counts <- sims[[i]][[j]]
        
        logcpm <- apply(counts, 2, function(x) {
          lib.size <- 10^6 / sum(x)
          x <- x * lib.size
          return(x)
        })
        logcpm <- t(as.data.frame(log2(logcpm + 0.5)))
        
        # create covariates
        covar <- data.frame(true = sims[[i]]$meta.prop$Excitatory,
                            nonexc = 1 - sims[[i]]$meta.prop$Excitatory,
                            group = as.factor(rep(c(1, 0), each = 50)))
        
      ## Run cpm-level models
        mods <- lm(logcpm ~ group + nonexc + nonexc:group, data = covar)

        sum <- summary(mods)
        pvals <- sapply(sum, function(y)  { y$coefficients["group1", "Pr(>|t|)"] } )
        log2fc <- sapply(sum, function(y)  { y$coefficients["group1", "Estimate"] } )
        z <- data.frame(pvals = pvals, log2fc = log2fc)
        rownames(z) <- sapply(strsplit(rownames(z), " "), "[", 2)

        ct.specific[[i]][[j]]$lm <- z
        
      ## Run GLMs
        dds <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = covar,
                                      design = as.formula("~ nonexc + group + nonexc:group"))
        dds <- DESeq(dds)
        res <- results(dds, contrast = c("group", "1", "0"), alpha = 0.05)

        res <- as.data.frame(res@listData)
        rownames(res) <- rownames(counts)

        ct.specific[[i]][[j]]$glm <- data.frame(pvals = res$pvalue, log2fc = res$log2FoldChange)
    }
  }
  
  save(ct.specific, file = "Cell-type Specific Models (Final) (lm).rda")

  
################################################################################################################################ #
## Visualise false positives  ----


## Collect 
  fp <- lapply(tests, function(x) {
    x <- x$confound.p
    sapply(x, function(y) {
      y$FDR <- p.adjust(y$pvals, method = "BH")
      y <- data.frame(Nominal = length(which(y$pvals < 0.05)), 
                      FDR = length(which(y$FDR < 0.05)))
      return(y)
    })
  })  

## Lineplot
  # collect all fdr-level fps
  plot.data <- sapply(fp, function(x) x["FDR",]) # extract from list
  plot.data <- apply(plot.data, 1, as.numeric) # coerce to a numeric matrix
  plot.data <- as.data.frame(plot.data) # coerce to dataframe
  plot.data$Confound <- sapply(cf[grep("rep1", names(cf))], function(x) x["Excitatory"]) 
  
  # convert to plotting format
  plot.data <- melt(plot.data, id.vars = "Confound")
  colnames(plot.data) <- c("Confound", "Model", "FP")
  plot.data <- plot.data[-which(plot.data$Model == "lm.corrected.e"),]
  plot.data$Model <- factor(plot.data$Model)
  levels(plot.data$Model) <- c("Uncorrected LM", "Corrected LM", "Corrected LM Quadratic", "Corrected Spline", "Uncorrected DESeq2", "Corrected DESeq2", "Corrected DESeq2 Quadratic")
  plot.data$Model <- factor(plot.data$Model, levels = levels(plot.data$Model)[c(1, 5, 2, 6, 3, 7, 4)])
  # model.colours <- c("#5F4690","#1D6996","#38A6A5","#0F8554","#73AF48","#EDAD08","#E17C05", "#CC503E") # original, pre-removal of est
  model.colours <- c("#EDAD08", "#CC503E","#1D6996","#38A6A5","#0F8554","#73AF48","#5F4690")
  names(model.colours) <- levels(plot.data$Model)
  
  # shuffle the row order. this reduces the overplotting issue
  plot.data <- plot.data[sample(rownames(plot.data),nrow(plot.data)),]
  

  # line plot, all models
  pdf(file = "FP - Count Across Confound Percent, Line Plot (All Models).pdf", height = 2.5, width = 8)
  ggplot(plot.data, aes(x = Confound, y = FP, colour = Model, fill = Model)) +
    geom_point(size = 1) +
    geom_smooth(se = FALSE, span = 0.3, alpha = 0.3) +
    theme_bw() +
    scale_x_continuous(labels = scales::percent, limits = c(-0.41, 0.4), expand = c(0, 0), breaks = c(-0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 1, 10, 100, 1000, 10000), limits = c(NA, NA)) +
    scale_colour_manual(values = model.colours) +
    labs(y = "Number of False Positives", x = "Difference in Excitatory Neuronal Proportion Between Groups") +
    theme(panel.border = invis, axis.line.y = element_line(), axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right", legend.background = element_rect(fill = alpha("white", 0)),
          panel.grid.minor = invis) 
  dev.off()
  

## Enrichment in markers
  # dataframe identifying markers
    marks <- nTopFeatures(signature = HCA.sig, n = 100, alg = "diff")
  
  # collect all overlaps
  res <- lapply(tests, function(x) {
    # collect up- and down-regulated genes
    x <- data.frame(P = x$confound.p$lm.uncorrected$pvals, Dir = sign(x$confound.p$lm.uncorrected$log2fc), EnsID = rownames(x$confound.p$lm.uncorrected))
    # x <- x[unperturbed,]
    x$Sig <- p.adjust(x$P, method = "BH") < 0.05
    sig <- list(up = x$EnsID[which(x$Dir == 1 & x$Sig)], down = x$EnsID[which(x$Dir == -1 & x$Sig)])
    
    # test enrichment in each of these lists
    y <- lapply(sig, function(z) {
      
      e <- list()
      for (j in levels(as.factor(marks$forCelltype))) {
        e[[j]] <- overEnrich(z, marks$EnsID[which(marks$forCelltype == j)], x$EnsID)
      }
      e <- do.call("rbind", e)
    })
    
    # extract information
    y <- lapply(y, function(z) {
      z1 <- as.numeric(z[,5]) # the column for p-value, one-sided for overEnrichment
      names(z1) <- rownames(z)
      return(z1)
    })
    
    # output
    as.data.frame(do.call("cbind", y))
  })
  
  # augment this list with the confound and cell-type
  m <- match(rownames(res$dn195.rep1), names(cf$dn195.rep1))
  for (j in names(res)) {
    res[[j]]$confound <- cf[[j]][m] # magnitude of confound
    res[[j]]$confound.sig <- cf.p[[j]][m] < 0.05 # significance of confound
    res[[j]]$confound.nDir <- sign(cf[[j]]["Excitatory"]) # direction of confound (neuronal proportion), but the same value is used for all cell-types
    res[[j]]$ct <- rownames(res[[j]]) # cell-type being confounded
    res[[j]]$group <- j # allows grouping of all enrichments that stemmed from a single set of FPs
    res[[j]]$nFP <- as.numeric(fp[[j]][2,1]) # uncorrected model, FDR < 0.05
    
  } 
  
  # convert overlaps into plottable dataframe
  res <- do.call("rbind", res)
  
  # notes:
  # the change in composition is given as mean(alpha) - mean(beta), i.e. positive means increased neuronal composition
  # therefore, when applying sign(), "1" returns when neuronal compostion is increased, and -1 when it's decreased
  
  # notes 2:
  # the sign of differential expression is calculated by sign() on the linear model's coefficient
  # based on a manual check:
  # w <- which.max(tests$ten.sim1$Neurons$uncorrected$coefs)
  # 
  # tests$ten.sim1$Neurons$uncorrected$coefs[w] # +0.78
  # mean(log2(as.numeric(sims$ten.sim1$Neurons[w,alpha]))) # 8.51 is its expression in group alpha
  # mean(log2(as.numeric(sims$ten.sim1$Neurons[w,beta]))) # 5.63 is its expression in group beta
  # it seems that a positive sign() means increased expression in alpha
  
  # colours
  ct.colours <- c("#009392", "#eeb479", "#e88471", "#e9e29c", "#fdfbe4")
  names(ct.colours) <- substr(colnames(HCA.sig), 1, 3)
  
  # plot
  plot.data <- res[which(res$nFP > 100),] # collect cases with > 100 FPs, to give the test ample power
  plot.data$ct <- substr(plot.data$ct, 1, 3)
  
  # plot confound.nDir == -1, i.e. samples with DECREASED neuronal composition
  pA <- melt(plot.data[which(plot.data$confound.nDir == -1),c(1,2,6)])
  pA$variable <- factor(pA$variable, levels = c("down", "up")) # converting into a factor, and reordering the factor
  levels(pA$variable) <- c("Downregulated FPs", "Upregulated FPs") # Renaming the levels of the factor, for better labelling on the plot
  
  pA <- ggplot(pA, aes(x = variable, y = -log10(value), fill = ct)) +
    geom_violin(colour = "black", scale = "width", position = position_dodge(width = 0.7), draw_quantiles = 0.5) +
    # geom_boxplot(fill = "white", width = 0.2, position = position_dodge(width = 0.7)) +
    geom_point(colour = "black", position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.7), size = 0.5, alpha = 0.1) +
    theme_bw() +
    scale_fill_manual(values = ct.colours) +
    scale_colour_manual(values = ct.colours) +
    theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis, axis.title.x = invis) +
    labs(y = "Enrichment for cell-type markers\nwithin FPs (-log10(P))", title = "Decreased Excitatory %") 
  
  # plot confound.nDir == 1, i.e. samples with INCREASED neuronal composition
  pB <- melt(plot.data[which(plot.data$confound.nDir == 1),c(1,2,6)])
  pB$variable <- factor(pB$variable, levels = c("down", "up")) # converting into a factor, and reordering the factor
  levels(pB$variable) <- c("Downregulated FPs", "Upregulated FPs") # Renaming the levels of the factor, for better labelling on the plot
  pB <- ggplot(pB, aes(x = variable, y = -log10(value), fill = ct)) +
    geom_violin(colour = "black", scale = "width", position = position_dodge(width = 0.7), draw_quantiles = 0.5) +
    geom_point(colour = "black", position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.7), size = 0.5, alpha = 0.1) +
    theme_bw() +
    scale_fill_manual(values = ct.colours) +
    theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis, axis.title.y = invis, axis.title.x = invis) +
    guides(fill = guide_legend("Celltype")) +
    labs(y = "Enrichment for cell-type markers\nwithin FPs (-log10(P))", title = "Increased Excitatory %") 
  
  pdf(file = "FP - Marker Enrichment Within FPs (LM Uncorrected).pdf", height = 2.3, width = 8)
  plot_grid(pA + theme(legend.position = "none"), pB, nrow = 1, rel_widths = c(1, 1.1))
  dev.off()   


################################################################################################################################ #
## Visualise false negatives ----

## First, get the vectors of perturbed genes
  load("Meta/Perturbed Genes (Final).rda")
  
## Plotting functions
  ## Function to extract plotting data
    get.tp.data <- function(fc) {
      
      plot.data <- lapply(tests, function(x) {
        # get genes which have an adjusted p-value of < 0.05
        sig <- sapply(x[[fc]], function(y) { p.adjust(y$pvals, method = "BH")  * sign(y$log2fc)  })
        rownames(sig) <- rownames(x$confound.pe.exc.1.3$lm.uncorrected)
        
        
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
    
  ## Function for plotting style 1
    # model.colours <- c("#5F4690","#1D6996","#38A6A5","#0F8554","#73AF48","#EDAD08","#E17C05")
    plot.tp.1 <- function(fc) {
      # setup
      filename <- paste0("TP Plot, Style 1 (", fc, ").pdf")
      print(paste0("Saving plot at ", filename))
      x <- plot.data[[fc]]
      x <- x[sample(rownames(x), nrow(x)),] # this makes the overplotting of points less of an issue
      # levels(x$Model) <- c("LM", "LM Correction", "Quadratic Correction", "LM Correction (Est Prop)", "Spline", "GLM", "GLM Corrected")
      levels(x$Model) <- names(model.colours)
      
      # create plots
      # pA <- ggplot(x[which(x$Metric == "nFP"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
      #   geom_point(size = 1) +
      #   theme_bw() +
      #   geom_smooth(se = FALSE, span = 0.3) +
      #   scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
      #   scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 1, 10, 100, 1000, 10000)) +
      #   scale_colour_manual(values = model.colours) +
      #   labs(y = "Number of False Positives", x = "") +
      #   theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
      #         legend.background = element_rect(fill = alpha("white", 0)),
      #         panel.grid.minor = invis, legend.position = "none", axis.title.x = invis)
      
      pB <- ggplot(x[which(x$Metric == "nTP.exc.dn"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
        geom_point(size = 1) +
        geom_smooth(se = FALSE, span = 0.3) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_colour_manual(values = model.colours) +
        labs(y = "TPR for Downregulated\nExcitatory Markers", x = "Difference in Excitatory Proportion") +
        theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
              legend.background = element_rect(fill = alpha("white", 0)),
              panel.grid.minor = invis, legend.position = "none", axis.title.x = invis)
      
      pC <- ggplot(x[which(x$Metric == "nTP.exc.up"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
        geom_point(size = 1) +
        geom_smooth(se = FALSE, span = 0.3) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_colour_manual(values = model.colours) +
        labs(y = "TPR for Upregulated\nExcitatory Markers", x = "Difference in Excitatory Proportion") +
        theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
              legend.background = element_rect(fill = alpha("white", 0)),
              panel.grid.minor = invis, legend.position = "none", axis.title.x = invis)
      
      pD <- ggplot(x[which(x$Metric == "nTP.rand.dn"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
        geom_point(size = 1) +
        geom_smooth(se = FALSE, span = 0.3) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_colour_manual(values = model.colours) +
        labs(y = "TPR for Downregulated\nNon-markers", x = "") +
        theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
              legend.background = element_rect(fill = alpha("white", 0)),
              panel.grid.minor = invis, legend.position = "none")
      
      pE <- ggplot(x[which(x$Metric == "nTP.rand.up"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
        geom_point(size = 1) +
        geom_smooth(se = FALSE, span = 0.3) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent) +
        scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
        scale_colour_manual(values = model.colours) +
        labs(y = "TPR for Upregulated\nNon-markers", x = "Difference in Excitatory Proportion") +
        theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
              legend.background = element_rect(fill = alpha("white", 0)),
              panel.grid.minor = invis, legend.position = "none")
      
      
      pF <- ggplot(x[which(x$Metric == "nTop"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
        geom_point(size = 1) +
        geom_smooth(se = FALSE, span = 0.3) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_colour_manual(values = model.colours) +
        guides(colour = guide_legend(ncol = 2)) +
        labs(y = "Fraction of Perturbed Genes\n in 200 Most Significant", x = "Difference in Excitatory Proportion") +
        theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
              legend.background = element_rect(fill = alpha("white", 0)),
              panel.grid.minor = invis, legend.position = "right")
      
      
      leg <- get_legend(pF)
      pF <- pF + theme(legend.position = "none")
      
      
      ## Combined plot
      
      pdf(file = filename, height = 5.5, width = 8)
      print(plot_grid(pB, pC, pD, pE, pF, leg, nrow = 3))
      # print(plot_grid(leg, p, ncol = 1, rel_heights = c(1,10)))
      dev.off()
      
    }
    
  ## Below is a reworked version of plot.tp.1, with a different focus...
    plot.tp.2 <- function(fc) {
      # setup
      filename <- paste0("TP Plot, Style 2 (", fc, ").pdf")
      print(paste0("Saving plot at ", filename))
      x <- plot.data[[fc]]
      x <- x[sample(rownames(x), nrow(x)),] # this makes the overplotting of points less of an issue
      # levels(x$Model) <- c("LM", "LM Correction", "Quadratic Correction", "LM Correction (Est Prop)", "Spline", "GLM", "GLM Corrected")
      levels(x$Model) <- names(model.colours)

      # x <- x[which(x$Metric ==),]

      # create plots
      pC <- ggplot(x[which(x$Metric == "nTP.exc.up"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
        geom_point(size = 0.5, alpha = 0.5) +
        geom_smooth(se = FALSE, span = 0.3) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_colour_manual(values = model.colours) +
        labs(y = "TPR for Upregulated\nExcitatory Markers", x = "Difference in Excitatory Proportion") +
        theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
              legend.background = element_rect(fill = alpha("white", 0)),
              panel.grid.minor = invis, legend.position = "none", axis.title.x = invis)

      pE <- ggplot(x[which(x$Metric == "nTP.rand.up"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
        geom_point(size = 0.5, alpha = 0.5) +
        geom_smooth(se = FALSE, span = 0.3) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent) +
        scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
        scale_colour_manual(values = model.colours) +
        labs(y = "TPR for Upregulated\nNon-markers", x = "Difference in Excitatory Proportion") +
        theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
              legend.background = element_rect(fill = alpha("white", 0)),
              panel.grid.minor = invis, legend.position = "none")

      line <- 0.95 * x$Value[which(x$Simulation == "up0.rep1" & x$Model == "LM" & x$Metric == "nTop")]

      pF <- ggplot(x[which(x$Metric == "nTop"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
        geom_point(size = 0.5, alpha = 0.5) +
        geom_smooth(se = FALSE, span = 0.3) +
        geom_hline(yintercept =  line, linetype = 2) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_colour_manual(values = model.colours) +
        labs(y = "Fraction of Perturbed Genes\n in 200 Most Significant", x = "") +
        theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
              legend.background = element_rect(fill = alpha("white", 0)),
              panel.grid.minor = invis, legend.position = "top", axis.title.x = invis)


      leg <- get_legend(pF)
      pF <- pF + theme(legend.position = "none")

      ## Combined plot
      p <- (plot_grid(pF, pC, pE + theme(axis.title.x = invis), nrow = 1))
      pdf(file = filename, height = 3, width = 8)
      print(plot_grid(leg, p, ncol = 1, rel_heights = c(1.5,10)))
      dev.off()

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
  use <- names(tests$dn195.rep1)
  plot.data <- list()
  
  for(j in use) { 
    # print(j)
    x <- get.tp.data(j)
    x <- x[-which(x$Model == "lm.corrected.e"),]
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
    use <- names(tests$dn195.rep1)
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
    
    pdf(file = "Discrimination < 0.95, Negative Confound, Excitatory-specific.pdf", height = 2.5, width = 5)
    discrimination.plot("Neg", "Excitatory-specific")
    dev.off()
    
    pdf(file = "Discrimination < 0.95, Positive Confound, Excitatory-specific.pdf", height = 2.5, width = 5)
    discrimination.plot("Pos", "Excitatory-specific")
    dev.off()
    
    
    
################################################################################################################################ #
## Cell-type-specific expression! ----

## First, get the vectors of perturbed genes
load("Meta/Perturbed Genes (Final).rda")    
  
## Number of false positives!!
    fp <- lapply(ct.specific, function(x) {
      x <- x$confound.p
      sapply(x, function(y) {
        y$FDR <- p.adjust(y$pvals, method = "BH")
        y <- data.frame(Nominal = length(which(y$pvals < 0.05)), 
                        FDR = length(which(y$FDR < 0.05)))
        return(y)
      })
    })  
    
    ## Lineplot
    # collect all fdr-level fps
    # plot.data <- sapply(fp, function(x) x["FDR",]) # extract from list
    plot.data <- sapply(fp, function(x) cbind(x["FDR",], x["FDR",])) # extract from list
    plot.data <- apply(plot.data, 1, as.numeric) # coerce to a numeric matrix
    plot.data <- as.data.frame(plot.data) # coerce to dataframe
    plot.data$Confound <- sapply(cf[grep("rep1", names(cf))], function(x) x["Excitatory"]) 
    
    # convert to plotting format
    plot.data <- melt(plot.data, id.vars = "Confound")
    colnames(plot.data) <- c("Confound", "Model", "FP")
    
    # shuffle the row order. this reduces the overplotting issue
    # plot.data <- plot.data[sample(rownames(plot.data),nrow(plot.data)),]
    
    
    # reannotate levels
    
    #CC503E,#94346E,#6F4070,#994E95,#666666
    
    # model.colours <- c("#1D6996", "#38A6A5","#E17C05", "#EDAD08", "#94346E")
    # 
    
    # levels(plot.data$Model) <- c("Uncorrected", "LM-corrected (True Prop)", "LM-corrected (Est Prop)", "Spline-corrected (True Prop)")
    
    #5F4690,,#,#0F8554,#73AF48,#EDAD08,#E17C05,#CC503E,#94346E,#6F4070,#994E95,#666666
    
    # line plot, all models
    pdf(file = "Ct-specific - False Positives Across Confound Percent, Line Plot (lm Only).pdf", height = 2.5, width = 8)
    ggplot(plot.data, aes(x = Confound, y = FP)) +
      geom_point(size = 1) +
      # geom_smooth(se = FALSE, span = 0.3, alpha = 0.3) +
      theme_bw() +
      scale_x_continuous(labels = scales::percent, limits = c(-0.41, 0.4), expand = c(0, 0), breaks = c(-0.4, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.4)) +
      scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 1, 10, 100, 1000, 10000), limits = c(NA, NA)) +
      # scale_colour_manual(values = model.colours) +
      labs(y = "Number of False Positives", x = "Difference in Excitatory Proportion") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
            legend.position = "right", legend.background = element_rect(fill = alpha("white", 0)),
            panel.grid.minor = invis) 
    dev.off()
    
## True positives rates!!
  ## Write function
  analyse.ct.specific <- function(fc) {
    
    plot.data <- lapply(ct.specific, function(x) {
      # get genes which have an adjusted p-value of < 0.05
      sig <- sapply(x[[fc]], function(y) { p.adjust(y$pvals, method = "BH")  * sign(y$log2fc)  })
      rownames(sig) <- rownames(x$confound.pe.exc.1.3$lm.uncorrected)
      
      
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
  
  plot.cts.1 <- function(fc) {
    # setup
    filename <- paste0("Ct-specific TP Plot, Style 1 (", fc, ").pdf")
    print(paste0("Saving plot at ", filename))
    x <- plot.data[[fc]]
    x <- x[sample(rownames(x), nrow(x)),] # this makes the overplotting of points less of an issue
    # levels(x$Model) <- c("LM", "LM Correction", "Quadratic Correction", "LM Correction (Est Prop)", "Spline", "GLM", "GLM Corrected")
    # levels(x$Model) <- names(model.colours)
    
    # create plots
    # pA <- ggplot(x[which(x$Metric == "nFP"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
    #   geom_point(size = 1) +
    #   theme_bw() +
    #   geom_smooth(se = FALSE, span = 0.3) +
    #   scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
    #   scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 1, 10, 100, 1000, 10000)) +
    #   scale_colour_manual(values = model.colours) +
    #   labs(y = "Number of False Positives", x = "") +
    #   theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
    #         legend.background = element_rect(fill = alpha("white", 0)),
    #         panel.grid.minor = invis, legend.position = "none", axis.title.x = invis)
    
    pB <- ggplot(x[which(x$Metric == "nTP.exc.dn"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
      geom_point(size = 1) +
      geom_smooth(se = FALSE, span = 0.3) +
      theme_bw() +
      scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1)) +
      # scale_colour_manual(values = model.colours) +
      labs(y = "TPR for Downregulated\nExcitatory Markers", x = "Difference in Excitatory Proportion") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
            legend.background = element_rect(fill = alpha("white", 0)),
            panel.grid.minor = invis, legend.position = "none", axis.title.x = invis)
    
    pC <- ggplot(x[which(x$Metric == "nTP.exc.up"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
      geom_point(size = 1) +
      geom_smooth(se = FALSE, span = 0.3) +
      theme_bw() +
      scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1)) +
      # scale_colour_manual(values = model.colours) +
      labs(y = "TPR for Upregulated\nExcitatory Markers", x = "Difference in Excitatory Proportion") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
            legend.background = element_rect(fill = alpha("white", 0)),
            panel.grid.minor = invis, legend.position = "none", axis.title.x = invis)
    
    pD <- ggplot(x[which(x$Metric == "nTP.rand.dn"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
      geom_point(size = 1) +
      geom_smooth(se = FALSE, span = 0.3) +
      theme_bw() +
      scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1)) +
      # scale_colour_manual(values = model.colours) +
      labs(y = "TPR for Downregulated\nNon-markers", x = "") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
            legend.background = element_rect(fill = alpha("white", 0)),
            panel.grid.minor = invis, legend.position = "none")
    
    pE <- ggplot(x[which(x$Metric == "nTP.rand.up"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
      geom_point(size = 1) +
      geom_smooth(se = FALSE, span = 0.3) +
      theme_bw() +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
      # scale_colour_manual(values = model.colours) +
      labs(y = "TPR for Upregulated\nNon-markers", x = "Difference in Excitatory Proportion") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
            legend.background = element_rect(fill = alpha("white", 0)),
            panel.grid.minor = invis, legend.position = "none")
    
    
    pF <- ggplot(x[which(x$Metric == "nTop"),], aes(x = Confound, y = Value, colour = Model, fill = Model)) +
      geom_point(size = 1) +
      geom_smooth(se = FALSE, span = 0.3) +
      theme_bw() +
      scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1)) +
      # scale_colour_manual(values = model.colours) +
      guides(colour = guide_legend(ncol = 2)) +
      labs(y = "Fraction of Perturbed Genes\n in 200 Most Significant", x = "Difference in Excitatory Proportion") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.ticks.x = invis,
            legend.background = element_rect(fill = alpha("white", 0)),
            panel.grid.minor = invis, legend.position = "right")
    
    
    leg <- get_legend(pF)
    pF <- pF + theme(legend.position = "none")
    
    
    ## Combined plot
    
    pdf(file = filename, height = 5.5, width = 8)
    print(plot_grid(pB, pC, pD, pE, pF, leg, nrow = 3))
    # print(plot_grid(leg, p, ncol = 1, rel_heights = c(1,10)))
    dev.off()
    
  }
  
    
  ## Apply function
    use <- names(ct.specific$dn195.rep1)[-1]
    plot.data <- list()
    for(j in use) { 
      print(j)
      # plot.data[[j]] <- analyse.ct.specific(j)
      plot.cts.1(j)
    }
    
    
