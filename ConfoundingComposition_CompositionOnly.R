## In this (revised) script, I analyse those simulate mixtures for testing differential expression

################################################################################################################################ #
## Generic setup  ----

## Start!
  rm(list = ls()); gc()

## Directory
  # root.dir <- "INSERT.YOUR.ROOT.DIRECTORY.PLEASE/"
  root.dir <- "/Volumes/Data1/PROJECTS/BrainCellularComposition/BrainCellularComposition_GitHub/"
  setwd(paste0(root.dir, "Results/ConfoundingComposition/CompositionOnly/"))

## Functions
source("../../../Scripts/Fun_Composition.R")
source("../../../Scripts/Fun_Preprocessing.R")

## Packages
library(DESeq2)
library(splines)

## Load simulations
load("../../../Data/Preprocessed/ConfoundingComposition_Simulations.rda")


################################################################################################################################ #
## Characterisation of these simulated composition(s) ----

## Group definition
  alpha <- 1:50
  beta <- 51:100

## Boxplot true proportions
  pdf(file = "True Proportion Boxplots.pdf", height = 2, width = 3)
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
  save(cf, cf.p, file = "Percent Confounds.rda") # load("Percent Confounds.rda")
  
  
## Extract the list of perturbed genes
  pert.genes <- sims[[1]]$meta.exp$confound.pe.1.1$genes
  all.genes <- rownames(sims[[1]]$confound.pe.1.1)
  
  all.pert <- all.genes %in% pert.genes$EnsID
  rand.up <- all.genes %in% pert.genes$EnsID[which(pert.genes$forCelltype == "Random" & pert.genes$Direction == "Up")]
  rand.dn <- all.genes %in% pert.genes$EnsID[which(pert.genes$forCelltype == "Random" & pert.genes$Direction == "Down")]
  exc.up <- all.genes %in% pert.genes$EnsID[which(pert.genes$forCelltype == "Excitatory" & pert.genes$Direction == "Up")]
  exc.dn <- all.genes %in% pert.genes$EnsID[which(pert.genes$forCelltype == "Excitatory" & pert.genes$Direction == "Down")]
  
  save(pert.genes, all.genes, all.pert, rand.up, rand.dn, exc.up, exc.dn, file = "Perturbed Gene List.rda")


################################################################################################################################ #
## Run differential expression ----

 ## Setup
  tests <- list()
  runs <- names(sims[[1]])[-c(2,3,5)]
  
  glm.formulae <- c("~ group", "~ true + group", "~ true + quad + group")
  names(glm.formulae) <- c("glm.uncorrected", "glm.corrected", "glm.q")
  
  
## Loop
  sims.to.run <- names(sims)
  sims.to.run <- sims.to.run[grep("rep1", sims.to.run)] # just doing rep1! Feel free to change this 
  
  # expect ~2minutes per run (j), so ~20min per sims.to.run (i)
  
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
                             group = as.factor(rep(c(1, 0), each = 50)))

        # create spline knots
        knots <- quantile(covar$true, p = c(0.25, 0.5, 0.75))


    # ## Run cpm-level models
      mods$lm.uncorrected <- lm(logcpm ~ group, data = covar) # lm without correction
      mods$lm.corrected <- lm(logcpm ~ group + true, data = covar)
      mods$lm.corrected.q <- lm(logcpm ~ group + true + quad, data = covar)
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
    save(tests, file = "Models.rda") # load("Linear Models.rda")


  
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
  plot.data$Model <- factor(plot.data$Model)
  # levels(plot.data$Model) <- c("Uncorrected LM", "Corrected LM", "Corrected LM Quadratic", "Corrected Spline", "Uncorrected DESeq2", "Corrected DESeq2", "Corrected DESeq2 Quadratic")
  # plot.data$Model <- factor(plot.data$Model, levels = levels(plot.data$Model)[c(1, 5, 2, 6, 3, 7, 4)])
  # model.colours <- c("#EDAD08", "#CC503E","#1D6996","#38A6A5","#0F8554","#73AF48","#5F4690")
  # names(model.colours) <- levels(plot.data$Model)
  
  # shuffle the row order. this reduces the overplotting issue
  plot.data <- plot.data[sample(rownames(plot.data),nrow(plot.data)),]
  

  # line plot, all models
  pdf(file = "FP Count Across Confound Percent, Line Plot (All Models).pdf", height = 2.5, width = 8)
  ggplot(plot.data, aes(x = Confound, y = FP, colour = Model, fill = Model)) +
    geom_point(size = 1) +
    geom_smooth(se = FALSE, span = 0.3, alpha = 0.3) +
    theme_bw() +
    scale_x_continuous(labels = scales::percent, limits = c(-0.41, 0.4), expand = c(0, 0), breaks = c(-0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 1, 10, 100, 1000, 10000), limits = c(NA, NA)) +
    # scale_colour_manual(values = model.colours) +
    labs(y = "Number of False Positives", x = "Difference in Excitatory Neuronal Proportion Between Groups") +
    theme(panel.border = invis, axis.line.y = element_line(), axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right", legend.background = element_rect(fill = alpha("white", 0)),
          panel.grid.minor = invis) 
  dev.off()
  

## Enrichment in markers
  # dataframe identifying markers
    CA.sig <- read.csv("../../../Data/Preprocessed/ConfoundingComposition_Signature.csv", row.names = 1)
    marks <- nTopFeatures(signature = CA.sig, n = 100, alg = "diff")
  
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
  m <- match(rownames(res[[1]]), names(cf[[1]]))
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
  names(ct.colours) <- substr(colnames(CA.sig), 1, 3)
  
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
  
  pdf(file = "Marker Enrichment Within FPs (LM Uncorrected).pdf", height = 2.3, width = 8)
  plot_grid(pA + theme(legend.position = "none"), pB, nrow = 1, rel_widths = c(1, 1.1))
  dev.off()   



 