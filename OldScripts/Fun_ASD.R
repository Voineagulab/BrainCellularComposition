##################################################################################################################################
## Libraries

# Analysis
  require(limma)
  require(nlme)
  require(DESeq2)
  require(goseq)
  require(WGCNA)
  require(multtest)
  require(doParallel)
  require(variancePartition)
  
# Plotting
  require(gridExtra)
  require(data.table)
  require(cowplot)
  require(reshape2)
  require(ggplot2)
  require(gplots)
  require(VennDiagram)
  require(viridis)
  require(pals)
  require(stringr)
  

##################################################################################################################################
## Analysis

## Algorithm/signature combination to HIGHLIGHT in all analyses; results from other algorithms will also be shown, albeit less
## prominently, and seen through a lens of replication
main.alg <- "CIB.IP"


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

## Plotting theme
theme_gjs <- theme_bw() + theme(panel.border = element_blank(), axis.line = element_line())

## Plotting colours
asd.colours <- c("darkorange", "dodgerblue")
names(asd.colours) <- c("ASD", "CTL")

comp.colours <- alphabet()[c(13, 14)]
names(comp.colours) <- c("CI", "CD")

## Function for fraction of overlaps between two lists
pctOverlap <- function(interestingList, inList) { 
  round(length(which(interestingList %in% inList)) / length(interestingList), digits = 4) 
} 
    
## Function for jaccard similarity
jaccard <- function(list1, list2) { 
  x <- length(intersect(list1, list2))
  y <- length(union(list1, list2))
  z <- x/y
  round(z, digits = 2)
}

  # ## Function to apply sLED to genes in a module
  # run.sLED <- function(condition, module) {
  #   data <- residuals[[condition]]
  #   genes <- im[[condition]][[module]]
  # 
  #   data <- data[which(rownames(data) %in% genes),]
  #   data <- t(data)
  # 
  #   res <- sLED(X = data[which(meta$ASD.CTL == "ASD"),], Y = data[which(meta$ASD.CTL == "CTL"),],
  #                  adj.beta = softPowers[[condition]], rho = 1000, sumabs.seq = 0.1, npermute = 1000,
  #                  useMC = FALSE, mc.cores = 1, seeds = NULL, verbose = TRUE, niter = 20, trace = FALSE)
  # 
  #   leverage <- as.numeric(res$leverage)
  #   names(leverage) <- colnames(data)
  #   leverage <- leverage[order(leverage, decreasing = TRUE)]
  #   leverage <- leverage[which(leverage >= 0.01)] # collect genes accounting for >= 1% of leverage
  #   leverage <- round(leverage, 3)
  # 
  #   output <- list()
  #   output$leverage <- leverage
  #   output$pValue <- round(res$pVal, 3)
  # 
  #   return(output)
  # }
  
  ## Function to plot GOSeq output
    plot.GOseq <- function(x, terms = FALSE, ontology = "BP", topN = 5, legend.position = "none", text.size, y.wrap = FALSE) {
      # x should be a NAMED list, where each element is a different run of GOseq
      stopifnot(class(x) == "list")
      
      # convert x to a filtered dataframe
      for (j in names(x)) x[[j]]$Run <- j # annotate each level
      x <- lapply(x, function(y) {
        y <- y[which(y$ontology %in% ontology),] # filter to the desired ontologies (any of BP, CC, MF)
        y <- y[1:topN,] # filter to the desired number of top hits
        y <- y[,c("over_represented_pvalue", "term", "ontology", "Run")] # filter to relevant columns
        colnames(y) <- c("FDR", "Term", "Ontology", "Run")
        y$FDR <- -log10(y$FDR)
        return(y)
      })
      
      # merge
      x <- do.call("rbind", x) # merge into a single dataframe
      
      # reorder
      z <- c()
      for(j in names(table(x$Run))) {
        k <- x[which(x$Run == j),]
        l <- k[order(k$FDR),]
        z <- rbind(z, l)
      }
      
      x <- z
      if (y.wrap != FALSE) x$Term <- str_wrap(x$Term, width = y.wrap)
      label <- x$Term
      x$Term <- as.factor(paste0(x$Term, x$Run))
      x$Term <- factor(x$Term, levels = x$Term[1:nrow(x)])
      
      # plot
      ggplot(x, aes(x = Term, y = FDR, fill = Run)) +
        geom_col(col = "black", position = position_dodge()) +
        coord_flip() + 
        theme_bw() +
        theme(axis.title.y = element_blank(), panel.border = element_blank(), 
              axis.line.x = element_line(), panel.grid.major.y = element_blank(),
              axis.text.y = element_text(size = text.size), legend.title = element_blank(), legend.position = legend.position) +
        geom_hline(yintercept = c(0, -log10(0.05)), col = c("black", "black"), linetype = c(1,2)) +
        labs(y = "-log10(FDR)") +
        scale_x_discrete(labels = label)
    }
