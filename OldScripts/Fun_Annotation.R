################################################################################################################################ #
## Generate annotation files ----

setwd("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/")
source("../Scripts/Fun_Preprocessing.R")

## Gene information files
  load(file = "Preprocessed/geneInfo.rda") # its generation is documented below
  # gtf <- read.table("/Volumes/MacintoshHD_RNA/Users/rna/REFERENCE/HUMAN/Ensembl_GRCh37_hg19/ANNOTATION/GTF/genes.gtf", sep = "\t")
  # geneInfo <- strsplit(as.character(gtf[,9]), split = ";") # split column 9, which contains useful annotation information
  # geneInfo <- sapply(geneInfo, `[`, c(2,3,4,5))
  # geneInfo <- as.data.frame(t(geneInfo))
  # geneInfo <- unique(geneInfo) # from 2,280,612 to 154,486 entries
  # 
  # # I need information on gene_biotype, gene_id, and gene_name. However, these 3 bits of information fall sequentially in eithers columns 1-3 or 2-4; this is how I fix it to always be columns 2-4, after which I remove column 1
  # g <- grep("gene_biotype", geneInfo$V1)
  # geneInfo$V4[g] <- geneInfo$V3[g]
  # geneInfo$V3[g] <- geneInfo$V2[g]
  # geneInfo$V2[g] <- geneInfo$V1[g] # g <- grep("gene_biotype", geneInfo$V2); f <- grep("gene_id", geneInfo$V3); h <- grep("gene_name", geneInfo$V4) # sanity checks: should equal nrow of geneInfo
  # geneInfo <- geneInfo[,-1]
  # colnames(geneInfo) <- c("Biotype", "ensID", "Gene.Symbol")
  # geneInfo <- unique(geneInfo)
  # 
  # # Each entry is formatted as, say, "columnName information"; here, I remove the columnName from each entry
  # geneInfo$Biotype <- sapply(strsplit(geneInfo$Biotype, split = " "), `[`, 3)
  # geneInfo$ensID <- sapply(strsplit(geneInfo$ensID, split = " "), `[`, 3)
  # geneInfo$Gene.Symbol <- sapply(strsplit(geneInfo$Gene.Symbol, split = " "), `[`, 3)
  # 
  # save(geneInfo , file = "Preprocessed/geneInfo.rda") # all ensIDs are unique, but only 52775/62069 Gene.Symbols are unique. Using table(geneInfo$Biotype[which(duplicated(geneInfo$Gene.Symbol))]), only 2226 / 9294 duplicated genes are protein coding. However, there are only 22834 protein coding genes...
  
## Gene length file
  load("Preprocessed/exonicLength.rda") # its generation is documented below
  #txdb=makeTxDbFromGFF("/Volumes/MacintoshHD_RNA/Users/rna/ANNOTATION//GENCODE_human/gencode.v19.annotation.gtf", format="gtf")
  #exons.list.per.gene <- exonsBy(txdb,by="gene")
  #exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
  #exonicLength=as.matrix(exonic.gene.sizes )
  #colnames(exonicLength)="exoniclength_bp"
  #rownames(exonicLength) <- transpose(strsplit(rownames(exonicLength), ".", fixed = TRUE))[[1]]
  #save(exonicLength, file = "Preprocessed/exonicLength.rda")   
  
