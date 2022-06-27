# Overview
This directory contains source code for reproducing the analysis in:

Sutton, G.J., Poppe D., Simmons R.K., Walsh K., Nawaz U., Lister R., Gagnon-Bartsch J.A., and Voineagu I. Comprehensive evaluation of deconvolution methods for human brain gene expression. Nat Commun 13, 1358 (2022). https://doi.org/10.1038/s41467-022-28655-4

File naming conventions are as follows, and should be followed in order:
- Group 1: scripts that preprocess public and generated data 
- Group 2: estimating composition in simulated mixtures, and analyses of its output. Broadly corresponds to Figures 1-3. RME: In vitro RNA mixtures of cultured neuronal and astrocyte RNA. SCME: Mixtures derived by simulating pseudobulks from Darmanis et al.'s single cell RNA-seq. SNME: Mixtures derived by simulating pseudobulks from Velmeshev et al.'s or the Human Cell Atlas' single nucleus RNA-seq (VL and CA, respectively).
- Group 3: simulation of datasets with confounded composition between groups, exploring DE driven by this phenomenon. Broadly corresponds to Figures 4-5. CIBx: analyses using CIBERSORTx to impute cell-type-specific expression
- Group 4: analyses of public bulk brain RNA-seq from Parikshak et al., the GTEx consortium, and an autism spectrum disorder (ASD) dataset found in Parikshak et al. Broadly corresponds to Figure 6.
- Group 5: miscellaneous analysis scripts. CrossSignatureCorrelations: for plotting the relationship between signatures. CrossTissue: deconvolution of heart and pancreas data. Immunopanned: deconvolution of Zhang et al.'s immunopanned pure brain cell-types. RNA_Content: exploring whether deconvolution estimates cellular or RNA proportion
- Fun: contains general custom functions for calling in other scripts

# Datasets
The sequencing data generated in this study have been deposited in the GEO database under accession code GSE175772. 

### For generating signatures
Processed signature data can be accessed as [Supplementary Data 5](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28655-4/MediaObjects/41467_2022_28655_MOESM8_ESM.xlsx) in our study. Links to access the raw data underlying this can be found in the table below.

| Name  | Species | Description | Reference and Link | Used In | 
| ----- | ------- | --------- | ---- | ------- |
| CA | Human | Post-mortem brain tissue, smartSeq2 snRNA-seq | Hodge, R. D. et al. Conserved cell types with divergent features in human versus mouse cortex. Nature 573, 61–68 (2019). [Human Cell Atlas website](http://portal.brain-map.org/)| ------- |
| DM | Human | Surgically-resected brain tissue, smartSeq scRNA-seq | Darmanis, S. et al. A survey of human brain transcriptome diversity at the single cell level. Proc. Natl Acad. Sci. USA 112, 7285–7290 (2015). [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835), or processed version on [GitHub](https://github.com/VCCRI/CIDR-examples/tree/master/Brain) | ------- |
| F5 | Human | Cultured brain cells, bulk CAGE-seq | Forrest, A. R. R. et al. A promoter-level mammalian expression atlas. Nature 507, 462–470 (2014). [FANTOM5 Website](http://fantom.gsc.riken.jp/5/data/)  | ----|
| NG | Human | Post-mortem brain tissue, 10X snRNA-seq | Nagy, C. et al. Single-nucleus transcriptomics of the prefrontal cortex in major depressive disorder implicates oligodendrocyte precursor cells and excitatory neurons. Nat. Neurosci. 23, 771–781 (2020). [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144136) | ------- |
| VL | Human | Post-mortem brain tissue, 10X snRNA-seq | Velmeshev, D. et al. Single-cell genomics identifies cell type–specific molecular changes in autism. Science (80-.) 364, 685–689 (2019). [Data browser](https://autism.cells.ucsc.edu/) | ------- |
| LK | Human | Post-mortem brain tissue, 10X snRNA-seq | Lake, B. B. et al. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain. Nat. Biotechnol. 36, 70–80 (2018). [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97942) | ---- |
| IP | Human | Surgically-resected brain tissue, bulk RNA-seq | Zhang, Y. et al. Purification and characterization of progenitor and mature human astrocytes reveals transcriptional and functional differences with mouse. Neuron 89, 37–53 (2016). [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73721) | ------- |
| MM | Mouse | Mouse brain, bulk RNA-seq | Zhang, Y. et al. An RNA-sequencing transcriptome and splicing database of glia, neurons, and vascular cells of the cerebral cortex. J. Neurosci. 34, 11929–11947 (2014). [Paper website](https://web.stanford.edu/group/barres_lab/brain_rnaseq.html)  | ---- |
| TS | Mouse | Mouse brain, SmartSeq2 scRNA-seq | Tasic, B. et al. Shared and distinct transcriptomic cell types across neocortical areas. Nature 563, 72–78 (2018). [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115746) | ------- |

### Bulk RNA-seq datasets deconvolved as exemplars
| Name  | Samples | Details | Reference and Link | 
| ----- | ------- | --------- | ---- |
| Parikshak | Human | Ribo-depleted RNA-seq from 251 post-mortem samples including frontal cortex, temporal cortex, and cerebellar vermis samples from 48 ASD and 49 control individuals, aged 2–67 | Parikshak, N. N. et al. Genome-wide changes in lncRNA, splicing, and regional gene expression patterns in autism. Nature 540, 423–427 (2016). [GitHub](https://github.com/dhglab/Genome-wide-changes-in-lncRNA-alternative-splicing-and-cortical-patterning-in-autism/releases) |
| GTEx V7 Brain | Human | Poly-A+ RNA-seq of 1671 post-mortem samples, including 13 subregions which we classified as cortical, sub-cortical, cerebellar, or spinal | Consortium, Gte. Genetic effects on gene expression across human tissues. Nature 550, 204–213 (2017). [GTEx website](https://gtexportal.org/home/datasets) |



### Non-brain data
We also generated cell-type-specific signatures for two non-brain tissues, to test the tissue-specificity of our results.

| Name  | Tissue | Description | Reference and Link | Used In | 
| ----- | ------- | --------- | ---- | ------- |
| EN | Pancreas | --------- | Enge, M. et al. Single-cell analysis of human pancreas reveals transcriptional signatures of aging and somatic mutation patterns. Cell 171, 321–330e14 (2017). [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81547) | ----  |
| BL | Pancreas | --------- | Blodgett, D. M. et al. Novel observations from next-generation RNA sequencing of highly purified human adult and fetal islet cell subsets. Diabetes 64, 3172–3181 (2015). [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67543) | ------- |
| FS | Pancreas | --------- | Furuyama, K. et al. Diabetes relief in mice by glucose-sensing insulin-secreting human α-cells. Nature 567, 43–48 (2019). [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117454) | ------- |
| FG | Pancreas | --------- | Furuyama, K. et al. Diabetes relief in mice by glucose-sensing insulin-secreting human α-cells. Nature 567, 43–48 (2019). [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117454) | ------- |
| F5 | Heart | --------- | Forrest, A. R. R. et al. A promoter-level mammalian expression atlas. Nature 507, 462–470 (2014). [FANTOM5 Website](http://fantom.gsc.riken.jp/5/data/)  | ------- |
| EN | Heart | --------- | Djebali, S. et al. Landscape of transcription in human cells. Nature 489, 101–108 (2012). [ENCODE](https://www.encodeproject.org/publication-data/ENCSR590RJC/) | ------- |
| SC | Heart | --------- | Wang, L. et al. Single-cell reconstruction of the adult human heart during heart failure and recovery reveals the cellular landscape underlying cardiac function. Nat. Cell Biol. 22, 108–119 (2020). [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109816) | ------- |

# Deconvolution Algorithms

| Method  | Classification | Reference and installation | 
| ------- | -------------- | -------------------------- |
| DeconRNASeq  | Partial deconvolution  | Gong, T. & Szustakowski, J. D. DeconRNASeq: a statistical framework for deconvolution of heterogeneous tissue samples based on mRNA-Seq data. Bioinformatics 29, 1083–1085 (2013). [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/DeconRNASeq.html) |
| CIBERSORT  | Partial deconvolution  | Newman, A. M. et al. Robust enumeration of cell subsets from tissue expression profiles. Nat. Methods 12, 453–457 (2015). Note: the R source for this package is available by request at the [CIBERSORT website](https://cibersortx.stanford.edu/)  |
| dtangle  | Partial deconvolution  | Hunt, G. J., Freytag, S., Bahlo, M. & Gagnon-Bartsch, J. A. dtangle: accurate and robust cell type deconvolution. Bioinformatics 290262 (2018). [CRAN](https://cran.r-project.org/web/packages/dtangle/index.html)|
| MuSiC  | Partial deconvolution  | Wang, X., Park, J., Susztak, K., Zhang, N. R. & Li, M. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nat. Commun. 10, 380 (2019). [GitHub](https://xuranw.github.io/MuSiC/articles/MuSiC.html)|
| Linseed  | Complete deconvolution  | Zaitsev, K., Bambouskova, M., Swain, A. & Artyomov, M. N. Complete deconvolution of cellular mixtures based on linearity of transcriptional signatures. Nat. Commun. 10, 2209 (2019). [GitHub](https://github.com/ctlab/LinSeed)|
| Coex  | Complete enrichment  | Kelley, K. W., Nakao-Inoue, H., Molofsky, A. V. & Oldham, M. C. Variation among intact tissue samples reveals the core transcriptional features of human CNS cell classes. Nat. Neurosci. 21, 265397 (2018). Our own implementation is provided in the scripts. |
| BrainInABlender  | Enrichment  | Hagenauer, M. H. et al. Inference of cell type content from human brain transcriptomic datasets illuminates the effects of age, manner of death, dissection, and psychiatric diagnosis. PLoS ONE 13, 89391 (2018). [GitHub](https://github.com/hagenaue/BrainInABlender)|
| xCell  | Enrichment  | Aran, D., Hu, Z. & Butte, A. J. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biol. 18, 1–14 (2017). [GitHub](https://github.com/dviraran/xCell)|

# Other packages
### For core analyses
In progress

### For auxiliary work, e.g. data wrangling and visualisation
In progress

# Directory structure
These scripts assume the following directory structure:

```
.
├── Data
│   ├── Preprocessed
│   └── Raw
├── Results
│   ├── ConfoundingComposition
│   │   ├── CompositionOnly
│   │   ├── CompositionAndExpression
│   │   └── CIBERSORTx
│   ├── BulkTissue
│   │   ├── Autism
│   │   ├── GTEx
│   │   └── Parikshak
│   ├── Other
│   │   ├── CompositionOnly
│   │   ├── CompositionAndExpression
│   │   └── CIBERSORTx
│   ├── BenchmarkingDatasets
│   │   ├── CA
│   │   ├── DM
│   │   ├── Immunopanned
│   │   ├── RNA
│   │   └── VL
└── Scripts
```



# Other recommended readings
For those interested in further reading about deconvolution benchmarking, we recommend the following as a starting point. These papers tend to focus on algorithm choice and data normalisation in non-brain tissues.

Avila Cobos, F., Alquicira-Hernandez, J., Powell, J.E. et al. Benchmarking of cell type deconvolution pipelines for transcriptomics data. Nat Commun 11, 5650 (2020). https://doi.org/10.1038/s41467-020-19015-1

Jin, H., Liu, Z. A benchmark for RNA-seq deconvolution analysis under dynamic testing environments. Genome Biol 22, 102 (2021). https://doi.org/10.1186/s13059-021-02290-6

Mohammadi, S., Zuckerman, N., Goldsmith, A., & Grama, A. A critical survey of deconvolution methods for separating cell types in complex tissues. Proceedings of the IEEE, 105(2), 340-366 (2017). https://doi.org/10.1109/JPROC.2016.2607121.
