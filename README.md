# BrainCellularComposition

The scripts folder contains the final R Code used for analysis in:

Sutton, Gavin J., Daniel Poppe, Rebecca K. Simmons, Kieran Walsh, Urwah Nawaz, Ryan Lister, Johann A. Gagnon-Bartsch, and Irina Voineagu. "Comprehensive evaluation of deconvolution methods for human brain gene expression." Nature Communications 13, no. 1 (2022): 1-18.

Files are as follows:
- Group 1: scripts that preprocess raw public and generated data 
- Group 2: estimating composition in simulated mixtures, and analyses of its output. Broadly corresponds to Figures 1-3. RME: In vitro RNA mixtures of cultured neuronal and astrocyte RNA. SCME: Mixtures derived by simulating pseudobulks from Darmanis et al.'s single cell RNA-seq. SNME: Mixtures derived by simulating pseudobulks from Velmeshev et al.'s or the Human Cell Atlas' single nucleus RNA-seq (VL and CA, respectively).
- Group 3: simulation of datasets with confounded composition between groups, exploring DE driven by this phenomenon. Broadly corresponds to Figures 4-5. CIBx: analyses using CIBERSORTx to impute cell-type-specific expression
- Group 4: analyses of public bulk brain RNA-seq from Parikshak et al., the GTEx consortium, and an autism spectrum disorder (ASD) dataset found in Parikshak et al. Broadly corresponds to Figure 6.
- Group 5: miscellaneous analysis scripts. CrossSignatureCorrelations: for plotting the relationship between signatures. CrossTissue: deconvolution of heart and pancreas data. Immunopanned: deconvolution of Zhang et al.'s immunopanned pure brain cell-types. RNA_Content: exploring whether deconvolution estimates cellular or RNA proportion
- Fun: contains general custom functions called in other scripts
