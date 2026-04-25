# MRTrios_analysis
This repo consists of R scripts that aim in causal network inference with Mendelian randomization for trios consisting of copy number alteration, gene expression and DNA methylation.  The R scripts here perform preprocessing, causal inference (using functions in the R package MRTrios: https://github.com/audreyqyfu/MRTrios), and downstream analyses. The user will need to change the paths and names of the files mentioned in these scripts to suit their own needs.

The scripts should be used in the following order:

- `DataProcessing.R`: this file performs logit transformation of the methylation data and extracts ER+ and ER- patients from the clinical dataset.    
- `mainTrioMatch.R`: this file generates the trio data matrix by integrating the CNA, methylation, and gene expression data, with each trio in a separate line and each line containing the row numbers of the probe or gene in the input data.
- `trio.gene.type.sep.R`: this file takes the trios data matrix as input and keeps only "protein coding genes" or "lncRNAs".
- `main.findPCs.R`: this file calculates the principal component (PC) score matrix for methylation and separately for gene expression data, and identifies PCs that are significantly associated each trio.    
- `main.analyzeTrios.R`: this file performs causal network inference for trios and their associated confounders (i.e., PCs) to infer the causal models.   
- `HumanMethProbeInfo.R`: this file extracts the probe information like mapinfo, gene start/end, etc for a specific model type.
- `trio.location.Prcoding.lncRNA.R`: this file extracts the location of methylation probe in the individual genes in each trio.   
- `GOanalysis.R`: this file performs the gene ontology (GO) enrichment analysis for the mediation trios.
- `Methyl_Figures.R`: R script for generating publication-quality visualizations of DNA methylation data, including distance distributions, genomic feature analysis, and comparative plots across breast cancer subtypes and statistical models.
- `get_trio_data.R`: Extracts molecular data (methylation, CNA, and gene expression) for any specified trio of ER+/ER- patients and infers its causal model using infer.trio() from the MRGN package.    
- `get_trio_models.R`: Calculate all trios causal models, Loop over all trios and collect results( combined the trios with models info).

Additional scripts may be used for debugging:
- `model loc data extract.R`: this file extracts row numbers in the trio data matrix when a specific model type and location (e.g., TSS1500, body, 3' UTR, etc.) of a gene is provided.


  
