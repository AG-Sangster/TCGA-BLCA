# TCGA-BLCA

Run scripts in order of 
1. TCGA-BLCA-download.R
2. TCGA-BLCA-dataPrep.R
3. TCGA-BLCA-analysis.R

Check preferences at head of file before running. 


## DDR_CM_geneList.csv
Contains list of genes (hugo symbols) and the corresponding DNA Damage Response (DDR) and/or Chromatin Modification (CM) classification.

This file is required by TCGA-BLCA-dataPrep.R to identify DDR and CM genes. 

It must be in the same directory as the R script.



## TCGA-BLCA-download.R

This script downloads The Cancer Genome Atlas (TCGA) data for bladder cancer (BLCA). This data is in the format of a Mutation Annotation File (MAF) and there are 4 different MAF's provided by TCGA. This data is large and could take some time to download all 4 files. 

See comments in TCGA-BLCA-download.R to set preferences before running.

Please download and install required packages before running the script. 

Required packages: TCGAbiolinks



## TCGA-BLCA-dataPrep.R

This script prepares the downloaded data as it was done for the publication. The preprocessing is broken down into steps and each step is performed by an individual function. Running the script will do the preprocessing and could take some time to complete. 

This script can produce the figures from the publication that describe the data (coefficient of variation and venn diagram that compares the MAF mutation calls).

See comments in TCGA-BLCA-dataPrep.R to set preferences before running.

Note: the MAF files must be loaded into the current working environment and be named according to TCGA-BLCA-download.R

Please download and install required packages before running the script. 

Required packages: VennDiagram



## TCGA-BLCA-analysis.R

This script performs a statistical analysis of the TCGA-BLCA data. This process is also brokendown into functions and running this script will call them to perform the analysis, this could take some time to complete.

This script can also produce figures from the publication that describe the random sampling and observed values.

see comments in TCGA-BLCA-analysis.R to set preferences before running. 



## Environment Details

R version 4.0.0 (2020-04-24) -- "Arbor Day"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin17.0 (64-bit)
