# PBMC_AH_LPS_scRNA-seq

Code and additional data for single-cell RNA-seq studies from patient PBMCs performed in the Nagy Lab.

This repository contains all additional code and files used in the following manuscript:

A secondary immune surveillance pathway in human peripheral monocytes activated by co-regulated gene cassettes

Note:
For all scripts and files, the pathway is labeled as "XXXXXX"
This can be replaced with your pathway, ideally a directory that contains everything.
Within this directory, I have pre-labeled everything based on the sub-directories I used.
It is not required for your subdirectories to be the same, but you may need to create them.

List of scripts with descriptions:

Seurat_Cluster_1.R - R script used to create Seurat objects for all samples including publically available single cell RNA-seq PBMC data from 8 different platforms. Script also normalizes the data following the SCTransform algorithm and filters the data.
Input: Path to directory containing all outputs from cellranger (alignment files). Also path to directory containing publically available PBMC data.
Output: PBMC_Immune.integrated_scT_every_filt.Robj

Seurat_Cluster_2.R - R script used to perform clustering analyses, data integration (across all samples and public data), and cell cluster labeling.
Input: R object from Seurat_Cluster_1.R (PBMC_Immune.integrated_scT_every_filt.Robj)
Output: PBMC_Immune.combined_scT_every_filt.Robj

Seurat_Cluster_3.R - R script used to create all figures and perform differntial expression analyses.
Input: R object from Seurat_Cluster_2.R (PBMC_Immune.combined_scT_every_filt.Robj)
Output: All violin plots and clustering figures

List of files with descriptions:

PBMC_Immune.integrated_scT_every_filt_20191016.Robj
This is the R object output from Seurat_Cluster_1.R

PBMC_Immune.combined_scT_every_filt_20191017.Robj
This is the R object output from Seurat_Cluster_2.R - This file can be readily used to run Seurat_Cluster_3.R for any specific gene expression analyses and XXXXXX.R for any correlation analyses.

pbmc_ssc_mat.rds
pbmc_ssc_metadata.rds
Publically available PBMC data from https://satijalab.org/seurat/v3.0/integration.html


For sake of reducing redundancy, scripts and data for Bulk-RNA-seq analysis of Liver are available in a different Github Repository:
atomadam2/Liver_Deconvolution
LiverDeconv_BulkRNASleuth.R
