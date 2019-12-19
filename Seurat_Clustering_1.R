#!/usr/bin/Rscript

# Script formatted to run on HPC

# Creates Seurat Objects for all samples including
# 16 PBMC single cell datasets
# Publically available single cell data from 8 platforms

# Normalizes data by the SCTransform algorithm and filters
# Combines all data into a Seurat object

########################################
##                                    ##
##   1. Install programs needed       ##
##                                    ##
########################################

# If your HPC is already set up with all these packages, this is not necessary
# But it depends on your system

### You will need this to update R ###
#install.packages('Seurat')

### Need python?? ###
#library(reticulate)
#py_config()
# create a new environment
#conda_create("r-reticulate")

#install.packages("installr")
#install.packages("stringr")

#library(installr)

#install.packages("dplyr")
#install.packages("cowplot")


########################################
##                                    ##
##         2. Load packages           ##
##                                    ##
########################################

library(Seurat)
library(cowplot)
library(dplyr)
library(sctransform)

packageVersion("Seurat")

# Replace XXXXXX with path to a folder to store all final output files
# NOTE "XXXXXX" will be the same throughout this script, so all folders within are standardized

setwd("~/XXXXXX/all_clust")
getwd()

########################################
##                                    ##
##          3. Upload Data            ##
##                                    ##
########################################

## Import publically available single cell RNA-seq data from PBMCs performed with 8 different platforms

# Data are available here:
# https://satijalab.org/seurat/v3.0/integration.html
# But within this repository, it is also available (to ensure data are preserved)

# Replace XXXXXX with path to data
# Prior to running this script, download these files and put it into the correct directory

pbmc.data <- readRDS("~/XXXXXX/PBMC_Compiled_Data/pbmc_ssc_mat.rds")
pbmc.metadata <- readRDS("~/XXXXXX/PBMC_Compiled_Data/pbmc_ssc_metadata.rds")

# Create Seurat Object of publically available data
pbmc <- CreateSeuratObject(counts = pbmc.data, meta.data = pbmc.metadata)

# Subset data to remove cells with low reads
pbmc <- subset(pbmc, subset = nFeature_RNA > 200)

# Create a list of all the different methods used in the publically available dataset (8 in total)
pbmc.list <- SplitObject(pbmc, split.by = "Method")
for (i in names(pbmc.list)) {
  pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = TRUE)
}

## Load the PBMC dataset
# For this, I created a function that performs all steps including:
# Labeling the samples based on treatment and patient disease
# Getting the data
# Creating Seurat object
# Filtering
# SCTransform algorithm for normalization

# Replace XXXXXX with path to directory that contains the data (output from cellranger)

SeuratAll <- function(x){
  print(x)
  spec_sample <- x
  # Generate the correct file name
  filename <- (paste("~/XXXXXX/ind_analysis",
                     spec_sample,
                     "outs/filtered_feature_bc_matrix/", sep = "/"))
  print(filename)
  
  # Determine what each sample treatment is based on file name
  testvec <- unlist(strsplit(spec_sample, "_"))
  patientval <- testvec[3]
  stimval <- testvec[4]
  if (patientval=="HC"){
      finalpat <- "Healthy"
  }
  else {
      finalpat <- "Alcohol"
  }
  if (stimval=="B"){
    finalstim <- "Basal"
  }
  else {
    finalstim <- "LPS"
  }
  
  # Import data
  SeuratSample <- Read10X(data.dir = filename)
  
  # Set up control
  Seurat_temp <- CreateSeuratObject(counts = SeuratSample , project = spec_sample)
  
  # Store mitochondrial percentage in object meta data
  Seurat_temp <- PercentageFeatureSet(Seurat_temp, pattern = "^MT-", col.name = "percent.mt")
  Seurat_temp$type <- finalpat
  Seurat_temp$stim <- finalstim
  Seurat_temp$CellType <- "none"
  Seurat_temp$Method <- "none"
  Seurat_temp$Experiment <- "none"
  
  # Look at numbers - seems similar +/- aggr 
  # If run on the HPC, you will not see these, but it will work in RStudio
  VlnPlot(Seurat_temp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot1 <- FeatureScatter(Seurat_temp, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(Seurat_temp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  
  # Run sctransform
  Seurat_temp <- SCTransform(Seurat_temp, vars.to.regress = "percent.mt", verbose = TRUE)
  # Filter cells with low reads, potential doublets, and high mitochondrial reads  
  Seurat_temp <- subset(Seurat_temp, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
  return(Seurat_temp)
}

## Run function to upload and normailze for all PBMC samples
HC_B_01_Seurat <- SeuratAll("pbmc_01_HC_B")
HC_L_02_Seurat <- SeuratAll("pbmc_02_HC_L")
AH_B_03_Seurat <- SeuratAll("pbmc_03_AH_B")
AH_L_04_Seurat <- SeuratAll("pbmc_04_AH_L")

HC_B_05_Seurat <- SeuratAll("pbmc_05_HC_B")
HC_L_06_Seurat <- SeuratAll("pbmc_06_HC_L")
AH_B_07_Seurat <- SeuratAll("pbmc_07_AH_B")
AH_L_08_Seurat <- SeuratAll("pbmc_08_AH_L")

HC_B_09_Seurat <- SeuratAll("pbmc_09_HC_B")
HC_L_10_Seurat <- SeuratAll("pbmc_10_HC_L")
AH_B_11_Seurat <- SeuratAll("pbmc_11_AH_B")
AH_L_12_Seurat <- SeuratAll("pbmc_12_AH_L")

HC_B_13_Seurat <- SeuratAll("pbmc_13_HC_B")
HC_L_14_Seurat <- SeuratAll("pbmc_14_HC_L")
AH_B_15_Seurat <- SeuratAll("pbmc_15_AH_B")
AH_L_16_Seurat <- SeuratAll("pbmc_16_AH_L")

# Combine into list
sample_all <- c(HC_B_01_Seurat,HC_L_02_Seurat,AH_B_03_Seurat,AH_L_04_Seurat,
                HC_B_05_Seurat,HC_L_06_Seurat,AH_B_07_Seurat,AH_L_08_Seurat,
                HC_B_09_Seurat,HC_L_10_Seurat,AH_B_11_Seurat,AH_L_12_Seurat,
                HC_B_13_Seurat,HC_L_14_Seurat,AH_B_15_Seurat,AH_L_16_Seurat, pbmc.list)

########################################
##                                    ##
##        4. Integrate Data           ##
##                                    ##
########################################

# Use the  reduction = "rpca" for all cells everywhere

Immune.features <- SelectIntegrationFeatures(object.list = sample_all, nfeatures = 3000)
options(future.globals.maxSize = 4800 * 1024^2)

Immune.list <- PrepSCTIntegration(object.list = sample_all, 
	       				      anchor.features = Immune.features, 
                                    	      verbose = TRUE)
Immune.list <- lapply(X = Immune.list, FUN = RunPCA, verbose = TRUE, features = Immune.features)
Immune.anchors <- FindIntegrationAnchors(object.list = Immune.list, 
                                         normalization.method = "SCT",
                                         reduction = "rpca",
                                         anchor.features = Immune.features,
                                         verbose = TRUE)


########################################
##                                    ##
##      5. Save R Objects Data        ##
##                                    ##
########################################

save(Immune.anchors,file="PBMC_Immune.anchors_scT_every_20191016.Robj")

Immune.integrated <- IntegrateData(anchorset = Immune.anchors, normalization.method = "SCT",
                     verbose = TRUE)

save(Immune.integrated,file="PBMC_Immune.integrated_scT_every_filt.Robj")

