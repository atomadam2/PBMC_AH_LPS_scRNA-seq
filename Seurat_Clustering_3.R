#!/usr/bin/Rscript

# Script formatted to run on either HPC or on computer

# Performs all clustering analysis and cell labeling

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
library(monocle3)
library(sctransform)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
sessionInfo()

setwd("~/XXXXXX/RPics")
getwd()

########################################
##                                    ##
##    3. Load and Organize Data       ##
##                                    ##
########################################

# Replace XXXXXX with path to a folder to store all final output files
# NOTE "XXXXXX" will be the same throughout this script, so all folders within are standardized
# The input for this script is the output from Seurat_Clustering_2.R also available in repository

load("I:/XXXXXX/PBMC_Immune.combined_scT_every_filt_20191017.Robj")

load("I:/Adam/Unix/scRNA-seq/PBMC_LPS_scRNA/every_clust/PBMC_Immune.combined_scT_every_filt_20191017.Robj")

Idents(immune.combined) <- "preserved"

# Organize the clusters the way I like
Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("CD4_T-cell1","CD4_T-cell2","CD4_T-cell3",
                                                                      "Cytotoxic_T-cell1", "Cytotoxic_T-cell2","Cytotoxic_T-cell3","Cytotoxic_T-cell4",
                                                                      "CD14_Monocyte1","CD14_Monocyte2","CD14_Monocyte3","CD16_Monocyte",
                                                                      "B-cell","NK-cell",
                                                                      "pDC","Megakaryocyte",
                                                                      "Unassigned"))
immune.combined$celltype <- factor(immune.combined$celltype, levels = c("CD4_T-cell1_Healthy","CD4_T-cell2_Healthy","CD4_T-cell3_Healthy",
                                                                        "Cytotoxic_T-cell1_Healthy", "Cytotoxic_T-cell2_Healthy",
                                                                        "Cytotoxic_T-cell3_Healthy","Cytotoxic_T-cell4_Healthy",
                                                                        "CD14_Monocyte1_Healthy","CD14_Monocyte2_Healthy",
                                                                        "CD14_Monocyte3_Healthy","CD16_Monocyte_Healthy",
                                                                        "B-cell_Healthy","NK-cell_Healthy",
                                                                        "pDC_Healthy","Megakaryocyte_Healthy",
                                                                        "Unassigned_Healthy",
                                                                        "CD4_T-cell1_Alcohol","CD4_T-cell2_Alcohol","CD4_T-cell3_Alcohol",
                                                                        "Cytotoxic_T-cell1_Alcohol", "Cytotoxic_T-cell2_Alcohol",
                                                                        "Cytotoxic_T-cell3_Alcohol","Cytotoxic_T-cell4_Alcohol",
                                                                        "CD14_Monocyte1_Alcohol","CD14_Monocyte2_Alcohol",
                                                                        "CD14_Monocyte3_Alcohol","CD16_Monocyte_Alcohol",
                                                                        "B-cell_Alcohol","NK-cell_Alcohol",
                                                                        "pDC_Alcohol","Megakaryocyte_Alcohol",
                                                                        "Unassigned_Alcohol",
                                                                        "CD4_T-cell1_NA","CD4_T-cell2_NA","CD4_T-cell3_NA",
                                                                        "Cytotoxic_T-cell1_NA", "Cytotoxic_T-cell2_NA",
                                                                        "Cytotoxic_T-cell3_NA","Cytotoxic_T-cell4_NA",
                                                                        "CD14_Monocyte1_NA","CD14_Monocyte2_NA",
                                                                        "CD14_Monocyte3_NA","CD16_Monocyte_NA",
                                                                        "B-cell_NA","NK-cell_NA",
                                                                        "pDC_NA","Megakaryocyte_NA",
                                                                        "Unassigned_NA"))

########################################
##                                    ##
##      3. Make Clustering Plot       ##
##                                    ##
########################################
# Create the UMAP Clustering Figure for

DimPlot(immune.combined, label = TRUE, label.size = 5)
# Define the number of colors you want
nb.cols <- 16
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

pdf("ClusteringFig1A.pdf", height = 5, width = 7, useDingbats=FALSE)
p <- DimPlot(immune.combined, label = TRUE, label.size = 3, cols = mycolors)
p
dev.off()

# Clustering check with public cell names
pdf("Clustering_extra.pdf", height = 5, width = 7, useDingbats=FALSE)
p <- DimPlot(immune.combined, reduction = "umap", group.by = "CellType", 
             label = TRUE, cols = mycolors, label.size = 3)
p
dev.off()


########################################
##                                    ##
##       3. Make Violin Plots         ##
##                                    ##
########################################
Idents(immune.combined) <- "preserved"

# If you want to look at the expression of your gene of interest, replace these genes here with your favorite.
# This will plot all cell types, organized by disease, with basal and LPS separated.
plots <- VlnPlot(immune.combined, features = c("IL1B", "CXCL8", "CXCL5"), split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

# Reorganized the lanes a different way 
Idents(immune.combined) <- "celltype.stim"
immune.combined$celltype <- factor(immune.combined$celltype, levels = c("CD4_T-cell1_Healthy","CD4_T-cell1_Alcohol",
                                                                        "CD4_T-cell2_Healthy","CD4_T-cell2_Alcohol",
                                                                        "CD4_T-cell3_Healthy","CD4_T-cell3_Alcohol",
                                                                        "Cytotoxic_T-cell1_Healthy","Cytotoxic_T-cell1_Alcohol",
                                                                        "Cytotoxic_T-cell2_Healthy","Cytotoxic_T-cell2_Alcohol",
                                                                        "Cytotoxic_T-cell3_Healthy","Cytotoxic_T-cell3_Alcohol",
                                                                        "Cytotoxic_T-cell4_Healthy","Cytotoxic_T-cell4_Alcohol",
                                                                        "CD14_Monocyte1_Healthy","CD14_Monocyte1_Alcohol",
                                                                        "CD14_Monocyte2_Healthy","CD14_Monocyte2_Alcohol",
                                                                        "CD14_Monocyte3_Healthy","CD14_Monocyte3_Alcohol",
                                                                        "CD16_Monocyte_Healthy","CD16_Monocyte_Alcohol",
                                                                        "B-cell_Healthy","B-cell_Alcohol",
                                                                        "NK-cell_Healthy","NK-cell_Alcohol",
                                                                        "pDC_Healthy","pDC_Alcohol",
                                                                        "Megakaryocyte_Healthy","Megakaryocyte_Alcohol",
                                                                        "Unassigned_Healthy","Unassigned_Alcohol",
                                                                        "CD4_T-cell1_NA",
                                                                        "CD4_T-cell2_NA",
                                                                        "CD4_T-cell3_NA",
                                                                        "Cytotoxic_T-cell1_NA",
                                                                        "Cytotoxic_T-cell2_NA",
                                                                        "Cytotoxic_T-cell3_NA",
                                                                        "Cytotoxic_T-cell4_NA",
                                                                        "CD14_Monocyte1_NA",
                                                                        "CD14_Monocyte2_NA",
                                                                        "CD14_Monocyte3_NA",
                                                                        "CD16_Monocyte_NA",
                                                                        "B-cell_NA",
                                                                        "NK-cell_NA",
                                                                        "pDC_NA",
                                                                        "Megakaryocyte_NA",
                                                                        "Unassigned_NA"))



DefaultAssay(immune.combined) <- "SCT"

pdf("ViolinFig2E_CXCL8CXCL5.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("CXCL8", "CXCL5"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

DefaultAssay(immune.combined) <- "SCT"
pdf("ViolinFig2E_CXCL2CXCL3.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("CXCL2", "CXCL3"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

DefaultAssay(immune.combined) <- "SCT"
pdf("ViolinFig2E_CXCL1EREG.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("CXCL1", "EREG"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()


pdf("ViolinFig2D_IL1BCCL2.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("IL1B", "CCL2"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()


pdf("ViolinFig2D_S100.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("S100A8", "S100A9"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()


pdf("ViolinFig2D_IFIT13.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("IFIT1", "IFIT3"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ViolinFig2D_IFIT3OAS1.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("IFIT3", "OAS1"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ViolinFig1E_CLEC4A7A.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("CLEC4A", "CLEC7A"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()


pdf("ViolinFig1DE_CLEC4E5A.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("CLEC4E", "CLEC5A"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ViolinFig1D_CD14CD16.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("CD14", "FCGR3A"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ViolinFig1DE_C3AR1C5AR1.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("C3AR1", "C5AR1"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ViolinFig1D_TLR2TLR4.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("TLR2", "TLR4"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ViolinFig1DE_CD86163.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("CD86", "CD163"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ViolinFig1E_Clec2B2D.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("CLEC2B", "CLEC2D"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ViolinFig1DE_TLR8Clec12A.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("TLR8", "CLEC12A"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE, idents = c('CD14_Monocyte1_Healthy_Basal','CD14_Monocyte2_Healthy_Basal',
                                                          'CD14_Monocyte3_Healthy_Basal','CD16_Monocyte_Healthy_Basal',
                                                          'CD14_Monocyte1_Healthy_LPS','CD14_Monocyte2_Healthy_LPS',
                                                          'CD14_Monocyte3_Healthy_LPS','CD16_Monocyte_Healthy_LPS',
                                                          'CD14_Monocyte1_Alcohol_Basal','CD14_Monocyte2_Alcohol_Basal',
                                                          'CD14_Monocyte3_Alcohol_Basal','CD16_Monocyte_Alcohol_Basal',
                                                          'CD14_Monocyte1_Alcohol_LPS','CD14_Monocyte2_Alcohol_LPS',
                                                          'CD14_Monocyte3_Alcohol_LPS','CD16_Monocyte_Alcohol_LPS'))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()


########################################
##                                    ##
##    4. Differential Expression      ##
##                                    ##
########################################

# DE measurements of all clusters comparing the following:
# AH Basal vs HC Basal
# HC LPS vs HC Basal
# AH LPS  vs AH Basal
# AH LPS vs HC LPS
# Differential expression was performed using a loop to run through all clusters
# Output are text files
# Measurements are made using logratio
# All samples are paired

# DE can be done by other algorithms including MAST, but the following is needed:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MAST")

Idents(immune.combined) <- "orig.ident"

# This definies each of the paired samples
immune.combined <- RenameIdents(immune.combined, 
                                `pbmc_01_HC_B` = "Pair1", `pbmc_02_HC_L` = "Pair1", `pbmc_03_AH_B` = "Pair1", `pbmc_04_AH_L` = "Pair1",
                                `pbmc_05_HC_B` = "Pair2", `pbmc_06_HC_L` = "Pair2", `pbmc_07_AH_B` = "Pair2", `pbmc_08_AH_L` = "Pair2",
                                `pbmc_09_HC_B` = "Pair3", `pbmc_10_HC_L` = "Pair3", `pbmc_11_AH_B` = "Pair3", `pbmc_12_AH_L` = "Pair3",
                                `pbmc_13_HC_B` = "Pair4", `pbmc_14_HC_L` = "Pair4", `pbmc_15_AH_B` = "Pair4", `pbmc_16_AH_L` = "Pair4")
immune.combined$samplesource <- Idents(immune.combined)

Idents(immune.combined) <- "celltype.stim"
DefaultAssay(immune.combined) <- "SCT"
for (spec_cell in c("CD4_T-cell1","CD4_T-cell2","CD4_T-cell3",
                    "Cytotoxic_T-cell1", "Cytotoxic_T-cell2","Cytotoxic_T-cell3","Cytotoxic_T-cell4",
                    "CD14_Monocyte1","CD14_Monocyte2","CD14_Monocyte3","CD16_Monocyte",
                    "B-cell","NK-cell",
                    "pDC","Megakaryocyte",
                    "Unassigned")){
  #for (spec_cell in c("CD14_Monocyte1","CD14_Monocyte2","CD14_Monocyte3","CD16_Monocyte"
  #                    )){
  Basal_H <- (paste(spec_cell, "_Healthy_Basal", sep = ""))
  LPS_H <- (paste(spec_cell, "_Healthy_LPS", sep = ""))
  Basal_A <- (paste(spec_cell, "_Alcohol_Basal", sep = ""))
  LPS_A <- (paste(spec_cell, "_Alcohol_LPS", sep = ""))
  
  Alcohol.response <- FindMarkers(immune.combined, ident.1 = Basal_A, ident.2 = Basal_H, verbose = FALSE, test.use = "LR", latent.vars = "samplesource")
  Basal_Comp_File <- (paste(spec_cell, "_Basal_PBMC_scRNA.txt", sep = ""))
  write.table(Alcohol.response, file=Basal_Comp_File)
  
  LPS.H.response <- FindMarkers(immune.combined, ident.1 = LPS_H, ident.2 = Basal_H, verbose = FALSE, test.use = "LR", latent.vars = "samplesource")
  HLPS_Comp_File <- (paste(spec_cell, "_HealthyLPS_PBMC_scRNA.txt", sep = ""))
  write.table(LPS.H.response, file=HLPS_Comp_File)
  
  LPS.A.response <- FindMarkers(immune.combined, ident.1 = LPS_A, ident.2 = Basal_A, verbose = FALSE, test.use = "LR", latent.vars = "samplesource")
  ALPS_Comp_File <- (paste(spec_cell, "_AlcoholLPS_PBMC_scRNA.txt", sep = ""))
  write.table(LPS.A.response, file=ALPS_Comp_File)
  
  LPSallresponse <- FindMarkers(immune.combined, ident.1 = LPS_A, ident.2 = LPS_H, verbose = FALSE, test.use = "LR", latent.vars = "samplesource")
  LPSall_Comp_File <- (paste(spec_cell, "_LPSall_PBMC_scRNA.txt", sep = ""))
  write.table(LPSallresponse, file=LPSall_Comp_File)
}


########################################
##                                    ##
##     5. Supplemental Figure 3       ##
##                                    ##
########################################

# The following are the figures for Supplemental Figure 3
# This shows the different sequencing stats for each cell type.

Idents(immune.combined) <- "celltype.type"
DefaultAssay(immune.combined) <- "SCT"

VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plots <- VlnPlot(immune.combined, features = c("nFeature_RNA","nCount_RNA","percent.mt"), 
                 split.by = "stim", group.by = "celltype", pt.size = 0, combine = FALSE) + 
  stat_summary(fun.y="median", geom="point") +
  geom_boxplot(width=.3) 
CombinePlots(plots = plots, ncol = 1)

VlnPlot(immune.combined, features = c("nCount_RNA"), split.by = "stim", group.by = "celltype", pt.size = 0,
        idents = c("CD4_T-cell1_Healthy","CD4_T-cell1_Alcohol",
                   "CD4_T-cell2_Healthy","CD4_T-cell2_Alcohol",
                   "CD4_T-cell3_Healthy","CD4_T-cell3_Alcohol",
                   "Cytotoxic_T-cell1_Healthy","Cytotoxic_T-cell1_Alcohol",
                   "Cytotoxic_T-cell2_Healthy","Cytotoxic_T-cell2_Alcohol",
                   "Cytotoxic_T-cell3_Healthy","Cytotoxic_T-cell3_Alcohol",
                   "Cytotoxic_T-cell4_Healthy","Cytotoxic_T-cell4_Alcohol",
                   "CD14_Monocyte1_Healthy","CD14_Monocyte1_Alcohol",
                   "CD14_Monocyte2_Healthy","CD14_Monocyte2_Alcohol",
                   "CD14_Monocyte3_Healthy","CD14_Monocyte3_Alcohol",
                   "CD16_Monocyte_Healthy","CD16_Monocyte_Alcohol",
                   "B-cell_Healthy","B-cell_Alcohol",
                   "NK-cell_Healthy","NK-cell_Alcohol",
                   "pDC_Healthy","pDC_Alcohol",
                   "Megakaryocyte_Healthy","Megakaryocyte_Alcohol",
                   "Unassigned_Healthy","Unassigned_Alcohol")) + 
  stat_summary(fun.y="median", geom="point") +
  #  stat_summary(fun.data = quantiles_95, geom="boxplot") +
  geom_boxplot(width=.3) 

VlnPlot(immune.combined, features = c("percent.mt"), split.by = "stim", group.by = "celltype", pt.size = 0,
        idents = c("CD4_T-cell1_Healthy","CD4_T-cell1_Alcohol",
                   "CD4_T-cell2_Healthy","CD4_T-cell2_Alcohol",
                   "CD4_T-cell3_Healthy","CD4_T-cell3_Alcohol",
                   "Cytotoxic_T-cell1_Healthy","Cytotoxic_T-cell1_Alcohol",
                   "Cytotoxic_T-cell2_Healthy","Cytotoxic_T-cell2_Alcohol",
                   "Cytotoxic_T-cell3_Healthy","Cytotoxic_T-cell3_Alcohol",
                   "Cytotoxic_T-cell4_Healthy","Cytotoxic_T-cell4_Alcohol",
                   "CD14_Monocyte1_Healthy","CD14_Monocyte1_Alcohol",
                   "CD14_Monocyte2_Healthy","CD14_Monocyte2_Alcohol",
                   "CD14_Monocyte3_Healthy","CD14_Monocyte3_Alcohol",
                   "CD16_Monocyte_Healthy","CD16_Monocyte_Alcohol",
                   "B-cell_Healthy","B-cell_Alcohol",
                   "NK-cell_Healthy","NK-cell_Alcohol",
                   "pDC_Healthy","pDC_Alcohol",
                   "Megakaryocyte_Healthy","Megakaryocyte_Alcohol",
                   "Unassigned_Healthy","Unassigned_Alcohol")) + 
  stat_summary(fun.y="median", geom="point") +
  #  stat_summary(fun.data = quantiles_95, geom="boxplot") +
  geom_boxplot(width=.3 ) 

VlnPlot(immune.combined, features = c("nFeature_RNA"), split.by = "stim", group.by = "celltype", pt.size = 0,
        idents = c("CD4_T-cell1_Healthy","CD4_T-cell1_Alcohol",
                   "CD4_T-cell2_Healthy","CD4_T-cell2_Alcohol",
                   "CD4_T-cell3_Healthy","CD4_T-cell3_Alcohol",
                   "Cytotoxic_T-cell1_Healthy","Cytotoxic_T-cell1_Alcohol",
                   "Cytotoxic_T-cell2_Healthy","Cytotoxic_T-cell2_Alcohol",
                   "Cytotoxic_T-cell3_Healthy","Cytotoxic_T-cell3_Alcohol",
                   "Cytotoxic_T-cell4_Healthy","Cytotoxic_T-cell4_Alcohol",
                   "CD14_Monocyte1_Healthy","CD14_Monocyte1_Alcohol",
                   "CD14_Monocyte2_Healthy","CD14_Monocyte2_Alcohol",
                   "CD14_Monocyte3_Healthy","CD14_Monocyte3_Alcohol",
                   "CD16_Monocyte_Healthy","CD16_Monocyte_Alcohol",
                   "B-cell_Healthy","B-cell_Alcohol",
                   "NK-cell_Healthy","NK-cell_Alcohol",
                   "pDC_Healthy","pDC_Alcohol",
                   "Megakaryocyte_Healthy","Megakaryocyte_Alcohol",
                   "Unassigned_Healthy","Unassigned_Alcohol")) + 
  stat_summary(fun.y="median", geom="point") +
  #  stat_summary(fun.data = quantiles_95, geom="boxplot") +
  geom_boxplot(width=.3 ) 


