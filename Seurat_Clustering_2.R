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
sessionInfo()

setwd("~/XXXXXX/all_clust")
getwd()

########################################
##                                    ##
##        3. Clustering only          ##
##                                    ##
########################################

# Replace XXXXXX with path to a folder to store all final output files
# NOTE "XXXXXX" will be the same throughout this script, so all folders within are standardized
# The input for this script is the output from Seurat_Clustering_2.R also available in repository

load("I:/XXXXXX/PBMC_Immune.integrated_scT_every_filt_20191016.Robj")

# Clustering

Immune.integrated <- RunPCA(Immune.integrated, verbose = FALSE)
Immune.integrated <- RunUMAP(Immune.integrated, dims = 1:30)
Immune.integrated <- FindNeighbors(Immune.integrated)
Immune.integrated <- FindClusters(Immune.integrated,  resolution = 0.5)

## Make plots
# If run on the HPC, you will not see these, but it will work in RStudio
plots <- DimPlot(Immune.integrated, combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x + 
                  theme(legend.position = "top") + 
                  guides(color = guide_legend(nrow = 3, 
                                              byrow = TRUE, 
                                              override.aes = list(size = 3))))
CombinePlots(plots)

# Plots separating based on stimulation (LPS or Basal)
p1 <- DimPlot(Immune.integrated, reduction = "umap", group.by = "stim")
p2 <- DimPlot(Immune.integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

# Plots separating based on disease (AH or Healthy)
p1 <- DimPlot(Immune.integrated, reduction = "umap", group.by = "type")
p2 <- DimPlot(Immune.integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

# Create new meta data labels
Immune.integrated$phenogroup <- paste(Immune.integrated$type, Immune.integrated$stim, sep = "_")
Immune.integrated$JoinedClusters <- paste(Immune.integrated$seurat_clusters, Immune.integrated$CellType, sep = "_")

# Plots separating based on Publically Available data
p1 <- DimPlot(Immune.integrated, reduction = "umap", group.by = "CellType", label = TRUE)
p2 <- DimPlot(Immune.integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(Immune.integrated, reduction = "umap", group.by = "CellType", label = TRUE, label.size = 3,split.by = "phenogroup")
DimPlot(Immune.integrated, reduction = "umap", split.by = "phenogroup")

## Identify each cluster
DefaultAssay(Immune.integrated) <- "RNA"

## Cluster Markers
# Find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(Immune.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Create a file containing all marker genes 
# Caveat: The "FindAllMarkers" function is actually not great for finding definitive markers for cell clsuters
write.table(pbmc.markers, file="ClusterMarkers_PBMCscRNA_final.txt")


########################################
##                                    ##
##      4. Plots to characterize      ##
##                                    ##
########################################

## Plots for different genes and their expression in cells
VlnPlot(Immune.integrated, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"), 
        pt.size = 0.2, ncol = 4)
FeaturePlot(Immune.integrated, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), pt.size = 0.2,ncol = 3)
FeaturePlot(Immune.integrated, features = c("CCR7", "IL7R","GNLY", "GZMK", "CD8B", "CCL2","LTB","LYZ", "IFI44L"), min.cutoff = "q9")

markers.to.plot <- c("LEF1","TCF7","ITK",
                     "LEPROTL1","MAL","TNFRSF25",
                     "MCUB","TUT4","SMIM26", 
                     "C12orf75","GZMM","CCL4",
                     "TRGC1","TRGC2","LYAR",
                     "GNAI2","HMGB2","ETS1",
                     "PHACTR2","RNF125","TGFB1",
                     "HSPA1A","FABP5","MRC1",
                     "CSF3R","S100A8","CAPG",
                     "NAGK","HLA-DRB5","HLA-DQB1",
                     "CSF1R","LST1","MS4A7",
                     "FCMR","BANK1","MARCH1",
                     "SH2D1B","CMC1","HOPX",
                     "IRF4","PLD4","SERPINF1",
                     "RGS18","C2orf88","TMEM40")

DotPlot(Immune.integrated, features = rev(markers.to.plot), cols = c("blue","green","red"), dot.scale = 6, 
        split.by = "type", assay = 'SCT' ) + RotatedAxis()

## Name each cluster based on cell type gene expression and publically available data
# There are two versions of this, either - or _ between terms.

immune.combined <- RenameIdents(Immune.integrated, `0` = "CD4_T-cell1", `1` = "CD4_T-cell2", 
                                `2` = "B-cell",`3` = "NK-cell", `4` = "Cytotoxic_T-cell1", 
                                `5` = "CD14_Monocyte1",`6` = "CD14_Monocyte2", `7` = "Cytotoxic_T-cell2",
                                `8` = "Cytotoxic_T-cell3",`9` = "CD14_Monocyte3",
                                `10` = "Cytotoxic_T-cell4", `11` = "CD4_T-cell3",`12` = "CD16_Monocyte",
                                `13` = "pDC",`14` = "Megakaryocyte",
                                `15` = "Unassigned")
DimPlot(immune.combined, label = TRUE, label.size = 4)

#immune.combined <- RenameIdents(Immune.integrated, `0` = "CD4-T-cell1", `1` = "CD4-T-cell2", 
#                                `2` = "B-cell",`3` = "NK-cell", `4` = "Cytotoxic-T-cell1", 
#                                `5` = "CD14-Monocyte1",`6` = "CD14-Monocyte2", `7` = "Cytotoxic-T-cell2",
#                                `8` = "Cytotoxic-T-cell3",`9` = "CD14-Monocyte3",
#                                `10` = "Cytotoxic-T-cell4", `11` = "CD4-T-cell3",`12` = "CD16-Monocyte",
#                                `13` = "pDC",`14` = "Megakaryocyte",
#                                `15` = "Unassigned")
#DimPlot(immune.combined, label = TRUE, label.size = 4)

#Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("CD4-T-cell1","CD4-T-cell2","CD4-T-cell3",
#                                                                      "Cytotoxic-T-cell1", "Cytotoxic-T-cell2","Cytotoxic-T-cell3","Cytotoxic-T-cell4",
#                                                                      "CD14-Monocyte1","CD14-Monocyte2","CD14-Monocyte3","CD16-Monocyte",
#                                                                      "B-cell","NK-cell",
#                                                                      "pDC","Megakaryocyte",
#                                                                      "Unassigned"))

Immune.integrated$CompleteType <- paste(Immune.integrated$type, Immune.integrated$Experiment, sep = "-")
immune.combined$CompleteType <- paste(immune.combined$type, immune.combined$Experiment, sep = "-")

markers.to.plot <- c("LRRN3","ITGA6","TMEM204",
                     "LTB","ARHGAP15","NR3C1",
                     "LEF1","PRCKA","CD27",
                     "CD8A","CD99","CALCOCO1",
                     "LYAR","TRGC1","CHST12",
                     "FYN","GZMK","CMC1",
                     "KLRB1","CXCR6","AQP3",
                     "CCNL1","MAP3K1","STXBP2",
                     "CTSB","CD14","S100A12",
                     "HLA-DPA1","HLA-DQB1","HLA-DPB1",
                     "S100A11","PTPN6","FCGR3A",
                     "CD79A","LINC00926","MS4A1",
                     "PRF1","SPON2","IFNGR1",
                     "STMN1","CCDC50","MCM5",
                     "TAGLN2","MAP3K7CL","OAZ1",
                     "TXN","SOX4","CDK6")

# Creates dotplot figure - Not used in publication
pdf("DotPlotFig.pdf", height = 10, width = 15, useDingbats=FALSE)
DotPlot(immune.combined, features = rev(markers.to.plot), 
        cols = c("blue","red", "black", "black"), 
        dot.scale = 8, 
        split.by = "CompleteType",
        assay = 'SCT'
) + RotatedAxis()
dev.off()

markers.to.plot <- c("LRRN3","ITGA6",
                     "LTB","ARHGAP15",
                     "LEF1","PRKCA",
                     "CD8A","CALCOCO1",
                     "TRGC1","CHST12",
                     "GZMK","CMC1",
                     "KLRB1","CXCR6",
                     "MAP3K1","STXBP2",
                     "CD14","S100A12",
                     "HLA-DPA1","HLA-DQB1",
                     "PTPN6","FCGR3A",
                     "CD79A","LINC00926",
                     "PRF1","SPON2",
                     "STMN1","CCDC50",
                     "TAGLN2","MAP3K7CL",
                     "TXN","SOX4")

pdf("DotPlotFig_2.pdf", height = 10, width = 15, useDingbats=FALSE)
DotPlot(immune.combined, features = rev(markers.to.plot), 
        cols = c("blue","red", "black", "black"), 
        dot.scale = 8, 
        split.by = "CompleteType",
        assay = 'SCT'
) + RotatedAxis()
dev.off()

par(mfrow=c(1,2))

genes.viz=c("CD14","FCGR3A")
VlnPlot(immune.combined, genes.viz, pt.size = 0)


## Number cells in cluster
# Calculate number of cells per cluster from object@ident
table (Idents(immune.combined))
table(Idents(object = immune.combined))
immune.combined$celltype_origin <- paste(immune.combined$orig.ident, Idents(immune.combined), sep = "_")
cell.num <- as.data.frame(table(immune.combined$celltype_origin))
write.table(cell.num, file="CellTypes_PBMCscRNA.txt")


## Create a few meta data labels to help later
immune.combined$preserved <- Idents(immune.combined)

immune.combined$celltype.type <- paste(Idents(immune.combined), immune.combined$type, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.type"

Alcohol.response <- FindMarkers(immune.combined, ident.1 = "CD14_Monocyte1_Alcohol", ident.2 = "CD14_Monocyte1_Healthy", verbose = FALSE)
head(Alcohol.response, n = 15)

immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"

########################################
##                                    ##
##      5. Save R Objects Data        ##
##                                    ##
########################################

save(immune.combined,file="PBMC_Immune.combined_scT_every_filt.Robj")
