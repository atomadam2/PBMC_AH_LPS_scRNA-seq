
# Script used for single-cell gene correlation analyses

# Requires R object output from Seurat_Clustering_2.R
# Requires Bigscale

########################################
##                                    ##
##        1. Install Packages         ##
##                                    ##
########################################

#BIGscale

#Install and open packages
install.packages('devtools')
library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BioQC")

devtools::install_github("iaconogi/bigSCale2")

library("bigSCale", lib.loc="~/R/win-library/3.6")

# Convert Seurat data to SingeCellExperiment Type
# install scater https://bioconductor.org/packages/release/bioc/html/scater.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scater")
library(scater)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')


#install.packages("igraph")
## Load package
library(igraph)


########################################
##                                    ##
##         2. Load PAckages           ##
##                                    ##
########################################

library("bigSCale", lib.loc="~/R/win-library/3.6")
library(scater)
library(loomR)
library(Seurat)
library(ggplot2)
library(reshape2)
library(igraph)
library(cowplot)

########################################
##                                    ##
##         3. Load R object           ##
##                                    ##
########################################

# Replace XXXX with path to Seurat object from Seurat_Clustering_2.R script
load("XXXX/PBMC_Immune.combined_scT_every_filt_20191017.Robj")

setwd("")
getwd()


########################################
##                                    ##
##      4. Subset Monocyte Data       ##
##                                    ##
########################################

HC_B_Mono_c <- subset(immune.combined, idents = c("CD14_Monocyte1_Healthy_Basal","CD14_Monocyte2_Healthy_Basal",
                                                  "CD14_Monocyte3_Healthy_Basal","CD16_Monocyte_Healthy_Basal"))
HC_L_Mono_c <- subset(immune.combined, idents = c("CD14_Monocyte1_Healthy_LPS","CD14_Monocyte2_Healthy_LPS",
                                                  "CD14_Monocyte3_Healthy_LPS","CD16_Monocyte_Healthy_LPS"))
AH_B_Mono_c <- subset(immune.combined, idents = c("CD14_Monocyte1_Alcohol_Basal","CD14_Monocyte2_Alcohol_Basal",
                                                  "CD14_Monocyte3_Alcohol_Basal","CD16_Monocyte_Alcohol_Basal"))
AH_L_Mono_c <- subset(immune.combined, idents = c("CD14_Monocyte1_Alcohol_LPS","CD14_Monocyte2_Alcohol_LPS",
                                                  "CD14_Monocyte3_Alcohol_LPS","CD16_Monocyte_Alcohol_LPS"))

HC_B_Mono_c <- FindVariableFeatures(HC_B_Mono_c, selection.method = "vst", nfeatures = 2000)
HC_L_Mono_c <- FindVariableFeatures(HC_L_Mono_c, selection.method = "vst", nfeatures = 2000)
AH_B_Mono_c <- FindVariableFeatures(AH_B_Mono_c, selection.method = "vst", nfeatures = 2000)
AH_L_Mono_c <- FindVariableFeatures(AH_L_Mono_c, selection.method = "vst", nfeatures = 2000)


########################################
##                                    ##
##       5. Gene Correlation          ##
##                                    ##
########################################

# The correlations are beasured here
# This takes a while
# This is what was used in the manuscript
# None of the more interesting network analyses were put into the manuscript
# But it is all interesting and work to do

DefaultAssay(HC_B_Mono_c) <- "SCT"
HC_B.data <- as.data.frame(as.matrix(GetAssayData(HC_B_Mono_c, slot = "counts")))
gene.names <- rownames(HC_B.data)
HC_B.results=compute.network(expr.data = HC_B.data,
                             gene.names = gene.names
                             #clustering = 'direct'
)
HC_B.results$graph
DT::datatable(HC_B.results$centrality)
#write.table(HC_B.results$centrality, file="Monocyte_Network_HC_Basal.txt")


DefaultAssay(HC_L_Mono_c) <- "SCT"
HC_L.data <- as.data.frame(as.matrix(GetAssayData(HC_L_Mono_c, slot = "counts")))
gene.names <- rownames(HC_L.data)
HC_L.results=compute.network(expr.data = HC_L.data,
                             gene.names = gene.names
                             #clustering = 'direct'
)
HC_L.results$graph
DT::datatable(HC_L.results$centrality)
#write.table(HC_L.results$centrality, file="Monocyte_Network_HC_LPS.txt")

DefaultAssay(AH_B_Mono_c) <- "SCT"
AH_B.data <- as.data.frame(as.matrix(GetAssayData(AH_B_Mono_c, slot = "counts")))
gene.names <- rownames(AH_B.data)
AH_B.results=compute.network(expr.data = AH_B.data,
                             gene.names = gene.names
                             #clustering = 'direct'
)
AH_B.results$graph
DT::datatable(AH_B.results$centrality)
#write.table(AH_B.results$centrality, file="Monocyte_Network_AH_Basal.txt")

DefaultAssay(AH_L_Mono_c) <- "SCT"
AH_L.data <- as.data.frame(as.matrix(GetAssayData(AH_L_Mono_c, slot = "counts")))
gene.names <- rownames(AH_L.data)
AH_L.results=compute.network(expr.data = AH_L.data,
                             gene.names = gene.names
                             #clustering = 'direct'
)
AH_L.results$graph
DT::datatable(AH_L.results$centrality)
#write.table(AH_L.results$centrality, file="Monocyte_Network_AH_LPS.txt")

######################################################
# The following section was not used but can be used for network analyses
# This is mostly untested code as I still need time to  work on it
nodeinterest <- c("CLEC5A", "CLEC7A", "CLEC4D", "CLEC6A", "CLEC4A", "KLRG1")
nodeinterest <- c("CXCL8", "CXCL5", "CXCL6", "CXCL1", "CXCL3", "CXCL2")
nodeinterest <- c("CXCL8", "SOCS1", "CLEC4E", "MALAT1", "IL1B", "FUCA1")

HC_B.results$name <- rownames(HC_B.results$centrality)
plot.igraph(HC_B.results$graph, axes = FALSE, add = FALSE, xlim = c(-1, 1),
            vertex.size = 1,
            vertex.label = ifelse(HC_B.results$name %in% nodeinterest, HC_B.results$name, NA),
            vertex.label.font= 1,  
            vertex.label.cex = 1,
            ylim = c(-1, 1))

HC_L.results$name <- rownames(HC_L.results$centrality)
plot.igraph(HC_L.results$graph, axes = FALSE, add = FALSE, xlim = c(-1, 1),
            vertex.size = 1,
            vertex.label = ifelse(HC_L.results$name %in% nodeinterest, HC_L.results$name, NA),
            vertex.label.font= 1,  
            vertex.label.cex = 1,
            ylim = c(-1, 1))

AH_B.results$name <- rownames(AH_B.results$centrality)
plot.igraph(AH_B.results$graph, axes = FALSE, add = FALSE, xlim = c(-1, 1),
            vertex.size = 1,
            vertex.label = ifelse(AH_B.results$name %in% nodeinterest, AH_B.results$name, NA),
            vertex.label.font= 1,  
            vertex.label.cex = 1,
            ylim = c(-1, 1))

AH_L.results$name <- rownames(AH_L.results$centrality)
plot.igraph(AH_L.results$graph, axes = FALSE, add = FALSE, xlim = c(-1, 1),
            vertex.size = 1,
            vertex.label = ifelse(AH_L.results$name %in% nodeinterest, AH_L.results$name, NA),
            vertex.label.font= 1,  
            vertex.label.cex = 1,
            ylim = c(-1, 1))
######################################################

########################################
##                                    ##
##         6. Prepare data for        ##
##            subsets/figures         ##
########################################

HCBmatrix<-as.data.frame(as.numeric(HC_B.results$correlations))
HCLmatrix<-as.data.frame(as.numeric(HC_L.results$correlations))
AHBmatrix<-as.data.frame(as.numeric(AH_B.results$correlations))
AHLmatrix<-as.data.frame(as.numeric(AH_L.results$correlations))

# You are my density
allHCB<-melt(HCBmatrix)
#HCB <- plot(density(allHCB$value))
allHCL<-melt(HCLmatrix)
#HCL <- plot(density(allHCL$value))
allAHB<-melt(AHBmatrix)
#AHB <- plot(density(allAHB$value))
allAHL<-melt(AHLmatrix)
#AHL <- plot(density(allAHL$value))


# NKC and CXC Figures

# These are lists of genes of interest
# I have provided these files to generate figures in the paper and more
# These lists are organized such that all genes are in genome position order

# Replace XXXX with path to .csv file available at github
targets <- read.csv("XXXX/NKC_CXC_combo.csv", header = FALSE)
NKCtargets <-targets$V1 
rm(targets)

# HC Basal Figure 

matrixtest <- HCBmatrix
matrixall_rank_all <- allHCB

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed3$gene <- rownames(matrix_subbed3)
NKC_all <- melt(matrix_subbed3)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed3$gene))

NKC_rank_all <- NKC_all

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

pdf("HCB_NKC_CXC_Fig6A.pdf", height = 7, width = 7, useDingbats=FALSE)
p
dev.off()

# HC LPS Figure 

matrixtest <- HCLmatrix
matrixall_rank_all <- allHCL

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed3$gene <- rownames(matrix_subbed3)
NKC_all <- melt(matrix_subbed3)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed3$gene))
NKC_rank_all <- NKC_all

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

pdf("HCL_NKC_CXC_Fig6A.pdf", height = 7, width = 7, useDingbats=FALSE)
p
dev.off()

# AH Basal Figure

matrixtest <- AHBmatrix
matrixall_rank_all <- allAHB

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed3$gene <- rownames(matrix_subbed3)
NKC_all <- melt(matrix_subbed3)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed3$gene))

NKC_rank_all <- NKC_all

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

pdf("AHB_NKC_CXC_Fig6A.pdf", height = 7, width = 7, useDingbats=FALSE)
p
dev.off()

# AH LPS Figure

matrixtest <- AHLmatrix
matrixall_rank_all <- allAHL

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed3$gene <- rownames(matrix_subbed3)
NKC_all <- melt(matrix_subbed3)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed3$gene))

NKC_rank_all <- NKC_all

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.95))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.95))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.05))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.05))]

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

pdf("AHL_NKC_CXC_Fig6A.pdf", height = 7, width = 7, useDingbats=FALSE)
p
dev.off()

# CC Figures

# These are lists of genes of interest
# I have provided these files to generate figures in the paper and more
# These lists are organized such that all genes are in genome position order


# Replace XXXX with path to .csv file available at github
targets <- read.csv("XXXX/CCL_genes.csv", header = FALSE)
NKCtargets <-targets$V1 
rm(targets)

# HC Basal Figure 

matrixtest <- HCBmatrix
matrixall_rank_all <- allHCB

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed3$gene <- rownames(matrix_subbed3)
NKC_all <- melt(matrix_subbed3)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed3$gene))

NKC_rank_all <- NKC_all

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

pdf("HCB_CC_SuppFig3.pdf", height = 7, width = 7, useDingbats=FALSE)
p
dev.off()

# HC LPS Figure 

matrixtest <- HCLmatrix
matrixall_rank_all <- allHCL

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed3$gene <- rownames(matrix_subbed3)
NKC_all <- melt(matrix_subbed3)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed3$gene))
NKC_rank_all <- NKC_all

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

pdf("HCL_CC_SuppFig3.pdf", height = 7, width = 7, useDingbats=FALSE)
p
dev.off()

# AH Basal Figure

matrixtest <- AHBmatrix
matrixall_rank_all <- allAHB

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed3$gene <- rownames(matrix_subbed3)
NKC_all <- melt(matrix_subbed3)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed3$gene))

NKC_rank_all <- NKC_all

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

pdf("AHB_CC_SuppFig3.pdf", height = 7, width = 7, useDingbats=FALSE)
p
dev.off()

# AH LPS Figure

matrixtest <- AHLmatrix
matrixall_rank_all <- allAHL

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed3$gene <- rownames(matrix_subbed3)
NKC_all <- melt(matrix_subbed3)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed3$gene))

NKC_rank_all <- NKC_all

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

pdf("AHL_CC_SuppFig3.pdf", height = 7, width = 7, useDingbats=FALSE)
p
dev.off()

