################
#load libraries#
################
library(Seurat)
library(ggplot2)
library(utils)
library(DESeq2)
library(grDevices)
library(cowplot)
library(dplyr)
library(patchwork)


###########################
#Basic pipeline for Seurat#
###########################
no_CDSUTR_Sporo_FUGI_R_D7465034=readRDS("no_CDSUTR_Sporo_FUGI_R_D7465034/QC_run2/QC_object_pass_only_cells.RDS")
no_CDSUTR_v9_Sporo_FUGI_R_D7465035=readRDS("no_CDSUTR_v9_Sporo_FUGI_R_D7465035/QC_run2/QC_object_pass_only_cells.RDS")
no_CDSUTR_Sporo_FUGI_R_D7465036=readRDS("no_CDSUTR_Sporo_FUGI_R_D7465036/QC_run2/QC_object_pass_only_cells.RDS")
no_CDSUTR_Sporo_FUGI_R_D7465037=readRDS("no_CDSUTR_Sporo_FUGI_R_D7465037/QC_run2/QC_object_pass_only_cells.RDS")

#merge
sporo.big <- merge(no_CDSUTR_Sporo_FUGI_R_D7465034, y = c(no_CDSUTR_v9_Sporo_FUGI_R_D7465035, no_CDSUTR_Sporo_FUGI_R_D7465036,no_CDSUTR_Sporo_FUGI_R_D7465037), add.cell.ids = c("FUGI_R_D7465034", "FUGI_R_D7465035", "FUGI_R_D7465036", "FUGI_R_D7465037"), project = "Sporos_FUGI")
sporo.big = NormalizeData(sporo.big)
sporo.big = FindVariableFeatures(sporo.big, selection.method = "vst", nfeatures = 2000) #there are different methods to find variable genes
#plot <- VariableFeaturePlot(sporo.big)

all.genes <- rownames(sporo.big)
sporo.big <- ScaleData(sporo.big, features = all.genes)
sporo.big <- RunPCA(sporo.big, features = VariableFeatures(object = sporo.big))
DimPlot(sporo.big, reduction = "pca")
ElbowPlot(sporo.big, ndims = 50)
#Check how many PCs to use for the data using Mathews Molecular cross validation
source("code.R")
mcv = molecularCrossValidation(sporo.big@assays$RNA@counts, varGenes = VariableFeatures(sporo.big),normalisation = minimalSeuratV3)

#So from this we have that the optimal number is 12
sporo.big <- FindNeighbors(sporo.big, dims = 1:12)
sporo.big <- FindClusters(sporo.big, resolution = 1)
sporo.big <- RunUMAP(sporo.big, dims = 1:12)
DimPlot(sporo.big, reduction = "umap")

#Find markers
sporo.big@misc$rocmarkers = FindAllMarkers(sporo.big,test.use='roc',only.pos=TRUE, min.pct = 0, return.thresh = 0)
sporo.big@misc$wilcoxmarkers = FindAllMarkers(sporo.big,only.pos=TRUE)



