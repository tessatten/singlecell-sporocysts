library(Seurat)
library(ggplot2)
sporo=readRDS("Downloads/sporo_dataset.RDS")
sporo@meta.data$batches = ifelse(grepl("2_", rownames(sporo@meta.data)), "batch2", ifelse(grepl("3_", rownames(sporo@meta.data)), "batch3", ifelse(grepl("4_", rownames(sporo@meta.data)), "batch4", "batch1")))
saveRDS(sporo, "sporo_dataset_final.RDS")

#This is for label transfer
para.anchors <- FindTransferAnchors(reference = somules, query = sporo, dims = 1:30, features = intersect(VariableFeatures(somules), VariableFeatures(sporo)))

predictions <- TransferData(anchorset = para.anchors, refdata = Idents(somules), dims = 1:30)

sporo <- AddMetaData(sporo, metadata = predictions)
DimPlot(sporo, group.by = 'predicted.id', label = T, label.size = 2)



#This is for integration

para.list = list(somules, sporo)
anchors <- FindIntegrationAnchors(object.list = para.list, anchor.features = intersect(VariableFeatures(somules), VariableFeatures(sporo)), reduction = "rpca", k.anchor = 20)
para.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

library(ggplot2)
library(cowplot)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(para.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
para.integrated <- ScaleData(para.integrated, verbose = FALSE)
para.integrated <- RunPCA(para.integrated, npcs = 30, verbose = FALSE)
para.integrated <- RunUMAP(para.integrated, reduction = "pca", dims = 1:30)

DimPlot(para.integrated)
saveRDS(para.integrated, "para.integrated.RDS")
