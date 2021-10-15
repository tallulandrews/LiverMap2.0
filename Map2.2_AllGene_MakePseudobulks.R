require("Seurat")
mergedobj <- readRDS("Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno.rds")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

obj_5pr <- mergedobj[,mergedobj@meta.data$assay_type == "5pr"]
obj_3pr <- mergedobj[,mergedobj@meta.data$assay_type == "3pr"]

dim(mergedobj)
cluster_5pr_pseudo <- get_pseudobulk(obj_5pr@assays$RNA@counts, obj_5pr@meta.data$Coarse_clusters, obj_5pr@meta.data$donor)
cluster_3pr_pseudo <- get_pseudobulk(obj_3pr@assays$RNA@counts, obj_3pr@meta.data$Coarse_clusters, obj_3pr@meta.data$donor)


