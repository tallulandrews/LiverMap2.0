all <- readRDS("../All_genes_Merged_obj_v2.rds")
#Get metadata from analyized integrated data
b <- readRDS("../harmony_integrated2_20livers_harmony_plus_analysis.rds")
mdata <- b@meta.data
rm(b)
# Check cell-IDs match
identical(mdata$cell_ID, all@meta.data$cell_ID)
# Remove temporary clustering labels
mdata2 <- mdata[,-grepl("^knn_", colnames(mdata))]
head(mdata2)
grepl("^knn", colnames(mdata))
mdata2 <- mdata[,!grepl("^knn_", colnames(mdata))]
head(mdata2)
mdata2 <- mdata2[,!grepl("^RNA_snn", colnames(mdata2))]
head(mdata2)
# Attach to file with all genes
all@meta.data <- mdata2
# Subset by type
a <- all[,all@meta.data$Fine_clusters %in% c(0,14,27)]
saveRDS(a,"All_genes_Subcluster_Hep-Erythroid.rds")
a <- all[,all@meta.data$Fine_clusters %in% c(1,2,3,4,6,7,17,18,24,26)]
saveRDS(a,"All_genes_Subcluster_Hep-Traj.rds")
a <- all[,all@meta.data$Fine_clusters %in% c(5,10,20)]
saveRDS(a,"All_genes_Subcluster_NKT.rds")
a <- all[,all@meta.data$Fine_clusters %in% c(8,12,15,19,21)]
saveRDS(a,"All_genes_Subcluster_Endo.rds")
a <- all[,all@meta.data$Fine_clusters %in% c(9,11,13,23)]
saveRDS(a,"All_genes_Subcluster_Mac.rds")
a <- all[,all@meta.data$Fine_clusters %in% c(9,11,13,21)]
saveRDS(a,"All_genes_Subcluster_Mac.rds")
a <- all[,all@meta.data$Fine_clusters %in% c(16,23)]
saveRDS(a,"All_genes_Subcluster_Bcell.rds")
savehistory("Setup_Subcluster.Rhistory")
