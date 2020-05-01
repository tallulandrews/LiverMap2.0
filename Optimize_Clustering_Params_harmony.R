require("Seurat")

proj_name="harmony_integrated2_20livers"
#liver.integrated <- readRDS("harmony_integrated.rds")
#liver.integrated <- readRDS("integration_v2.rds")
liver.integrated <- readRDS("All_merged_universal_genes_harmony_integrated_v2.rds")

set.seed(7742)
# Dimensionality Reduction
#require("sctransform")
#liver.integrated <- ScaleData(liver.integrated);
#liver.integrated <- RunPCA(liver.integrated, features = VariableFeatures(object = liver.integrated))
npcs <- 15 # was 20 - 29/04/2020
require("ggplot2")
png(paste(proj_name,"pca_elbow.png", sep="_"), width=4, height=4, units="in", res=100)
a <- ElbowPlot(liver.integrated)+geom_vline(aes(xintercept=npcs), linetype="dotted", color="grey35", size=2)
plot(a)
dev.off()


# Cluster with many different parameters
res <- seq(from=0.3, to=2, by=0.3)
nkNN <- seq(from=40, to=80, by=10)

for(res_param in res) {
for(nkNN_param in nkNN){
liver.integrated <- FindNeighbors(liver.integrated, reduction="harmony", dims = 1:npcs, k.param=nkNN_param)
#liver.integrated <- FindClusters(liver.integrated, reduction="harmony", resolution = res_param, k.param=nkNN_param)
liver.integrated <- FindClusters(liver.integrated, reduction="harmony", resolution = res_param, group.singletons=FALSE, method="igraph") # for UHN Seurat installation

name <- paste("knn_",nkNN_param,"_res_", res_param, sep="");

liver.integrated@meta.data[[name]] <- liver.integrated@meta.data$seurat_clusters;

}}

saveRDS(liver.integrated, paste(proj_name, "harmony_plus_clusters.rds", sep="_"))

# Compare all these clusterings
require(igraph)
require(gplots)
clust_table <- liver.integrated@meta.data[, grepl("^knn_", colnames(liver.integrated@meta.data))]
clust_table <- as.matrix(apply(clust_table,2,as.numeric))
require("proxy")
clust_dists <- proxy::dist(clust_table, method=function(x,y){igraph::compare(x,y,method="vi")}, by_rows=FALSE)
#clust_similr1 <- proxy::simil(clust_table, method=function(x,y){igraph::compare(x,y,method="nmi")}, by_rows=FALSE)
#clust_similr2 <- proxy::simil(clust_table, method=function(x,y){igraph::compare(x,y,method="adjusted.rand")}, by_rows=FALSE)


# Find robust exemplar clustering(s)
require("apcluster")
require("gplots")
set.seed(18371)

res1 <- apcluster(-1*as.matrix(clust_dists), p=-2.5)
#res2 <- apcluster(as.matrix(clust_similr1), p=-2)
#res3 <- apcluster(as.matrix(clust_similr1), p=-2)

#valid_clusterings <- res1@exemplars[which(res1@exemplars %in% res2@exemplars & res1@exemplars %in% res3@exemplars)]

coarse_lvl <- names(res1@exemplars)[2]
fine_lvl <- names(res1@exemplars)[3]
#coarse_lvl <- "knn_70_res_0.3" # basic_integration_analysis.rds = knn_70_res_0.3
#fine_lvl <- "knn_70_res_0.9" # basic_integration_analysis.rds = knn_90_res_1.5

#manually select which exemplar to use
liver.integrated@meta.data$Coarse_clusters <- liver.integrated@meta.data[[coarse_lvl]] 
liver.integrated@meta.data$Fine_clusters <- liver.integrated@meta.data[[fine_lvl]]

apcluster::heatmap(res1, -1*as.matrix(clust_dists))

png(paste(proj_name,"compare_clusterings_heatmap.png",sep="_"), width=6, height=6, units="in", res=300)
lab <- matrix("", ncol=ncol(clust_table), nrow=ncol(clust_table))
lab[colnames(clust_table)==fine_lvl, colnames(clust_table)==fine_lvl] <- "F"
lab[colnames(clust_table)==coarse_lvl, colnames(clust_table)==coarse_lvl] <- "C"
heatmap.2(as.matrix(clust_dists), trace="none", distfun=function(x){return(as.dist(clust_dists))}, cellnote=lab)
dev.off()


# Visualize the Chosen clusterings
#nkNN <- 70

#liver.integrated <- RunTSNE(liver.integrated, reduction="harmony",  dims = 1:npcs)
#liver.integrated <- RunUMAP(liver.integrated, reduction="harmony",  dims = 1:npcs, parallel=FALSE, n.neighbour=nkNN)
png(paste(proj_name,"coarse_umap.png", sep="_"), width=6, height=6, units="in", res=150)
DimPlot(liver.integrated, reduction = "umap", group.by="Coarse_clusters")
dev.off()
png(paste(proj_name,"coarse_tsne.png",sep="_"), width=6, height=6, units="in", res=150)
DimPlot(liver.integrated, reduction = "tsne", group.by="Coarse_clusters")
dev.off()
png(paste(proj_name,"fine_umap.png", sep="_"), width=6, height=6, units="in", res=150)
DimPlot(liver.integrated, reduction = "umap", group.by="Fine_clusters")
dev.off()
png(paste(proj_name,"fine_tsne.png", sep="_"), width=6, height=6, units="in", res=150)
DimPlot(liver.integrated, reduction = "tsne", group.by="Fine_clusters")
dev.off()

saveRDS(liver.integrated, paste(proj_name, "harmony_plus_analysis.rds", sep="_"))
#png("Coarse_harmony_clusters_by_donor.png", width=6, height=6, units="in", res=150)
#barplot(table(liver.integrated@meta.data$orig.ident, liver.integrated@meta.data$Coarse_clusters), col=rainbow(20))
#dev.off()
#png("Fine_harmony_clusters_by_donor.png", width=6, height=6, units="in", res=150)
#barplot(table(liver.integrated@meta.data$orig.ident, liver.integrated@meta.data$Fine_clusters), col=rainbow(20))
#dev.off()

# Auto-annotation
#source("Setup_autoannotation.R")

#all_anno <- readRDS("All20_automatedannotation.rds");

#liver.integrated@meta.data$scmap_anno <- rep("unknown", ncol(liver.integrated));
#liver.integrated@meta.data$scmap_anno2 <- rep("unknown", ncol(liver.integrated));
#for (donor in unique(liver.integrated@meta.data$orig.ident)) {
#        cell_ids <- liver.integrated@meta.data$cell_barcode[liver.integrated@meta.data$donor == donor]
#        anno <- all_anno[[donor]];
#        anno <- anno[anno$cell_barcode %in% cell_ids,]
#        anno <- anno[match(cell_ids, anno$cell_barcode),]
#        liver.integrated@meta.data$scmap_anno[liver.integrated@meta.data$orig.ident == donor] <- as.character(anno$scmap_cluster_anno$lm1)
#        liver.integrated@meta.data$scmap_anno2[liver.integrated@meta.data$orig.ident == donor] <- as.character(anno$scmap_cell_anno)
#}

#png("AutoAnno_harmony_integrated_umap.png", width=8, height=6, units="in", res=150)
#DimPlot(liver.integrated, reduction = "umap", group.by="scmap_anno")
#dev.off()
#png("AutoAnno_harmony_integrated_tsne.png", width=8, height=6, units="in", res=150)
#DimPlot(liver.integrated, reduction = "tsne", group.by="scmap_anno")
#dev.off()
#png("AutoAnno2_harmony_integrated_umap.png", width=8, height=6, units="in", res=150)
#DimPlot(liver.integrated, reduction = "umap", group.by="scmap_anno2")
#dev.off()
#png("AutoAnno2_harmony_integrated_tsne.png", width=8, height=6, units="in", res=150)
#DimPlot(liver.integrated, reduction = "tsne", group.by="scmap_anno2")
#dev.off()

#saveRDS(liver.integrated, "integration_harmony_plus_analysis.rds")

# Would this be helpful?
# create a heatmap where: cell = average similarity of this clustering to all other clusterings
# distance between clusterings is measured using igraph::compare(c1, c2, method="vi")

