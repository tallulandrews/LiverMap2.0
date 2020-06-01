require("Seurat")
source("My_R_Scripts.R")

proj_name = "Tcell"

liver.integrated <- readRDS("integration_harmony_plus_analysis.rds")

cluster_labs <- liver.integrated@meta.data$Fine_clusters
names(cluster_labs) <- paste(liver.integrated@meta.data$donor, liver.integrated@meta.data$cell_barcode, sep="_")
harmony_redim <- liver.integrated@reductions$harmony
meta.data <- liver.integrated@meta.data
meta.data <- meta.data[,!grepl("^knn", colnames(meta.data))]


rm(liver.integrated)
liver.integrated <- readRDS("all_gene_merged_obj.rds")

harmony_redim@cell.embeddings <- harmony_redim@cell.embeddings[match(colnames(liver.integrated), rownames(harmony_redim@cell.embeddings)),]
liver.integrated@reductions$harmony <- harmony_redim;


meta.data <- meta.data[match(colnames(liver.integrated), rownames(meta.data)),]
liver.integrated@meta.data <- meta.data


liver.integrated <- liver.integrated[,liver.integrated@meta.data$Fine_clusters %in% c(8,22)];

# require all genes to be detected in at least one cell in a majority of donors
gene_filter <- group_rowmeans(liver.integrated@assays$RNA@counts, liver.integrated@meta.data$donor);
gene_filter <- apply(gene_filter, 1, median)
gene_filter <- gene_filter > 0;

liver.integrated <- liver.integrated[gene_filter,]

set.seed(7742)

# Dimensionality Reduction
### Probably need to re-run harmony on just these cell for best results....####
liver.integrated <- NormalizeData(liver.integrated)
liver.integrated <- ScaleData(liver.integrated);
liver.integrated <- FindVariableFeatures(liver.integrated, selection.method = "vst", nfeatures = 2500)
liver.integrated <- RunPCA(liver.integrated, features = VariableFeatures(object = liver.integrated))
ElbowPlot(liver.integrated)



require("harmony")
set.seed(10131)
liver.integrated <- RunHarmony(liver.integrated, "donor", plot_convergence = TRUE)



npcs <- 15

# Cluster with many different parameters
res <- seq(from=0.1, to=1.5, by=0.2)
nkNN <- seq(from=30, to=60, by=10)

for(res_param in res) {
for(nkNN_param in nkNN){
liver.integrated <- FindNeighbors(liver.integrated, reduction="harmony", dims = 1:npcs, k.param=nkNN_param)
liver.integrated <- FindClusters(liver.integrated, reduction="harmony", resolution = res_param, k.param=nkNN_param)
name <- paste("knn_",nkNN_param,"_res_", res_param, sep="");

liver.integrated@meta.data[[name]] <- liver.integrated@meta.data$seurat_clusters;

}}

# Compare all these clusterings
require(igraph)
require(gplots)
clust_table <- liver.integrated@meta.data[, grepl("^knn_", names(liver.integrated@meta.data))]
clust_table <- as.matrix(apply(clust_table,2,as.numeric))
require("proxy")
clust_dists <- proxy::dist(clust_table, method=function(x,y){igraph::compare(x,y,method="vi")}, by_rows=FALSE)
clust_similr1 <- proxy::simil(clust_table, method=function(x,y){igraph::compare(x,y,method="nmi")}, by_rows=FALSE)
clust_similr2 <- proxy::simil(clust_table, method=function(x,y){igraph::compare(x,y,method="adjusted.rand")}, by_rows=FALSE)


# Find robust exemplar clustering(s)
require("apcluster")
require("gplots")
set.seed(18371)

res1 <- apcluster(-1*as.matrix(clust_dists), p=-5)
#res2 <- apcluster(as.matrix(clust_similr1), p=-2)
#res3 <- apcluster(as.matrix(clust_similr1), p=-2)

#valid_clusterings <- res1@exemplars[which(res1@exemplars %in% res2@exemplars & res1@exemplars %in% res3@exemplars)]

coarse_lvl <- "knn_50_res_0.9" # <- Update
fine_lvl <- "knn_70_res_1.1" # <- Update

#manually select which exemplar to use
liver.integrated@meta.data$Coarse_clusters <- liver.integrated@meta.data[[coarse_lvl]] 
liver.integrated@meta.data$Fine_clusters <- liver.integrated@meta.data[[fine_lvl]]

apcluster::heatmap(res1, -1*as.matrix(clust_dists))

png(paste(proj_name, "compare_clusterings_heatmap.png", sep="_"), width=6, height=6, units="in", res=300)
lab <- matrix("", ncol=ncol(clust_table), nrow=ncol(clust_table))
lab[colnames(clust_table)==fine_lvl, colnames(clust_table)==fine_lvl] <- "F"
lab[colnames(clust_table)==coarse_lvl, colnames(clust_table)==coarse_lvl] <- "C"
heatmap.2(as.matrix(clust_dists), trace="none", distfun=function(x){return(as.dist(clust_dists))}, cellnote=lab)
dev.off()


# Visualize the Chosen clusterings
nkNN <- 50

set.seed(07219)

liver.integrated <- RunTSNE(liver.integrated, reduction="harmony",  dims = 1:npcs)
liver.integrated <- RunUMAP(liver.integrated, reduction="harmony",  dims = 1:npcs, parallel=FALSE, n.neighbour=nkNN)
png(paste(proj_name, "coarse_umap.png", sep="_"), width=6, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "umap", group.by="Coarse_clusters")
dev.off()
png(paste(proj_name, "coarse_tsne.png", sep="_"), width=6, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "tsne", group.by="Coarse_clusters")
dev.off()
png(paste(proj_name, "coarse_umap_donor.png", sep="_"), width=6, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "umap", group.by="orig.ident")
dev.off()
DimPlot(liver.integrated, reduction = "tsne", group.by="orig.ident")

barplot(table(liver.integrated@meta.data$orig.ident, liver.integrated@meta.data$Coarse_clusters), col=rainbow(20))

tab <- table(liver.integrated@meta.data$orig.ident, liver.integrated@meta.data$Coarse_clusters)

DimPlot(liver.integrated, reduction = "umap", group.by="scmap_anno")
DimPlot(liver.integrated, reduction = "tsne", group.by="scmap_anno")
png(paste(proj_name, "coarse_umap_autoanno.png", sep="_"), width=6, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "umap", group.by="scmap_anno2")
dev.off()
DimPlot(liver.integrated, reduction = "tsne", group.by="scmap_anno2")

saveRDS(liver.integrated, paste(proj_name, "_analysis.rds", sep=""))

