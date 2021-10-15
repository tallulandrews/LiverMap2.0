require("Seurat")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")


#dim_reduce <- mergedobj@reductions
# Coarse_cluster sets

cluster_sets <- list(
	Hepatocyte1 = c(0, 1, 4, 5, 9, 19), 
	Hepatocyte2=c(2, 17, 18), 
	NKT=c(3, 11, 12, 16), 
	Stellate=c(13), 
	Cholangiocyte=c(15), 
	Macrophage=c(7, 8), 
	Endo=c(6,10), 
	AntiB=c(14))
key_marker_sets <- list(
	Hepatocyte1 = c("CYP1A2", "CYP2E1", "CYP3A4", "GLUL", "CYP2A7", "FABP1", "HAL", "AGT", "SDS"), 
	Hepatocyte2=c("CYP1A2", "CYP2E1", "CYP3A4", "GLUL", "CYP2A7", "FABP1", "HAL", "AGT", "SDS"), 
	NKT=c("IGKC", "TRAC", "TRDC", "TRBC1", "TRBC2", "TRGC1", "TRGC2"), 
	Stellate=c(), 
	Cholangiocyte=c(), 
	Endo=c("TIMP3", "VWF", "ID1", "PECAM1", "DNASE1L3", "LIFR", "STAB1"), 
	Mac=c(), 
	AntiB=c("IGKC"))

set_i = 6;
mergedobj <- readRDS("Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno.rds")
#for (set_i in 1:length(cluster_sets)) {
proj_name <- names(cluster_sets)[set_i]

if (!file.exists(paste(proj_name, "harmony_Subcluster.rds", sep="_"))) {
	# Subset Object #

	subset <- mergedobj[,mergedobj@meta.data[,"Coarse_clusters"] %in% cluster_sets[[set_i]]]

	# require all genes to be detected in at least one cell in a majority of donors
	gene_filter <- group_rowmeans(subset@assays$RNA@counts, subset@meta.data$donor);
	gene_filter <- apply(gene_filter, 1, median)
	gene_filter <- gene_filter > 0;

	saveRDS(subset, paste(proj_name, "harmony_Subcluster_Allgenes.rds", sep="_")) # created later, to help annotate clusters.
	subset <- subset[gene_filter,]
	saveRDS(subset, paste(proj_name, "harmony_Subcluster.rds", sep="_"))
} else {
	subset <- readRDS(paste(proj_name, "harmony_Subcluster.rds", sep="_"));
}
# Save old harmony for later



proj_name <- names(cluster_sets)[set_i]
subset <- readRDS(paste(proj_name, "harmony_Subcluster.rds", sep="_"));

# Subcluster
subset <- FindVariableFeatures(subset, selection.method="vst", nfeatures=2500)
VariableFeatures(subset) <- unique(c(key_marker_sets[[set_i]], VariableFeatures(subset)));
subset <- RunPCA(subset, features=VariableFeatures(object=subset))



set.seed(2009)
require(harmony)
subset <- RunHarmony(subset, "sample", plot_convergence=FALSE)

npcs <- 10

res <- seq(from=0.3, to=1.5, by=0.2)
nkNN <- seq(from=30, to=60, by=10)

for(res_param in res) {
for(nkNN_param in nkNN){
	subset <- FindNeighbors(subset, reduction="harmony", dims = 1:npcs, k.param=nkNN_param)
	subset <- FindClusters(subset, reduction="harmony", resolution = res_param, k.param=nkNN_param)
	name <- paste("knn_",nkNN_param,"_res_", res_param, sep="");

	subset@meta.data[[name]] <- subset@meta.data$seurat_clusters;

}}

saveRDS(subset, paste(proj_name, "harmony_Subcluster.rds", sep="_"))

#}
exit();
# Compare all these clusterings
require(igraph)
require(gplots)
clust_table <- subset@meta.data[,which(colnames(subset@meta.data) == "Fine_Manual_Anno"):ncol(subset@meta.data)]
clust_table <- clust_table[, grepl("^knn_", colnames(clust_table))]
clust_table <- as.matrix(apply(clust_table,2,as.numeric))

require("proxy")
clust_dists <- proxy::dist(clust_table, method=function(x,y){igraph::compare(x,y,method="vi")}, by_rows=FALSE)

# Find robust exemplar clustering(s)
require("apcluster")
set.seed(18371)

res1 <- apcluster(-1*as.matrix(clust_dists), p=-1.5) #AntiB = -1.5, Cholangiocyte = -1.5, Endo = -1.5, Hepatocyte1=-1.5, Hepatocyte2=-2, NKT=-1.5, Stellate=-1.5, Macrophage = -1.5


# Cluster-Cell-type Purity
type_purity <- function(clusters, annotations) {
        tmp <- table(clusters, annotations)
        shan <- median(vegan::diversity(tmp, index="shannon", MARGIN=1))
        simp <- median(vegan::diversity(tmp, index="simpson", MARGIN=1))
        return(c(shan, simp))

}

score <- apply(clust_table, 2, type_purity, annotations=subset@meta.data$marker_labs)
#plot(score, xlab="clustering" ylab="type_purity")
score2 <- colMeans(score)

exemplars <- c();
for (c in res1@clusters) {
        clusterings <- names(c)
        scores <- score2[clusterings]
	best <- scores[which(scores == min(scores))]
        exemplars <- c(exemplars, best[1])
}

subset@misc$exemplars <- exemplars;

core_lvl <- names(exemplars)[1]
coarse_lvl <- names(exemplars)[2]
fine_lvl <- names(exemplars)[3]


#manually select which exemplar to use
subset@meta.data$Core_clusters <- subset@meta.data[[core_lvl]]
subset@meta.data$Coarse_clusters <- subset@meta.data[[coarse_lvl]]
subset@meta.data$Fine_clusters <- subset@meta.data[[fine_lvl]]

png(paste(proj_name,"Subcluster_compare_clusterings_heatmap.png",sep="_"), width=6, height=6, units="in", res=300)
lab <- matrix("", ncol=ncol(clust_table), nrow=ncol(clust_table))
lab[colnames(clust_table)==fine_lvl, colnames(clust_table)==fine_lvl] <- "3"
lab[colnames(clust_table)==coarse_lvl, colnames(clust_table)==coarse_lvl] <- "2"
lab[colnames(clust_table)==core_lvl, colnames(clust_table)==core_lvl] <- "1"
heatmap.2(as.matrix(clust_dists), trace="none", distfun=function(x){return(as.dist(clust_dists))}, cellnote=lab)
dev.off()


npcs <- 10
source("/cluster/home/tandrews/scripts/LiverMap2.0/Colour_Scheme.R")
# Visualize the Chosen clusterings
set.seed(34817)
#subset <- RunUMAP(subset, reduction="harmony",  dims = 1:npcs, parallel=FALSE)

png(paste(proj_name,"Subcluster_core_umap.png", sep="_"), width=6, height=6, units="in", res=150)
DimPlot(subset, reduction = "umap", group.by="Core_clusters", label=TRUE)
dev.off()

png(paste(proj_name,"Subcluster_coarse_umap.png", sep="_"), width=6, height=6, units="in", res=150)
DimPlot(subset, reduction = "umap", group.by="Coarse_clusters", label=TRUE)
dev.off()

#png(paste(proj_name,"Subcluster_fine_umap.png", sep="_"), width=6, height=6, units="in", res=150)
#DimPlot(subset, reduction = "umap", group.by="Fine_clusters")
#dev.off()

subset <- readRDS( paste(proj_name, "harmony_Subcluster.rds", sep="_"))

png(paste(proj_name,"Subcluster_sample_umap.png", sep="_"), width=6, height=6, units="in", res=150)
DimPlot(subset, reduction = "umap", group.by="sample")
dev.off()

png(paste(proj_name,"Subcluster_generalanno_umap.png", sep="_"), width=6, height=6, units="in", res=150)
Type_DimPlot(subset, type_col="Coarse_Manual_Anno", cluster_col="Coarse_Manual_Anno")
dev.off()

saveRDS(subset, paste(proj_name, "harmony_Subcluster.rds", sep="_"))


# Gather Stats

prev_anno <- table(subset@meta.data$Coarse_clusters, subset@meta.data$Coarse_Manual_Anno)
ncells <- table(subset@meta.data$Coarse_clusters)
by_sample <- table(subset@meta.data$Coarse_clusters, subset@meta.data$sample)
sample_diversity <- vegan::diversity(table(subset@meta.data$Coarse_clusters, subset@meta.data$sample), index="shannon", MARGIN=1)
CCPhase <-  table(subset@meta.data$Coarse_clusters, subset@meta.data$Phase)

out_table <- cbind(ncells, sample_diversity, CCPhase, prev_anno, by_sample)
write.table(out_table, paste(proj_name, "Subcluster_anno_info.txt", sep="_"), row.names=T, col.names=T)


