require("Seurat")

set.seed(3921)

# Which do we include in the integrated map?
dir <- 
seurfiles <- c("C37_Anno_SeurObj.rds", 
	"C39_Anno_SeurObj.rds",
	"C41_Anno_SeurObj.rds",
	"C41-NPC_Anno_SeurObj.rds",
	"C42_Anno_SeurObj.rds",
	"C42-NPC_Anno_SeurObj.rds",
	"C43_Anno_SeurObj.rds",
	"C43-NPC_Anno_SeurObj.rds",
	"C46_Anno_SeurObj.rds",
	"C48_Anno_SeurObj.rds",
	"C49_Anno_SeurObj.rds",
	"C50_Anno_SeurObj.rds",
	"C51_Anno_SeurObj.rds",
	"C52_Anno_SeurObj.rds",
	"C53_Anno_SeurObj.rds",
	"C54_Anno_SeurObj.rds",
	"C56_Anno_SeurObj.rds",
	"C58_Anno_SeurObj.rds",
	"C59_Anno_SeurObj.rds",
	"C61_Anno_SeurObj.rds",
	"C64_Anno_SeurObj.rds",
	"C68_Anno_SeurObj.rds",
	"C70_Anno_SeurObj.rds",
	"C72_Anno_SeurObj.rds");

samp_names <- unlist(lapply(strsplit(seurfiles, "_"), function(x){x[[1]]}))
obj_list <- list()
all_genes <- c();
for (i in 1:length(seurfiles)) {
	n <- samp_names[i];
	obj <- readRDS(seurfiles[i]);

	# fix some auto-annotation issues
	obj@meta.data$consistent_labs <- as.character(obj@meta.data$consistent_labs);
	anno1 <- obj@meta.data$scmap_cluster_anno
	anno2 <- obj@meta.data$scmap_cell_anno
	mac_groups <- c("Non-inflammatoryMacrophages", "InflamatoryMacrophages")
	g_anno <- as.character(obj@meta.data$general_labs)
	g_anno[g_anno == "ambiguous" & anno1 %in% mac_groups & anno2 %in% mac_groups] <- "Macrophage"
	obj@meta.data$general_labs <- g_anno
	null_vals <- is.na(obj@meta.data$consistent_labs);
	obj@meta.data$consistent_labs[null_vals] <- g_anno[null_vals]

	#Fix sample ID, and Donor ID
	obj@meta.data$sample <- obj@meta.data$orig.ident
	obj@meta.data$donor <- sapply(strsplit(as.character(obj@meta.data$sample), "_"), function(x){x[[1]]})
	# save sample specific clusters
	obj@meta.data$sample_specific_clusters <- obj@meta.data$seurat_clusters

	# get rid of factors
	metadata_classes <- sapply(1:ncol(obj@meta.data), function(i){class(obj@meta.data[,i])})
	for (j in which(metadata_classes == "factor")) {
		obj@meta.data[,j] <- as.character(obj@meta.data[,j]);
	}

	
	obj <- Seurat::NormalizeData(obj, verbose = FALSE, method="LogNormalize", scale.factor=10000) 
	obj <- Seurat::ScaleData(obj);
	obj@meta.data$cell_barcode <- colnames(obj);
	if (length(all_genes) == 0) {
		all_genes <- rownames(obj);
	} else {
		all_genes <- union(all_genes, rownames(obj));
	}
	obj@meta.data[[7]] <- obj@meta.data[[7]][[1]]
	obj@meta.data$donor <- rep(n, ncol(obj));
	obj_list[[n]] <- obj
}

all_genes <- sort(all_genes);
#all_genes <- all_genes[!grepl("^MT-", all_genes)]
#global_count_mat <- NULL
#global_meta_data <- NULL
#for (i in 1:length(obj_list)) {
#	obj_list[[i]] <- obj_list[[i]][match(all_genes, rownames(obj_list[[i]])),]
#}

# Merge Datasets
#### Merging does not merge individually scaled datasets!!

merged_obj <- NULL;
universal_genes <- c(-1)
all_genes <- c();
for (i in 1:length(obj_list)) {
        n <- samp_names[i];
	if (i == 1) {
		merged_obj <- obj_list[[i]]
		universal_genes <- as.character(rownames(obj_list[[i]]))
	} else {
		merged_obj <- merge(merged_obj, y=obj_list[[i]], add.cell.ids=c("", n), project="LiverMap")
		universal_genes <- intersect(universal_genes, as.character(rownames(obj_list[[i]])))
	}
}

merged_obj@misc$universal_genes <- universal_genes;
merged_obj@misc$creation_date <- date();
saveRDS(merged_obj, "All_genes_Merged_obj_v2.rds")

set.seed(9428)

merged_obj <- merged_obj[rownames(merged_obj) %in% universal_genes,]
#merged_obj <- Seurat::NormalizeData(merged_obj, verbose = FALSE)
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
merged_obj <- ScaleData(merged_obj, verbose = FALSE)
#merged_obj <- RunPCA(merged_obj, pc.genes = VariableGenes(merged_obj), npcs = 20, verbose = FALSE)
merged_obj <- RunPCA(merged_obj, pc.genes = VariableFeatures(merged_obj), npcs = 20, verbose = FALSE)
merged_obj <- RunTSNE(merged_obj, dims = 1:10, verbose = FALSE)
merged_obj <- RunUMAP(merged_obj, dims = 1:10, verbose = FALSE)

png("merged_not_integrated_tsne.png", width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="tsne", group.by="donor", pt.size=0.1)
dev.off();
png("merged_not_integrated_umap.png", width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="donor", pt.size=0.1)
dev.off();
png("merged_not_integrated_tsne_autoanno.png", width=12, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="tsne", group.by="consistent_labs", pt.size=0.1)
dev.off();
png("merged_not_integrated_umap_autoanno.png", width=12, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="consistent_labs", pt.size=0.1)
dev.off();

require("harmony")
set.seed(10131)
merged_obj <- RunHarmony(merged_obj, "donor", plot_convergence = TRUE)
#DimPlot(object = merged_obj, reduction = "harmony", pt.size = .1, group.by = "donor")

merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:20)
merged_obj <- RunTSNE(merged_obj, reduction = "harmony", dims = 1:20)

png("merged_harmony_integrated_umap.png", width=9, height =6, units="in", res=100)
DimPlot(merged_obj, reduction = "umap", group.by = "donor", pt.size = .1)
dev.off();
png("merged_harmony_integrated_tsne.png", width=9, height =6, units="in", res=100)
DimPlot(merged_obj, reduction = "tsne", group.by = "donor", pt.size = .1)
dev.off();
png("merged_harmony_integrated_umap_autoanno.png", width=12, height =6, units="in", res=100)
DimPlot(merged_obj, reduction = "umap", group.by = "consistent_labs", pt.size = .1)
dev.off();
png("merged_harmony_integrated_tsne_autoanno.png", width=12, height =6, units="in", res=100)
DimPlot(merged_obj, reduction = "tsne", group.by = "consistent_labs", pt.size = .1)
dev.off();
saveRDS(merged_obj, "All_merged_universal_genes_harmony_integrated_v2.rds");

# add harmony dimensions to integrated object?
