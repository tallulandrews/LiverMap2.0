require("Seurat")

set.seed(3921)

seurfiles <- c("C37_Anno_SeurObj2.rds", "C46_Anno_SeurObj2.rds", "C52_Anno_SeurObj2.rds", "C59_Anno_SeurObj2.rds", "C39_Anno_SeurObj2.rds", "C48_Anno_SeurObj2.rds", "C53_Anno_SeurObj2.rds", "C61_Anno_SeurObj2.rds", "C41_Anno_SeurObj2.rds", "C49_Anno_SeurObj2.rds", "C54_Anno_SeurObj2.rds", "C62_Anno_SeurObj2.rds", "C42_Anno_SeurObj2.rds", "C50_Anno_SeurObj2.rds", "C56_Anno_SeurObj2.rds", "C63_Anno_SeurObj2.rds", "C43_Anno_SeurObj2.rds", "C51_Anno_SeurObj2.rds", "C58_Anno_SeurObj2.rds", "C64_Anno_SeurObj2.rds");

samp_names <- unlist(lapply(strsplit(seurfiles, "_"), function(x){x[[1]]}))
obj_list <- list()
all_genes <- c();
for (i in 1:length(seurfiles)) {
	n <- samp_names[i];
	obj <- readRDS(seurfiles[i]);
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

merged_obj <- NULL;
all_genes <- c();
for (i in 1:length(obj_list)) {
        n <- samp_names[i];
	if (i == 1) {
		merged_obj <- obj_list[[i]]
	} else {
		merged_obj <- merge(merged_obj, y=obj_list[[i]], add.cell.ids=c("", n), project="LiverMap")
	}
}
saveRDS(merged_obj, "all_gene_merged_obj.rds")

exit;


merged_obj <- Seurat::NormalizeData(merged_obj, verbose = FALSE) 
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
merged_obj <- ScaleData(merged_obj, verbose = FALSE)
merged_obj <- RunPCA(merged_obj, pc.genes = VariableGenes(merged_obj), npcs = 20, verbose = FALSE)
merged_obj <- RunTSNE(merged_obj, dims = 1:10, verbose = FALSE)
merged_obj <- RunUMAP(merged_obj, dims = 1:10, verbose = FALSE)

png("merged_not_integrated_tsne.png", width=6, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="tsne", group.by="donor", pt.size=0.1)
dev.off();
png("merged_not_integrated_umap.png", width=6, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="donor", pt.size=0.1)
dev.off();

require("harmony")
set.seed(10131)
merged_obj <- RunHarmony(merged_obj, "donor", plot_convergence = TRUE)
DimPlot(object = merged_obj, reduction = "harmony", pt.size = .1, group.by = "donor")

merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:20)
merged_obj <- RunTSNE(merged_obj, reduction = "harmony", dims = 1:20)

png("merged_harmony_integrated_umap.png", width=6, height =6, units="in", res=50)
DimPlot(merged_obj, reduction = "umap", group.by = "donor", pt.size = .1)
dev.off();
png("merged_harmony_integrated.png", width=6, height =6, units="in", res=50)
DimPlot(merged_obj, reduction = "tsne", group.by = "donor", pt.size = .1)
dev.off();
saveRDS(merged_obj, "harmony_integrated.rds");

# add harmony dimensions to integrated object?
