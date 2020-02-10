require("Seurat")

set.seed(3921)

seurfiles <- c("C37_Anno_SeurObj.rds", "C46_Anno_SeurObj.rds", "C52_Anno_SeurObj.rds", "C59_Anno_SeurObj.rds", "C39_Anno_SeurObj.rds", "C48_Anno_SeurObj.rds", "C53_Anno_SeurObj.rds", "C61_Anno_SeurObj.rds", "C41_Anno_SeurObj.rds", "C49_Anno_SeurObj.rds", "C54_Anno_SeurObj.rds", "C62_Anno_SeurObj.rds", "C42_Anno_SeurObj.rds", "C50_Anno_SeurObj.rds", "C56_Anno_SeurObj.rds", "C63_Anno_SeurObj.rds", "C43_Anno_SeurObj.rds", "C51_Anno_SeurObj.rds", "C58_Anno_SeurObj.rds", "C64_Anno_SeurObj.rds");

samp_names <- unlist(lapply(strsplit(seurfiles, "_"), function(x){x[[1]]}))
obj_list <- list()
common_genes <- c();
for (i in 1:length(seurfiles)) {
	n <- samp_names[i];
	obj <- readRDS(seurfiles[i]);
	obj@meta.data$cell_barcode <- colnames(obj);
	if (length(common_genes) == 0) {
		common_genes <- rownames(obj);
	} else {
		common_genes <- intersect(common_genes, rownames(obj));
	}
	obj@meta.data$donor <- rep(n, ncol(obj));
	obj_list[[n]] <- obj
}

common_genes <- sort(common_genes);
common_genes <- common_genes[!grepl("^MT-", common_genes)]
for (i in 1:length(obj_list)) {
	obj_list[[i]] <- obj_list[[i]][match(common_genes, rownames(obj_list[[i]])),]
}

#sctransform workflow - not available is installed version of seurat
#liver.features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 5000)
#obj_list <- PrepSCTIntegration(object.list=obj_list, anchor.features = liver.features, verbose=F)

#anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT", anchor.features = liver.features, verbose = FALSE)

#liver.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

#saveRDS(liver.integrated, "sctintegration.rds")

#alternat workflow
options(future.globals.maxSize = 10000 * 1024^2)
set.seed(10129)

liver.features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 1000)
anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = liver.features)

liver.integrated <- IntegrateData(anchorset = anchors)
saveRDS(liver.integrated, "integration_v2.rds")

# original workflow

#anchors <- FindIntegrationAnchors(object.list=obj_list, dims=1:20)

#liver.integrated <- IntegrateData(anchorset = anchors, dims=1:20)

#saveRDS(liver.integrated, "basicintegration.rds")


# Harmony Workflow

merged_obj <- NULL;
common_genes <- c();
for (i in 1:length(obj_list)) {
        n <- samp_names[i];
	if (i == 1) {
		merged_obj <- obj_list[[i]]
	} else {
		merged_obj <- merge(merged_obj, y=obj_list[[i]], add.cell.ids=c("", n), project="LiverMap")
	}
}

merged_obj <- Seurat::NormalizeData(merged_obj, verbose = FALSE) 
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
merged_obj <- ScaleData(merged_obj, verbose = FALSE)
merged_obj <- RunPCA(merged_obj, pc.genes = VariableGenes(merged_obj), npcs = 20, verbose = FALSE)
merged_obj <- RunTSNE(merged_obj, dims = 1:10, verbose = FALSE)

png("merged_not_integrated.png", width=6, height =6, units="in", res=50)
DimPlot(merged_obj, reduction="tsne", group.by="donor", pt.size=0.1)
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
