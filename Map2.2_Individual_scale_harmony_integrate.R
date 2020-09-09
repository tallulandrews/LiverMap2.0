require("Seurat")
source("~/scripts/LiverMap2.0/Colour_Scheme.R")

set.seed(3921)

# Which do we include in the integrated map?
dir <- "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned"
seurfiles <- c(
	"C37_EmptyOnly.rds",
	"C39_EmptyOnly.rds",
	"C41_EmptyOnly.rds", 
	"C41-NPC_EmptyOnly.rds",
	"C42_EmptyOnly.rds",
	"C42-NPC_EmptyOnly.rds",
	"C43_EmptyOnly.rds",
	"C43-NPC_EmptyOnly.rds",
	"C46_RESEQ_EmptyOnly.rds",
	"C48_EmptyOnly.rds",
	"C49_EmptyOnly.rds",
	"C50_EmptyOnly.rds",
	"C51_EmptyOnly.rds",
	"C51_Flush_EmptyOnly.rds",
	"C52_EmptyOnly.rds",
	"C53_RESEQ_EmptyOnly.rds",
	"C54_EmptyOnly.rds",
	"C56_RESEQ_EmptyOnly.rds",
	"C58_5pr_EmptyOnly.rds",
	"C58_RESEQ_EmptyOnly.rds",
	"C59_5pr_EmptyOnly.rds",
	"C59_EmptyOnly.rds",
	"C61_5pr_EmptyOnly.rds",
	"C61_RESEQ_EmptyOnly.rds",
	"C63_5pr_reseq_EmptyOnly.rds",
	"C63_reseq_EmptyOnly.rds",
	"C64_5pr_EmptyOnly.rds",
	"C64_RESEQ_EmptyOnly.rds",
	"C66_RESEQ_EmptyOnly.rds",
	"C68_RESEQ_EmptyOnly.rds",
	"C69_EmptyOnly.rds",
	"C70_5pr_reseq_EmptyOnly.rds",
	"C70_RESEQ_EmptyOnly.rds",
	"C72_RESEQ_EmptyOnly.rds"
	);

samp_names <- unlist(lapply(strsplit(seurfiles, "_"), function(x){x <- x[c(-length(x))]; return(paste(x, collapse="_"))}))

samp_type <- unlist(lapply(strsplit(seurfiles, "[_-]"), function(x){x[2]}))
samp_type[grepl(".rds", samp_type)] <- ""

obj_list <- list()
for (i in 1:length(seurfiles)) {
	n <- samp_names[i];
	obj <- readRDS(paste(dir, seurfiles[i], sep="/"));

	# Filter
	numi <- Matrix::colSums(obj@assays$RNA@counts)
	ngene <- Matrix::colSums(obj@assays$RNA@counts > 0)
	obj <- obj[,numi > 875 & ngene > 500]
	obj <- NormalizeData(obj, scale.factor=1500)

	#Fix sample ID, and Donor ID
	obj@meta.data$sample <- obj@meta.data$orig.ident
	
	obj@meta.data$donor <- sapply(strsplit(as.character(obj@meta.data$sample), "_"), function(x){x[[1]]})
	obj@meta.data$assay_type <- rep(samp_type[i], nrow(obj@meta.data));

	# save sample specific clusters
	obj@meta.data$sample_specific_clusters <- paste(n, obj@meta.data$seurat_clusters, sep="_")

	# get rid of factors
	metadata_classes <- sapply(1:ncol(obj@meta.data), function(i){class(obj@meta.data[,i])})
	for (j in which(metadata_classes == "factor")) {
		obj@meta.data[,j] <- as.character(obj@meta.data[,j]);
	}

	
	obj@meta.data$cell_barcode <- colnames(obj);
	obj@meta.data$sample <- rep(n, ncol(obj));
	obj@meta.data$cell_ID <- paste(obj@meta.data$sample, obj@meta.data$cell_barcode, sep="_")
	obj_list[[n]] <- obj
}

# Merge Datasets
#### Merging does not merge individually scaled datasets!!

# Find common HVGs and detected genes
merged_obj <- NULL;
universal_genes <- c(-1)
hvgs <- c();
for (i in 1:length(obj_list)) {
        n <- samp_names[i];
	if (i == 1) {
		merged_obj <- obj_list[[i]]
		universal_genes <- as.character(rownames(obj_list[[i]]))
		hvgs <- VariableFeatures(obj_list[[i]]);
	} else {
		merged_obj <- merge(merged_obj, y=obj_list[[i]], add.cell.ids=c("", n), project="LiverMap")
		universal_genes <- intersect(universal_genes, as.character(rownames(obj_list[[i]])))
		hvgs <- c(hvgs, VariableFeatures(obj_list[[i]]));
	}
}

fix_names <- paste(merged_obj@meta.data$orig.ident, merged_obj@meta.data$cell_barcode, sep="_")
merged_obj <- RenameCells(merged_obj, new.names=fix_names)

# Keep HVGs seen in at least 2 datasets
hvgs <- unique(hvgs[duplicated(hvgs)])
hvgs <- hvgs[!grepl("^MT-", hvgs)]
merged_obj@misc$repeated_hvgs <- hvgs;
hvgs <- hvgs[ hvgs %in% universal_genes] # might want to get rid of this.
merged_obj@misc$repeated_universal_hvgs <- hvgs;


# Scale within datasets
all_Scaled <- c();
scaled_cell_ids <- c()
for (i in 1:length(obj_list)) {
      n <- samp_names[i];
	obj <- obj_list[[i]]
	obj <- Seurat::ScaleData(obj, features=rownames(obj));
	
	scaled <- obj@assays$RNA@scale.data;
	scaled <- scaled[match(rownames(merged_obj), rownames(scaled)),]
	scaled[is.na(scaled)] <- 0;
	rownames(scaled) <- rownames(merged_obj);
	scaled_cell_ids <- c(scaled_cell_ids, obj_list[[i]]@meta.data$cell_ID);
	if (i == 1) {
		all_Scaled <- scaled;
	} else {
		# This might cause problems with genes not present in all datasets
		scaled <- scaled[match(rownames(all_Scaled), rownames(scaled)),]
		all_Scaled <- cbind(all_Scaled, scaled);
	}
}

colnames(all_Scaled) <- scaled_cell_ids;
merged_obj@assays$RNA@scale.data <- all_Scaled

merged_obj@misc$universal_genes <- universal_genes;
merged_obj@misc$creation_date <- date();
VariableFeatures(merged_obj) <- hvgs;
merged_obj@meta.data$seurat_clusters <- paste(merged_obj@meta.data$orig.ident,
		as.character(merged_obj@meta.data$seurat_clusters), sep="_")
saveRDS(merged_obj, "Merged_EmptyOnly_obj_Map2.2.rds")

set.seed(9428)

merged_obj <- merged_obj[rownames(merged_obj) %in% universal_genes,]
merged_obj@assays$RNA@scale.data <- all_Scaled
merged_obj <- RunPCA(merged_obj, pc.genes = hvgs, 
			npcs = 20, verbose = FALSE)
merged_obj <- RunTSNE(merged_obj, dims = 1:10, verbose = FALSE)
merged_obj <- RunUMAP(merged_obj, dims = 1:10, verbose = FALSE)

png("Map2.2_scaled_merged_not_integrated_tsne.png", width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="tsne", group.by="sample", pt.size=0.1)
dev.off();
png("Map2.2_scaled_merged_not_integrated_umap.png", width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="sample", pt.size=0.1)
dev.off();
png("Map2.2_scaled_merged_not_integrated_tsne_autoanno.png", width=12, height =6, units="in", res=300)
Type_DimPlot(merged_obj, reduction="tsne", type_col="marker_labs", cluster_col="marker_labs")
dev.off();
png("Map2.2_scaled_merged_not_integrated_umap_autoanno.png", width=12, height =6, units="in", res=300)
Type_DimPlot(merged_obj, reduction="umap", type_col="marker_labs", cluster_col="marker_labs")
dev.off();

#rescale across datasets
obj <- Seurat::ScaleData(merged_obj, features=hvgs);

obj <- RunPCA(obj, pc.genes = hvgs, 
			npcs = 20, verbose = FALSE)
obj <- RunTSNE(obj, dims = 1:10, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:10, verbose = FALSE)

png("Map2.2_rescaled_merged_not_integrated_tsne.png", width=9, height =6, units="in", res=300)
DimPlot(obj, reduction="tsne", group.by="sample", pt.size=0.1)
dev.off();
png("Map2.2_rescaled_merged_not_integrated_umap.png", width=9, height =6, units="in", res=300)
DimPlot(obj, reduction="umap", group.by="sample", pt.size=0.1)
dev.off();
png("Map2.2_rescaled_merged_not_integrated_tsne_autoanno.png", width=12, height =6, units="in", res=300)
Type_DimPlot(obj, reduction="tsne", type_col="marker_labs", cluster_col="marker_labs")
dev.off();
png("Map2.2_rescaled_merged_not_integrated_umap_autoanno.png", width=12, height =6, units="in", res=300)
Type_DimPlot(obj, reduction="umap", type_col="marker_labs", cluster_col="marker_labs")
dev.off();


merged_obj@meta.data$assay_type <- rep("single_cell", nrow(merged_obj@meta.data))
nuclei <- grepl("ST", merged_obj@meta.data$sample)
merged_obj@meta.data$assay_type[nuclei] <- "single_nuc"

require("harmony")
set.seed(10131)
merged_obj <- RunHarmony(merged_obj, c("sample", "assay_type", "donor"), plot_convergence = TRUE)

merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:20)
merged_obj <- RunTSNE(merged_obj, reduction = "harmony", dims = 1:20)

png("Map2.2_merged_harmony_integrated_umap.png", width=9, height =6, units="in", res=100)
DimPlot(merged_obj, reduction = "umap", group.by = "sample", pt.size = .1)
dev.off();
png("Map2.2_merged_harmony_integrated_tsne.png", width=9, height =6, units="in", res=100)
DimPlot(merged_obj, reduction = "tsne", group.by = "sample", pt.size = .1)
dev.off();
png("Map2.2_merged_harmony_integrated_umap_autoanno.png", width=12, height =6, units="in", res=100)
Type_DimPlot(merged_obj, reduction="umap", type_col="marker_labs", cluster_col="marker_labs")
dev.off();
png("Map2.2_merged_harmony_integrated_tsne_autoanno.png", width=12, height =6, units="in", res=100)
Type_DimPlot(merged_obj, reduction="tsne", type_col="marker_labs", cluster_col="marker_labs")
dev.off();
saveRDS(merged_obj, "Map2.2_EmptyOnly_merged_integrated.rds");

# add harmony dimensions to integrated object?
