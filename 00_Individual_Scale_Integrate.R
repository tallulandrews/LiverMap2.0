require("Seurat")

set.seed(3921)

# Datasets to Integrate:
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

# Extract sample name from the file name.
samp_names <- unlist(lapply(strsplit(seurfiles, "_"), function(x){x <- x[c(-length(x))]; return(paste(x, collapse="_"))}))

# Extract assay type from file name.
samp_type <- unlist(lapply(strsplit(seurfiles, "[_-]"), function(x){x[2]}))
samp_type[grepl(".rds", samp_type)] <- ""

# Get each object, filter it and Normalize it to the same total UMI/cell
obj_list <- list()
for (i in 1:length(seurfiles)) {
	n <- samp_names[i];
	obj <- readRDS(paste(dir, seurfiles[i], sep="/"));

	# Filter
	numi <- Matrix::colSums(obj@assays$RNA@counts)
	ngene <- Matrix::colSums(obj@assays$RNA@counts > 0)
	obj <- obj[,numi > 875 & ngene > 500]
	obj <- NormalizeData(obj, scale.factor=1500)

	#Add some metadata for this sample
	obj@meta.data$sample <- obj@meta.data$orig.ident
	obj@meta.data$donor <- sapply(strsplit(as.character(obj@meta.data$sample), "_"), function(x){x[[1]]})
	obj@meta.data$assay_type <- rep(samp_type[i], nrow(obj@meta.data));
	obj@meta.data$cell_barcode <- colnames(obj);
	obj@meta.data$cell_ID <- paste(obj@meta.data$sample, obj@meta.data$cell_barcode, sep="_")

	# save sample specific clusters
	obj@meta.data$sample_specific_clusters <- paste(n, obj@meta.data$seurat_clusters, sep="_")

	# get rid of factors in metadata to prevent errors
	metadata_classes <- sapply(1:ncol(obj@meta.data), function(i){class(obj@meta.data[,i])})
	for (j in which(metadata_classes == "factor")) {
		obj@meta.data[,j] <- as.character(obj@meta.data[,j]);
	}

	# Add it to the list to merge	
	obj_list[[n]] <- obj
}

# Merge Datasets
#### Merging does not merge individually scaled datasets!!

# Merge Datasets
merged_obj <- NULL;
universal_genes <- c(-1) # detected in all datasets
hvgs <- c(); #all highly variable genes across all datasets
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

# Fix the cell names to include original sample ID
fix_names <- paste(merged_obj@meta.data$orig.ident, merged_obj@meta.data$cell_barcode, sep="_")
merged_obj <- RenameCells(merged_obj, new.names=fix_names)

# Keep only HVGs seen in multiple datasets
hvgs <- unique(hvgs[duplicated(hvgs)]) # Keep HVGs seen in at least 2 datasets
hvgs <- hvgs[!grepl("^MT-", hvgs)] # exclude mitochondrial genes from HVGs
merged_obj@misc$repeated_hvgs <- hvgs;
hvgs <- hvgs[ hvgs %in% universal_genes] # might want to get rid of this.
merged_obj@misc$repeated_universal_hvgs <- hvgs;


# Scale within datasets
all_Scaled <- c();
scaled_cell_ids <- c()
for (i in 1:length(obj_list)) {
      n <- samp_names[i];
	obj <- obj_list[[i]]
	obj <- Seurat::ScaleData(obj, features=rownames(obj)); # Scale this sample
	
	scaled <- obj@assays$RNA@scale.data;
	scaled <- scaled[match(rownames(merged_obj), rownames(scaled)),] # match the gene order to the merged object.
	scaled[is.na(scaled)] <- 0; # missing genes from this sample are set to all 0s
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
VariableFeatures(merged_obj) <- hvgs; #Keep HVGs that were variable in at least 2 datasets, and detected in all datasets.

# Same the merged object
saveRDS(merged_obj, "Merged_EmptyOnly_obj_Map2.2.rds")
