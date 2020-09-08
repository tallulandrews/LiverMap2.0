
TCR_files <- list("C63_5pr" = "/cluster/projects/macparland/TA/LiverMap2.0/RawData/McGilvray_Sonya__C-63_Enriched_TCR/all_contig_annotations.csv",
		  "C70_5pr_reseq" = "/cluster/projects/macparland/TA/LiverMap2.0/RawData/MacParland_Sonya__C70_Caudate_TCR/all_contig_annotations.csv",
		  "C61_5pr" = "/cluster/projects/macparland/TA/LiverMap2.0/RawData/C61_5pr_TCR_20181121/all_contig_annotations.csv")

BCR_files <- list("C70_5pr_reseq" = "/cluster/projects/macparland/TA/LiverMap2.0/RawData/MacParland_Sonya__C70_Caudate_BCR/all_contig_annotations.csv",
		  "C63_5pr" = "/cluster/projects/macparland/TA/LiverMap2.0/RawData/McGilvray_Sonya__C-63_Enriched_BCR/all_contig_annotations.csv")


RNA_files <- list("C70_5pr_reseq"="/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C70_5pr_reseq_SoupX.rds", 
		   "C63_5pr"="/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C63_5pr_SoupX.rds",
		   "C61_5pr"="/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C61_5pr_SoupX.rds"
		)

MAIT <- read.delim("/cluster/home/tandrews/scripts/LiverMap2.0/TCR_Invariant.csv", sep=",", header=TRUE)

clone_type_diversity <- function(clonetype_vec) {
	#Simpson Index
	freq <- as.vector(unlist(table(clonetype_vec)))/length(clonetype_vec);
	simpson_raw <- sum(freq^2)
	# Adjust to account for different numbers of cells in the group
	theoretical_max <- length(clonetype_vec)*(1/length(clonetype_vec))^2
	return(simpson_raw/theoretical_max);
}

fix_cell_id <- function(barcodes, name) {
	barcodes <- sub("-1","", barcodes);
	barcodes <- paste(name, barcodes, sep="_");
	return(barcodes)
}

get_clone_types <- function(clone_type_mat, seur_obj) {
		# Match clone types to RNAseq
		clone_types <- unique(clone_type_mat[,c("cell", "raw_clonotype_id")])
		clone_types <- clone_types[match(seur_obj@meta.data$cell_ID, clone_types[,1]),]
		return(clone_types)
}

get_receptors <- function(clone_type_mat, seur_obj) {
		clone_type_mat$full_receptor <- paste(clone_type_mat$v_gene, clone_type_mat$d_gene, clone_type_mat$j_gene, clone_type_mat$c_gene, sep="_")
		cell_to_clone_type_mat <- aggregate(clone_type_mat$full_receptor, by=list(cell=clone_type_mat$cell), paste)
		cell_to_clone_type_mat[,2] <- unlist(lapply(cell_to_clone_type_mat[,2], function(x){paste(x, collapse=";")}))
		cell_to_clone_type_mat <- cell_to_clone_type_mat[match(seur_obj@meta.data$cell_ID, cell_to_clone_type_mat$cell),]
		return(cell_to_clone_type_mat);
}

obj_list <- list()

for (sample in names(RNA_files)) {
	seur_obj <- readRDS(RNA_files[[sample]])

	if (sample %in% names(TCR_files)) {
		# Read in TCR files
		TCRs <- read.delim(TCR_files[[sample]], sep=",", header=T)
		# Fix cell-ids
		TCRs$cell <- fix_cell_id(TCRs$barcode, sample)
		# Add clone-type metadata
		seur_obj@meta.data$TCR_clonetype <- get_clone_types(TCRs, seur_obj)[,2]

		# Look for Invariant
		TCRs$V_J_C <- paste(TCRs$v_gene, TCRs$j_gene, TCRs$c_gene, sep = "_")
		MAIT$V_J_C <- paste(MAIT$V, MAIT$J, MAIT$C, sep = "_")
		Invar <- MAIT[match(TCRs$V_J_C, MAIT$V_J_C),]

		Invar_labs <- cbind(TCRs$cell, as.character(Invar[,1]))
		Invar_labs <- Invar_labs[!is.na(Invar_labs[,2]),]
		Invar_labs <- Invar_labs[match(seur_obj@meta.data$cell_ID, Invar_labs[,1]),]
		seur_obj@meta.data$Invar_T_Cell <- Invar_labs[,2];
		
		# Add concatenated receptors to metadata
		seur_obj@meta.data$TCRs <- get_receptors(TCRs, seur_obj)[,2]
	}
	if (sample %in% names(BCR_files)) {
		#do the same for BCRs
		# Read in BCR files
		BCRs <- read.delim(BCR_files[[sample]], sep=",", header=T)
		# Fix cell-ids
		BCRs$cell <- fix_cell_id(BCRs$barcode, sample)
		# Add clone-type metadata
		seur_obj@meta.data$BCR_clonetype <- get_clone_types(BCRs, seur_obj)[,2]
		# Add concatenated receptors to metadata
		seur_obj@meta.data$BCRs <- get_receptors(BCRs, seur_obj)[,2]
	}
	obj_list[[sample]] <- seur_obj
}


