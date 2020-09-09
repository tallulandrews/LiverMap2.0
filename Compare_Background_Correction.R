require(Seurat)
require(stringr)

files <- Sys.glob("/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/SoupX/*emptyDrops_table_SeurObj_SoupX.rds")
PREFIX = "SoupXEmpty"; # SoupX = CR+SoupX

merged_obj <- NULL;
universal_genes <- c(-1)
all_genes <- c();
obj_list <- list()
all_genes <- c();
for (i in 1:length(files)) {
	print(i);
	print(files[i])
        obj <- readRDS(files[i]);

	soupcorrect <- obj@assays$SoupCorrected
	#soupcorrect <- obj@assays$RNA@counts
	#soupcorrect <- obj@assays$My_Corrected
	# Delete unnecessary
	obj@assays <- list(RNA=obj@assays$RNA);
	
	soupcorrect <- soupcorrect[!is.na(rownames(soupcorrect)),]
	obj <- obj[rownames(obj) %in% rownames(soupcorrect),]
	obj <- obj[match(rownames(soupcorrect), rownames(obj)),]
	obj <- obj[, match(colnames(soupcorrect), colnames(obj))]
	
	tmp <- paste("C", substring(files[[i]], regexpr("C[0-9]",files[[i]])+1), sep="")
	this_name <- gsub("_Anno_SeurObj_SoupX.rds", "", tmp)

	obj@assays$RNA@counts <- soupcorrect;
	obj <- obj[,Matrix::colSums(soupcorrect >0)>250]
	obj <- obj[Matrix::rowSums(soupcorrect >0)>10,]

        # fix some auto-annotation issues
        #obj@meta.data$consistent_labs <- as.character(obj@meta.data$consistent_labs);
        #anno1 <- obj@meta.data$scmap_cluster_anno
        #anno2 <- obj@meta.data$scmap_cell_anno
        #mac_groups <- c("Non-inflammatoryMacrophages", "InflamatoryMacrophages")
        #g_anno <- as.character(obj@meta.data$general_labs)
        #g_anno[g_anno == "ambiguous" & anno1 %in% mac_groups & anno2 %in% mac_groups] <- "Macrophage"
        #obj@meta.data$general_labs <- g_anno
        #null_vals <- is.na(obj@meta.data$consistent_labs);
        #obj@meta.data$consistent_labs[null_vals] <- g_anno[null_vals]

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


        obj@meta.data$cell_barcode <- colnames(obj);
        if (length(all_genes) == 0) {
                all_genes <- rownames(obj);
        } else {
                all_genes <- union(all_genes, rownames(obj));
        }


	# merge them together
        if (i == 1) {
                merged_obj <- obj
                universal_genes <- as.character(rownames(obj))
        } else {
                merged_obj <- merge(merged_obj, y=obj, add.cell.ids=c("", this_name), project="LiverMap")
                universal_genes <- intersect(universal_genes, as.character(rownames(obj)))
        }
}


# Merge Datasets
#### Merging does not merge individually scaled datasets!!

merged_obj@misc$universal_genes <- universal_genes;
merged_obj@misc$creation_date <- date();

print(PREFIX)
print(dim(merged_obj))
saveRDS(merged_obj, paste(PREFIX, "merged_obj.rds", sep="_"))

# Quick Cluster
res=1
npcs=20
nkNN=20
set.seed(9428)

merged_obj <- merged_obj[rownames(merged_obj) %in% universal_genes,]
merged_obj <- Seurat::NormalizeData(merged_obj, verbose = FALSE)
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
merged_obj <- ScaleData(merged_obj, verbose = FALSE)
merged_obj <- RunPCA(merged_obj, pc.genes = VariableFeatures(merged_obj), npcs = npcs, verbose = FALSE)

# Clustering
merged_obj <- Seurat::FindNeighbors(merged_obj, dims = 1:npcs)
merged_obj <- Seurat::FindClusters(merged_obj, resolution = res, k.param=nkNN)

tab <- table(merged_obj$donor, merged_obj$seurat_clusters)

shannon <- function(vec) {
	prop <- vec/sum(vec);
	prop <- prop[prop>0];
	return(-sum(prop*log2(prop)))
}

simpson <- function(vec) {
	prop <- vec/sum(vec)
	return(1-sum(prop*prop))
}

overall <- rowSums(tab)

cluster_simp <- apply(tab, 2, simpson)
all_simp <- simpson(overall)

png(paste(PREFIX,"merged_barplot.png", sep="_"), width=5, height=5, units="in", res=150)
barplot(c(all_simp, cluster_simp), names=c("all", as.character(1:ncol(tab))), xlab="cluster", ylab="Simpson diversity", main="");
dev.off()

print(all_simp)
print(mean(cluster_simp));
