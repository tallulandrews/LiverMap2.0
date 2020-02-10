
int.obj <- readRDS("integration_harmony_plus_analysis.rds")


cluster_assign <- function(x) {
	freqs <- table(x)/length(x);
	return(freqs[freqs==max(freqs)])
}

general_cell_type_ids <- as.character(int.obj@meta.data$scmap_anno2)
general_cell_type_ids[general_cell_type_ids %in% c("CD3abTcells","gdTcells2")] <- "Tcells"
general_cell_type_ids[general_cell_type_ids %in% c("MatureBcells","AntibodysecretingBcells")] <- "Bcells"
general_cell_type_ids[general_cell_type_ids %in% c("PeriportalLSECs","CentralvenousLSECs")] <- "LSECs"
general_cell_type_ids[general_cell_type_ids %in% c("interzonalHep","PeriportalHep", "UnidentifiedHep", "PericentralHep")] <- "Hepatocyte"
general_cell_type_ids[general_cell_type_ids %in% c("Non-inflammatoryMacrophages", "InflamatoryMacrophages")] <- "Macrophage"

assign_specific <- sapply(split(int.obj@meta.data$scmap_anno2, int.obj@meta.data$Fine_clusters), cluster_assign)
assign_general <- sapply(split(general_cell_type_ids, int.obj@meta.data$Fine_clusters), cluster_assign)


Cluster_guesses <- paste(0:(length(assign_general)-1),".",rep("Unknown", length(assign_general)), sep="");
Cluster_guesses[assign_general > 0.5] <- names(assign_general)[assign_general>0.5]
Cluster_guesses[assign_specific > 0.5] <- names(assign_specific)[assign_specific>0.5]


# Get markers for fine & coarse clusters
# problem: seurat integrated only has 1000 genes....
# harmony one only has ~3,000 - background! for GSEA!!!
source("Setup_autoannotation.R")

## Ranked list
cluster_means <- my_row_mean_aggregate(int.obj@assays$RNA@data, int.obj@meta.data$Fine_clusters)
cluster_detect <- my_row_mean_aggregate(int.obj@assays$RNA@counts > 0, int.obj@meta.data$Fine_clusters)

mark_mean <- my_markers(cluster_means)
mark_detect <- my_markers(cluster_detect);

# Output:

# Auto-anno table:
# % cells/ref cluster - scmap cluster, scmap cell
# marker enrichment for markers from reference clusters.




## Significant with MAST / hypergeometric

require("MAST")

sca <- FromMatrix(as.matrix(int.obj@assays$RNA@data)) ##### <---- not sparse this is very Mem heavy!!

ngenes <- colSums(int.obj@assays$RNA@counts > 0)
ncounts <- colSums(int.obj@assays$RNA@counts)
donor <- int.obj@meta.data$donor
cluster <- int.obj@meta.data$Fine_clusters

mod <- zlm(~cluster + donor + ngenes + ncounts, sca)


# Matrix of relative expression levels across clusters
# per cluster get means per donor and weight by overall donor freqs
donor_freqs <- table(donor)/length(donor)

# Get auto-annotation for fine & coarse clusters


