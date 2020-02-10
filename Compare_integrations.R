require("Seurat")
int.obj <- readRDS("integration_harmony_plus_analysis.rds")

int.obj <- readRDS("integration_v2_plus_analysis.rds")

cluster_entropy <- function(x) {
	freqs <- table(x)/length(x);
	Shannon <- -sum(freqs * log2(freqs))
	return(Shannon);
}


# Entropy of Auto cell-type annotations (at two levels of hierarchy) across fine clusters.
general_cell_type_ids <- as.character(int.obj@meta.data$scmap_anno2)
general_cell_type_ids[general_cell_type_ids %in% c("CD3abTcells","gdTcells2")] <- "Tcells"
general_cell_type_ids[general_cell_type_ids %in% c("MatureBcells","AntibodysecretingBcells")] <- "Bcells"
general_cell_type_ids[general_cell_type_ids %in% c("PeriportalLSECs","CentralvenousLSECs")] <- "LSECs"
general_cell_type_ids[general_cell_type_ids %in% c("interzonalHep","PeriportalHep", "UnidentifiedHep", "PericentralHep")] <- "Hep"
general_cell_type_ids[general_cell_type_ids %in% c("Non-inflammatoryMacrophages", "InflamatoryMacrophages")] <- "Mac"



anno_entropy_fine <- sapply(split(int.obj@meta.data$scmap_anno2, int.obj@meta.data$Fine_clusters), cluster_entropy)
anno_entropy_coarse <- sapply(split(int.obj@meta.data$scmap_anno2, int.obj@meta.data$Coarse_clusters), cluster_entropy)

# Entropy of donors across fine clusters.
donor_entropy_fine <- sapply(split(int.obj@meta.data$donor, int.obj@meta.data$Fine_clusters), cluster_entropy)
donor_entropy_coarse <- sapply(split(int.obj@meta.data$donor, int.obj@meta.data$Coarse_clusters), cluster_entropy)

png("harmony_cluster_entropy_fine.png")
barplot(rbind(anno_entropy_fine, donor_entropy_fine), beside=T, horiz=T, las=1, xlab="Entropy", ylab="Cluster", col=c("salmon", "navy"), bty="L")
legend("bottom", c("Annotation", "Donor"), fill=c("salmon", "navy"), bty="n", horiz=T, inset=c(0,1), xpd=T)
dev.off()

png("harmony_cluster_entropy_coarse.png")
barplot(rbind(anno_entropy_coarse, donor_entropy_coarse), beside=T, horiz=T, las=1, xlab="Entropy", ylab="Cluster", col=c("salmon", "navy"), bty="L")
legend("bottom", c("Annotation", "Donor"), fill=c("salmon", "navy"), bty="n", horiz=T, inset=c(0,1), xpd=T)
dev.off()
print(c(median(anno_entropy_fine), median(anno_entropy_coarse)))
print(c(median(donor_entropy_fine), median(donor_entropy_coarse)))

rm(int.obj)
