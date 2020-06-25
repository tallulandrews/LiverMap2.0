## Manually Annotated Global Clustering ##
source("~/scripts/LiverMap2.0/Colour_Scheme.R") # custom colour scheme

liver.integrated <- readRDS("All_genes_Merged_obj_v2_with_Analysis2.rds")
cluster_manual_Anno <- cbind(0:27, c("Hepatocyte", "Hepatocyte","Hepatocyte","Hepatocyte","Hepatocyte", "CD3abTcells", "Hepatocyte","Hepatocyte","cvLSECs", "InflammatoryMacrophages", "CD3abTcells", "Non-inflammatoryMacrophages", "cvLSECs", "Non-inflammatoryMacrophages", "Eryth", "PeriLSECs", "AntiBcell", "Hepatocyte","Hepatocyte", "Stellate", "gdTcells2", "Stellate", "Cholangiocyte", "MatBcell", "Hepatocyte", "Macrophage", "Hepatocyte", "Eryth"))

cluster_manual_Anno <- cbind(cluster_manual_Anno, map_cell_types(cluster_manual_Anno[,2]))

# Cluster ID labelled figures
agg_coord_by_cluster <- function(coords, clusters) {
	x <- split(seq(nrow(coords)), clusters)
	result <- sapply(x, function(a) apply(coords[a,],2,median))
	return(result)
}

tsne_lab_pos <- agg_coord_by_cluster(liver.integrated@reductions$tsne@cell.embeddings, liver.integrated@meta.data$Fine_clusters)
umap_lab_pos <- agg_coord_by_cluster(liver.integrated@reductions$umap@cell.embeddings, liver.integrated@meta.data$Fine_clusters)
lab_id <- colnames(tsne_lab_pos)

liver.integrated@meta.data$short_cluster_anno <- factor(cluster_manual_Anno[match(liver.integrated@meta.data$Fine_clusters, cluster_manual_Anno[,1]),3])

require("ggplot2")
new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
#liver.integrated@meta.data$short_cluster_anno <- factor(map_cell_types(liver.integrated@meta.data$cluster_quickanno), levels=new_colour_scheme[,1]);
new_colour_scheme <- new_colour_scheme[new_colour_scheme[,1] %in% liver.integrated@meta.data$short_cluster_anno,]

print("plotting")
#png("Autoanno_label_harmony_umap.png", width=7.5, height=6, units="in", res=300)
png("Manual_label_harmony_umap.png", width=7.5, height=6, units="in", res=300)
DimPlot(liver.integrated, reduction="umap", group.by="short_cluster_anno", pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=lab_id, colour="grey35")
dev.off()


#png("Autoanno_label_harmony_tsne.png", width=7.5, height=6, units="in", res=300)
png("Manual_label_harmony_tsne.png", width=7.5, height=6, units="in", res=300)
DimPlot(liver.integrated, reduction="tsne", group.by="short_cluster_anno", pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=tsne_lab_pos[1,], y=tsne_lab_pos[2,], label=lab_id, colour="grey35")
dev.off()

# T/NK = 10, 5, + 20
# Endo = 15, 12, 8, 19, 21
# Mac = 9, 11, 13, 25
# Hep Traj = 4, 1, 2, 24, 17, 3, 6, 7, 26
# Doublets? 20
# Hep / Erythroid = 0, 14, 27

source("~/scripts/LiverMap2.0/Map2.2_add_metadata.R")
metadata <- get_metadata(liver.integrated$donor)

tab <- table(metadata$AGE, liver.integrated@meta.data$short_cluster_anno)

test <- fisher.test(tab[,c("AntiBcell", "MatBcell")])
test <- fisher.test(tab[,c("InfMac", "Macrophage", "NonInfMac")])

