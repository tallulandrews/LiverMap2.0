#a <- readRDS("All_genes_Merged_obj_v2_with_Analysis.rds")
require("Seurat")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")
#donor <- a@meta.data$sample
donor <- factor(rep("all", length(a@meta.data$sample)))
image_name <-"All_gene_merged_clustermeans_heatmap.png" 
image_name2 <-"All_gene_merged_clustermeans_tcell_heatmap.png" 

#bulks <- get_pseudobulk_means(a@assays$RNA@data, a@meta.data$Fine_clusters, donor)
bulks <- group_rowmeans(a@assays$RNA@data, a@meta.data$Fine_clusters)
liver_markers <- read.table("~/movingfiles/Celltype_markers_Liver.csv", sep=",", header=T, stringsAsFactors=FALSE)
liver_markers <- liver_markers[liver_markers[,2] != "",]
immune_markers <- read.table("~/movingfiles/Celltype_markers_Immune.csv", sep=",", header=T, stringsAsFactors=FALSE)
all_markes <- rbind(liver_markers[,1:2], immune_markers[,1:2])
all_markes <- unique(all_markes)
repeats <- all_markes[duplicated(all_markes[,1]),1]
all_markes$lab_length<- sapply(all_markes[,2], stringr::str_length)
for( i in repeats) {
	tmp <- all_markes[all_markes[,1] == i,]
	keep <- which(tmp[,3] == min(tmp[,3]))[1]
	all_markes[all_markes[,1]==i,2] <- tmp[keep,2]
	#print(paste(i, tmp[keep,2]));
}

all_markes <- unique(all_markes);
require("stringr")
all_markes[,2] <- stringr::str_to_lower(all_markes[,2])
all_markes <- unique(all_markes);

bulks <- bulks[,colSums(bulks) > 1]
all_markes <- all_markes[all_markes[,1] %in% rownames(bulks),]

colours <- rep("grey65", nrow(all_markes));
colours[grepl("hep", all_markes[,2])] <- "firebrick"
colours[grepl("stellate", all_markes[,2])] <- "forestgreen"
colours[grepl("mac", all_markes[,2])] <- "dodgerblue"
colours[grepl("b cell", all_markes[,2])] <- "goldenrod1"
colours[grepl("plasma", all_markes[,2])] <- "goldenrod1"
colours[grepl("t cell", all_markes[,2])] <- "violet"
colours[grepl("t-cell", all_markes[,2])] <- "violet"
colours[grepl("tcell", all_markes[,2])] <- "violet"
colours[grepl("lsec", all_markes[,2])] <- "darkorange"
colours[grepl("non", all_markes[,2])] <- "navy"
colours[grepl("kupffer", all_markes[,2])] <- "navy"
colours[grepl("endo", all_markes[,2])] <- "salmon"

heatdata <- bulks[match(all_markes[,1], rownames(bulks)),]
require(gplots)

png(image_name, width=11, height=11, units="in", res=200)
heatmap.2(heatdata, scale="row", RowSideColors=colours, trace="none")
dev.off()

#t-cell receptors

tcell_recept <- immune_markers[immune_markers$X != "",]
tcell_marks <- c(as.character(tcell_recept[,1]),"CD3D", "CD3G", "CD3E", "CD247")

heatdata <- bulks[match(tcell_marks, rownames(bulks)),]
png(image_name2, width=11, height=11, units="in", res=200)
heatmap.2(heatdata, scale="row", trace="none")
dev.off()
