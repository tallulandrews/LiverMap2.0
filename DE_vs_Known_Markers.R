proj_name = "Harmony_Integrated"
set.seed(82710)
n_clusters = 23
cluster_names <- 1:n_clusters -1

anno_genes <- read.table("Celltype_markers_Liver.csv", sep=",", header=T)
source("~/R-Scripts/Ensembl_Stuff.R")
anno_genes$human <- as.character(General_Map(as.character(anno_genes[,1]), in.org="Mmus", in.name="symbol", out.org="Hsap", out.name="symbol"))
anno_genes$human[anno_genes$Species == "Human"] <- as.character(anno_genes$Gene[anno_genes$Species == "Human"])

output_list <- list();

for (c in cluster_names) {
	## for each cluster
	file <- paste("DEGeneFiles/", proj_name, "_", c, "_DE.rds", sep="")
	if (!file.exists(file)) {next;}
	obj <- readRDS(file)

	wilcox_res <- obj[["seur_wilcox"]]
	wilcox_res$dir <- sign(wilcox_res$avg_logFC)
	res_anno <- anno_genes[match(rownames(wilcox_res), anno_genes$human),]
	all_out <- cbind(wilcox_res, res_anno);
	all_out <- all_out[!is.na(all_out$human),]
	all_out$cluster <- rep(c, nrow(all_out));

	output_list[[as.character(c)]] <- all_out;
#outfile <- paste(proj_name, c, "DE_vs_marker.csv", sep="_")
#write.table(tab, outfile, sep=",", row.names=T, col.names=T)
}

saveRDS(output_list, file=paste(proj_name, "known_de_markers.rds", sep="_"));
