do_fgsea <- function(scored_genes, pathways=MSigAll, fdr=0.05, nmax=20){

	res <- fgsea(pathways, scored_genes, minSize=15, maxSize=1000)

	if (sum(!is.na(res$pval) & res$padj < fdr) == 0) {print("No significant enrichments"); return();}

	res <- res[!is.na(res$pval) & res$padj < fdr,]
	res <- res[order(res$NES),]
	res_full <- res;
	if (nrow(res) > nmax) {
		res_pos <- data.frame(res[unlist(res$NES) >0,])
		res_pos <- res_pos[!is.na(unlist(res_pos[,1])),]
		res_neg <- data.frame(res[unlist(res$NES) <0,])
		res_neg <- res_neg[!is.na(unlist(res_neg[,1])),]
		res_pos <- res_pos[order(abs(unlist(res_pos$NES))),]
		res_neg <- res_neg[order(abs(unlist(res_neg$NES))),]
		res <- rbind(res_pos[1:min(nrow(res_pos), nmax),], res_neg[1:min(nrow(res_neg), nmax),])
		res <- res[order(res$NES),]
	}

	size <- abs(res$NES)
	colour <- sign(res$NES)
	col_palette <- c("dodgerblue", "grey50", "firebrick")
	gene_lists <- res[,"leadingEdge"]
	sim_mat <- matrix(0, nrow=nrow(gene_lists), ncol=nrow(gene_lists))
	for (i in 1:nrow(gene_lists)) {
		for (j in i:nrow(gene_lists)) {
			int <- length(intersect(unlist(gene_lists[i,1]), unlist(gene_lists[j,1])))
			uni <- length(union(unlist(gene_lists[i]), unlist(gene_lists[j])))
			sim_mat[i,j] <- int/uni
			sim_mat[j,i] <- int/uni
			colnames(sim_mat) <- unlist(res[,1])
			rownames(sim_mat) <- unlist(res[,1])
		}
	}
	require(igraph)
	G <- simplify(graph_from_adjacency_matrix(sim_mat > 0.1, mode="undirected"))
	plot(G, vertex.color=col_palette[colour+2], vertex.size=size*5, edge.width=2)
	res$cluster <- components(G)$membership
	return(list(rich=res_full, graph=G, vertex_col = col_palette[colour+2], vertex_size = size*5))
}


#### My Varimax ####
RunVarimax <- function(seur_obj, npcs=50, do.scale=FALSE, use.projected = FALSE) {
	require(Seurat)

	# Project pca
	set.seed(1982)
	if (do.scale){
		obj <- ScaleData(obj, features=rownames(obj))
	}
	if (use.projected) {
		#obj <- ProjectDim(obj, reduction="pca", do.center=TRUE, verbose=F) # - always only 5000 genes?
		projected_pca <- obj@assays$SCT@scale.data %*% obj@reductions$pca@cell.embeddings
		projected_pca <- apply(projected_pca, 2, function(x){x/max(abs(x))})
	} else {
		projected_pca <- obj@reductions$pca@feature.loadings
	}

	pca_loading <- projected_pca[,1:npcs]
	#pca_loading <- pca_loading/rowSums(pca_loading)

	set.seed(1029)
	varimax_res <- varimax(pca_loading)
	
	rotated_pca <- obj@reductions$pca@cell.embeddings[,1:npcs] %*% varimax_res$rotmat
	rotated_loadings <- projected_pca[,1:npcs] %*% varimax_res$rotmat
	

	# Prevent inverted colours across components #
	for (i in 1:ncol(rotated_pca)) {
		# Find direction with most variability
		t <- quantile(abs(rotated_pca[,i]), 0.99)
		# count number of cells with extreme values of the component.
		if (sum(rotated_pca[,i] < -1*t) > sum(rotated_pca[,i] > t)) {
			# flip it!
			rotated_pca[,i] = -1 * rotated_pca[,i]
			rotated_loadings[,i] = -1 * rotated_loadings[,i]
		}
	}
	

	# Rank components by variance #
	varimax_var <- apply(rotated_pca, 2, var)
	varimax_var <- varimax_var/sum(varimax_var)
	#reorder <- order(varimax_var, decreasing=TRUE)
	#rotated_pca <- rotated_pca[,reorder]
	#rotated_loadings <- rotated_loadings[,reorder]
	#varimax_var <- varimax_var[reorder]

	# Add cell-scores to metadata for visualization #
	colnames(rotated_loadings) <- paste("RotPC", 1:ncol(rotated_loadings), sep="_")
	colnames(rotated_pca) <- paste("RotPC", 1:ncol(rotated_loadings), sep="_")
	
	obj@meta.data <- cbind(obj@meta.data, rotated_pca[,1:npcs])

	# Add Varimax as a Reduction to the object #
	varimax_pca <- CreateDimReducObject(
				embeddings = rotated_pca, 
				loadings=rotated_loadings, 
				key = "VM_", assay = DefaultAssay(obj), stdev=varimax_var)

	obj@reductions[["varimax"]] <- varimax_pca
	return(obj);
}

Interpret_varimax <- function(obj, cluster_col="Coarse_clusters", metadata_col=c("sample")) {

}

# Function to get the same colours as used by default in Seurat DimPlot
get_seurat_colours <- function(obj, group.by) {
	require(scales)
	identities <- obj@meta.data[,group.by]
	if (class(identities) != "factor") {
		identities <- factor(identities)
	}
	identities <- levels(identities)

	my_color_palette <- hue_pal()(length(identities))
	return(my_color_palette)
#usage:
# TSNEPlot(object = object, do.return = T) + 
# scale_color_manual(values = my_color_palette)

}


get_metadata_outcome <- function(obj){
	tmp <- read.delim("../../Caudate_recip_data_Dec 3_20.csv", sep=",")
	reject <- tmp[ match(obj@meta.data$donor, tmp[,1]), "Post.LT.Rejection"]
	reject <- factor(reject, levels=c("N", "?", "Y"))
	obj@meta.data$trans.rejected <- reject
	return(obj)
}

files <- c("AntiB_harmony_Subcluster.rds", 
		"Cholangiocyte_harmony_Subcluster.rds", 
		"Endo_harmony_Subcluster.rds", 
		"Hepatocyte1_harmony_Subcluster.rds", 
		"Hepatocyte2_harmony_Subcluster.rds", 
		"Macrophage_harmony_Subcluster.rds", 
		"NKT_harmony_Subcluster.rds", 
		"Stellate_harmony_Subcluster.rds")



top_up_gene_lists <- c()
top_down_gene_lists <- c()
gene_lists_names <- c();

ntop <- 50
npcs <- 50;

varimax_gene_loadings_matrix <- c();
obj_list <- list()

require(Polychrome)
donor_colourscheme <- glasbey.colors(25)
donor_colourscheme <- donor_colourscheme[-1]

for (f in files) {
	require(Seurat)
	set.seed(28210)
	tag <- unlist(strsplit(f, "\\/"))[1]
	tag <- unlist(strsplit(f, "_"))[1]
	print(tag)
	obj <- readRDS(f)
	obj <- get_metadata_outcome(obj)
	obj@meta.data$donor_sex <- factor(obj@meta.data$donor_sex)
	cluster_column <- "Coarse_clusters"

	cluster_cols <- get_seurat_colours(obj, group.by=cluster_column)
	cell_clust_col <-  cluster_cols[obj@meta.data[,cluster_column]]
	
	print(dim(obj))
	outfile <- sub("harmony", "varimax", f)


	obj <- RunPCA(obj, features=rownames(obj), npcs=50)

	obj <- RunVarimax(obj)

	saveRDS(obj, outfile)

#	cluster_col_palette <- c("#ebac23", "#b80058", "#008cf9", "#00bbad", 
#					"#d163e6", "#b24502", "#ff9287", "#5954d6", 
#					"#00c6f8", "#878500", "#00a76c", "#bdbdbd")
#	cell_cluster_col <- cluster_col_palette[obj@meta.data[,cluster_column]]

	cell_loadings <- obj@reductions$varimax@cell.embeddings
	gene_loadings <- obj@reductions$varimax@feature.loadings
	all_scores <- c();
	all_coeffs <- c();

	# ID important components.
	for (i in 1:ncol(cell_loadings)) {
		#sex_score <- t.test(cell_loadings[obj@meta.data$donor_sex=="F",i], cell_loadings[obj@meta.data$donor_sex=="M",i])
		res <- summary(glm(cell_loadings[,i] ~ obj@meta.data$donor_sex + obj@meta.data$donor_age_group))
		res2 <- summary(glm(cell_loadings[,i] ~ obj@meta.data[,cluster_column]))
		res3 <- summary(glm(cell_loadings[,i] ~ obj@meta.data$trans.rejected))
		all_scores <- rbind(all_scores, c(res$coefficients[,4], res3$coefficients[,4], res2$coefficients[,4])) # before 4,4,1
		all_coeffs <- rbind(all_coeffs, c(res$coefficients[,1], res3$coefficients[,1], res2$coefficients[,1] ))
		
	}
	rownames(all_scores) <- 1:ncol(cell_loadings)
	all_scores[order(all_scores[,2]),]
	rownames(all_coeffs) <- 1:ncol(cell_loadings)
	all_coeffs[order(all_coeffs[,2]),]


	sex_cols = c("pink", "dodgerblue")
	age_cols = c("grey45", "dodgerblue", "forestgreen")
	reject_cols = c("grey85", "grey65", "grey25")


	varimax_association_plot <- function(component1, component2) {
			layout(rbind(c(1,2,3), c(4,5,5)))
			par(mar=c(4,4,1,1))
			# Sex
			plot(cell_loadings[,component1], cell_loadings[,component2], pch=16,
				col=sex_cols[factor(obj@meta.data$donor_sex)],
				 xlab=paste("RotPC", component1), ylab=paste("RotPC", component2))
			legend("bottomright", levels(factor(obj@meta.data$donor_sex)), fill=sex_cols, bty="n")
			# Age
			plot(cell_loadings[,component1], cell_loadings[,component2], pch=16,
				col=age_cols[factor(obj@meta.data$donor_age_group)],
				xlab=paste("RotPC", component1), ylab=paste("RotPC", component2))
			legend("bottomright", levels(factor(obj@meta.data$donor_age_group)), fill=age_cols, bty="n")
			# Rejection
			plot(cell_loadings[,component1], cell_loadings[,component2], pch=16,
				col=reject_cols[obj@meta.data$trans.rejected],
				xlab=paste("RotPC", component1), ylab=paste("RotPC", component2))
			legend("bottomright", levels(obj@meta.data$trans.rejected), fill=reject_cols, bty="n")
			# Cluster
			plot(cell_loadings[,component1], cell_loadings[,component2], pch=16,
				col=cell_clust_col,
				xlab=paste("RotPC", component1), ylab=paste("RotPC", component2))
			legend("bottomright", levels(obj@meta.data[,cluster_column]), fill=cluster_cols, bty="n")
			# Donor
			plot(cell_loadings[,component1], cell_loadings[,component2], pch=16,
				col=donor_colourscheme[as.numeric(factor(obj@meta.data$donor))],
				xlab=paste("RotPC", component1), ylab=paste("RotPC", component2))
			legend("bottomright", levels(factor(obj@meta.data$donor)), 
					fill=donor_colourscheme[1:length(unique(obj@meta.data$donor))], ncol=2, bty="n")

	}
	varimax_association_boxplot <- function(component, meta_column, colours) {
		boxplot(cell_loadings[,component] ~ obj@meta.data[,cluster_column]+obj@meta.data[, meta_column],
			col=rep(colours, each=length(unique(obj@meta.data[, cluster_column]))),
			las=2, xlab="", ylab=paste("RotPC", component, sep="_"), outline=FALSE)
	}

	pdf(paste(tag, "associations_pval_varimax.pdf", sep="_"), width=9, height=6.5)
	# Sex
	sex_comp <- as.numeric(rownames(all_scores[order(all_scores[,2]),]))
	varimax_association_plot(sex_comp[1], sex_comp[2])
	
	# old
	old_comp <- as.numeric(rownames(all_scores[order(all_scores[,3]),]))
	varimax_association_plot(old_comp[1], old_comp[2])

	# young
	young_comp <- as.numeric(rownames(all_scores[order(all_scores[,4]),]))
	varimax_association_plot(young_comp[1], young_comp[2])

	
	# Rejection
	reject_comp <- as.numeric(rownames(all_scores[order(all_scores[,7]),]))
	varimax_association_plot(reject_comp[1], reject_comp[2])

	par(mfrow=c(4,2))
	par(mar=c(5,4,0,0.5))
	varimax_association_boxplot(sex_comp[1], "donor_sex", sex_cols)
	varimax_association_boxplot(sex_comp[2], "donor_sex", sex_cols)
	varimax_association_boxplot(old_comp[1], "donor_age_group", age_cols)
	varimax_association_boxplot(old_comp[2], "donor_age_group", age_cols)
	varimax_association_boxplot(young_comp[1], "donor_age_group", age_cols)
	varimax_association_boxplot(young_comp[2], "donor_age_group", age_cols)
	varimax_association_boxplot(reject_comp[1], "trans.rejected", reject_cols)
	varimax_association_boxplot(reject_comp[2], "trans.rejected", reject_cols)


	dev.off()
}


# Clusters
FeaturePlot(obj, paste("RotPC_", sex_comp[1], sep=""))
tail(sort(gene_loadings[,sex_comp[1]]), 20)


require(fgsea)
immune_path <- gmtPathways("../../ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("../../ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("../../ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("../../ExternalData/ReactomePathways.gmt")

richments <- do_fgsea(gene_loadings[,1], Hallmark_path, fdr=0.05)
	



## Cluster components:
colnames(top_up_gene_lists) <- gene_lists_names
colnames(top_down_gene_lists) <- gene_lists_names
colnames(varimax_gene_loadings_matrix) <- gene_lists_names



## Pairwise similarity
set.seed(1927)
require(proxy)

intersect_similarity <- function(x, y) {length(intersect(x, y))}
cor_similarity <- function(x, y) {
	keep <- !is.na(x) & !is.na(y)
	cor(x[keep], y[keep])
}



pairs <- t(combn(1:(npcs*length(files)), 2))

sim_mat_int <- matrix(0, nrow=npcs*length(files), ncol=npcs*length(files))
sim_mat_cor <- matrix(0, nrow=npcs*length(files), ncol=npcs*length(files))
direction <- rep(0, nrow(pairs))

for (i in 1:nrow(pairs)) {
	# Intersection
	a <- intersect_similarity(top_up_gene_lists[,pairs[i,1]], top_up_gene_lists[,pairs[i,2]]) # both up
	b <- intersect_similarity(top_down_gene_lists[,pairs[i,1]], top_down_gene_lists[,pairs[i,2]]) # both down
	c <- intersect_similarity(top_up_gene_lists[,pairs[i,1]], top_down_gene_lists[,pairs[i,2]]) # up vs down
	d <- intersect_similarity(top_down_gene_lists[,pairs[i,1]], top_up_gene_lists[,pairs[i,2]]) # down vs up
	if (a+b > c+d) {
		sim_mat_int[pairs[i,1], pairs[i,2]] <- a+b
		direction[i] = 1
	} else {
		sim_mat_int[pairs[i,1], pairs[i,2]] <- c+d
		direction[i] = -1
	}
	# correlation
	s <- cor_similarity(varimax_gene_loadings_matrix[,pairs[i,1]], varimax_gene_loadings_matrix[,pairs[i,2]])
	sim_mat_cor[pairs[i,1], pairs[i,2]] <- s
}

cor_loadings <- varimax_gene_loadings_matrix[rowSums(is.na(varimax_gene_loadings_matrix)) == 0,]

sim_mat_cor <- cor(cor_loadings)

diff_cor <- matrix(1, nrow=npcs*length(files), ncol=npcs*length(files))
diff_cor <- diff_cor-sim_mat_cor
#diff_cor <- diff_cor - t(sim_mat)
diag(diff_cor) <- 0
colnames(diff_cor) <- colnames(top_up_gene_lists)
rownames(diff_cor) <- colnames(top_up_gene_lists)


require("gplots")
hmap <- heatmap.2(diff_cor, scale="none", trace="none", col=colorRampPalette(c("black", "white"))(20), distfun=function(x){as.dist(x)})




diff_int <- matrix(ntop*2, nrow=npcs*length(files), ncol=npcs*length(files))
diff_int <- diff_int-sim_mat_int
diff_int <- diff_int - t(sim_mat_int)
diag(diff_int) <- 0
colnames(diff_int) <- colnames(top_up_gene_lists)
rownames(diff_int) <- colnames(top_up_gene_lists)


# MNN - no don't like this
#adj_mat <- matrix(0, nrow=nrow(diff_int), ncol=ncol(diff_int));
#for (j in 1:length(files)) {
#	columns <- 1:12 + 12*(j-1)
#	sub_mat <-  diff_int[, columns]
#	mnn <-apply(sub_mat, 1, function(x){which(x == min(x) & x < 80)[1]})
#
#	adj_mat[cbind(1:nrow(adj_mat), columns[mnn])] <- 1
#}
#mnn_mat <- adj_mat * t(adj_mat)


require("gplots")
hmap <- heatmap.2(diff_int, scale="none", trace="none", col=colorRampPalette(c("black", "white"))(20), distfun=function(x){as.dist(x)})


# Get clusters
min_size = 3;
cluster_lab <- rep("", ncol(diff_int))
cluster_id = 1;
for (k in ncol(diff_int):1) {
	clusters <- cutree(as.hclust(hmap$rowDendrogram), k = k) 
	c_size <- table(clusters)
	#print(k)
		for (c in names(c_size)[c_size >= min_size]) {
		curr_lab <- cluster_lab[clusters == c]
		if (sum(curr_lab == "") == length(curr_lab)) {
			cluster_lab[clusters == c] <- cluster_id;
			cluster_id <- cluster_id +1;
		} else if (length(unique(curr_lab[curr_lab != ""])) == 1) {
			cluster_lab[clusters == c] <- unique(cluster_lab[clusters == c & cluster_lab != ""]);
		} else {
			#print(curr_lab)
		}
	}
}

set.seed(2829)
library(RColorBrewer)
n <- max(as.numeric(cluster_lab))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cluster_cols <- col_vector[as.numeric(cluster_lab)]

png("heatmap_cluster_legend.png", width=4, height=4,units="in", res=300)
tmp <- cbind(cluster_lab, cluster_cols); tmp <- unique(tmp)
pie(rep(1, nrow(tmp)), col=tmp[,2], label=tmp[,1])
dev.off()

png("heatmap_cluster.png", width=8, height=8, units="in", res=300)
hmap <- heatmap.2(diff_int, RowSideColors=cluster_cols, ColSideColors=cluster_cols, scale="none", trace="none", col=colorRampPalette(c("black", "white"))(20))
dev.off()

names(cluster_lab) <- colnames(diff_int)
saveRDS(list(varimax_gene_loadings_matrix=varimax_gene_loadings_matrix, 
		diff_intersection=diff_int, 
		diff_correlation=diff_cor, 
		clusters_int=cluster_lab,
		clusters_cols_int=cluster_cols,
		heatmap_obj_int=hmap
		), file="Spatial_varimax_clustering.rds")

#

cluster_scores_mat <- c();

id = "12"

tmp <- varimax_gene_loadings_matrix[,cluster_lab == id]
cor(tmp[!is.na(rowSums(tmp)),])

names(cluster_lab)[cluster_lab == id]

tmp <- tmp[!is.na(rowSums(tmp) ),]
scores <- rowMeans(apply(tmp, 2, function(x){x/max(abs(x))}))
scores <- scores[!grepl("^RP", names(scores))]
scores <- scores[!grepl("^MT-", names(scores))]

head(sort(scores), 20)
tail(sort(scores), 20)

#require(gprofiler2)

#ge_out <- gost(names(sort(scores, decreasing=TRUE)), organism="hsapiens", ordered_query=TRUE, correction_method="fdr",
#			measure_underrepresentation=TRUE, sources=c("GO:BP", "TF", "REAC", "KEGG", "WP"), custom_bg = names(scores))


#res <- ge_out$result[ge_out$result$term_size < 500 & ge_out$result$term_size > 20,]

require(fgsea)
WP_gmt <- gmtPathways(gmt.file="../ExternalData/BaderLab25Aug2020/Human_WikiPathways_August_01_2020_symbol.gmt.txt")
kegg_gmt <- gmtPathways(gmt.file="../ExternalData/BaderLab25Aug2020/Human_KEGG_August_01_2020_symbol.gmt.txt")
react_gmt <- gmtPathways(gmt.file="../ExternalData/BaderLab25Aug2020/Human_Reactome_August_01_2020_symbol.gmt.txt")
msigdb_gmt <- gmtPathways(gmt.file="../ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")
iob_gmt <- gmtPathways(gmt.file="../ExternalData/BaderLab25Aug2020/Human_IOB_August_01_2020_symbol.gmt.txt")




out <- fgsea(pathway=WP_gmt, stats = sort(scores, decreasing=TRUE), minSize=15, maxSize=1000); out <- out[order(out$padj),]
out2 <- fgsea(pathway=kegg_gmt, stats = sort(scores, decreasing=TRUE), minSize=15, maxSize=1000); out2 <- out2[order(out2$padj),]
out3 <- fgsea(pathway=react_gmt, stats = sort(scores, decreasing=TRUE), minSize=15, maxSize=1000); out3 <- out3[order(out3$padj),]
out4 <- fgsea(pathway=msigdb_gmt, stats = sort(scores, decreasing=TRUE), minSize=15, maxSize=1000); out4 <- out4[order(out4$padj),]
out5 <- fgsea(pathway=iob_gmt, stats = sort(scores, decreasing=TRUE), minSize=15, maxSize=1000); out5 <- out5[order(out5$padj),]



out3[out3$padj < 0.05 & out3$NES > 0,]





########################



sort(table(clusters))
colnames(gene_lists)[clusters==8]

genes <- c(unlist(gene_lists[1:50,clusters==46]))

genes1 <- c(unlist(gene_lists[1:50, colnames(gene_lists) %in% c("PSC011_4_A1_2", "PSC011_4_B1_2", "PSC011_4_C1_2", "PSC011_4_D3")]))#,
#		unlist(gene_lists[51:100, colnames(gene_lists) %in% c("PSC011_4_D1_6")]))

genes2 <- c(unlist(gene_lists[51:100, colnames(gene_lists) %in% c("PSC011_4_A1_4", "PSC011_4_B1_4", "PSC011_4_C1_4")]),
		unlist(gene_lists[1:50, colnames(gene_lists) %in% c("PSC011_4_D1_6")]))

### Correlation sim/diff
tmp <- varimax_gene_
require(proxy)
