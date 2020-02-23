set.seed(1189)

int.obj <- readRDS("integration_harmony_plus_analysis.rds")


# Get auto-annotation for fine & coarse clusters

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

# Output:

# Auto-anno table:
# % cells/ref cluster - scmap cluster, scmap cell
# marker enrichment for markers from reference clusters.

# Cluster x Annotation table

## Cluster specific Marker genes ###
# Get markers for fine & coarse clusters
# problem: seurat integrated only has 1000 genes....
# harmony one only has ~3,000 - background for GSEA!!!

my_rowMeans <- function(x) {
	if (!is.null(ncol(x))) {
		if (ncol(x) > 1) {
			return(Matrix::rowMeans(x))
		}
	}
	return(x);
}
my_rowSums <- function(x) {
	if (!is.null(ncol(x))) {
		if (ncol(x) > 1) {
			return(Matrix::rowSums(x))
		}
	}
	return(x);
}

## Significant with MAST
set.seed(8817)

require("MAST")

sca <- FromMatrix(as.matrix(int.obj@assays$RNA@data)) ##### <---- not sparse this is very Mem heavy!!

ngenes <- colSums(int.obj@assays$RNA@counts > 0)
ncounts <- colSums(int.obj@assays$RNA@counts)
donor <- int.obj@meta.data$donor
cluster <- int.obj@meta.data$Fine_clusters

# This is slow with so many cells.... Perhaps move to H4h cluster.
mod <- zlm(~cluster + donor + ngenes + ncounts, sca)


saveRDS(mod, "harmony_MASTmodel.rds")

set.seed(3891)

res1 <- summary(mod, doLRT="cluster1")
res <- res1$datatable[res1$datatable$contrast=="cluster1" & res1$datatable$component %in% c("C","D"),]
res$fdr <- p.adjust(unlist(res[,4]), method="fdr")

con <- res[res$component =="C",]
dis <- res[res$component =="D",]

sig <- (con$fdr < 0.01 | dis$fdr < 0.01) & sign(con$coef) == sign(dis$coef)
out <- cbind(con[sig,], dis[sig,])
out$overall <- out[,7]+out[,16]
out <- out[order(out$overall, decreasing=T),]
write.table(out, file=paste(this_contrast,"_MAST_output.txt", sep=""), col.names=T, row.names=F)
 


res2 <- summary(mod, doLRT="cluster2")
res3 <- summary(mod, doLRT="cluster3")
res4 <- summary(mod, doLRT="cluster4")
res5 <- summary(mod, doLRT="cluster5")
res6 <- summary(mod, doLRT="cluster6")

## Pseudobulks ##
make_pseudobulk_table <- function(mat, clusters, indis) {
	c <- split(seq(ncol(mat)), clusters);
	donor_expr <- function(MAT, ds) {
		d <- split(seq(ncol(MAT)), ds);
		mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
		return(mus);
	}
	pseudobulk_expr <- lapply(c, function(clust) {
		d_expr <- donor_expr(mat[,clust], indis[clust]);
		return(d_expr);
	})
	return(pseudobulk_expr);
}

make_pseudobulk_v2 <- function(mat, clusters, indis) {
	tab <- apply(mat, 1, function(x){tapply(x, list(clusters, donors), mean)})
	rownames(tab) <- paste(rep(levels(as.factor(clusters)), each=length(levels(as.factor(donors)))), 
			rep(levels(as.factor(donors)), times=length(levels(as.factor(clusters)))), sep="_")
	return(t(tab))
}


## Ranked list using my method ##

# Matrix of relative expression levels across clusters
# per cluster get means per donor and weight by overall donor freqs

get_rel_expression <- function(mat, clusters, donors) {
	c <- split(seq(ncol(mat)), clusters);
	donor_freqs <- table(donors)/length(donors)
	# avg expression per donor in this cluster
	donor_expr <- function(MAT, ds) {
		d <- split(seq(ncol(MAT)), ds);
		mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
		return(mus);
	}
	clust_expr <- sapply(c, function(clust) {
		d_expr <- donor_expr(mat[,clust], donors[clust]);
		# weight by overall frequency of donors
		freqs <- donor_freqs[match(colnames(d_expr), names(donor_freqs))]
		freqs <- as.vector(freqs)/sum(freqs)
		c_expr <- my_rowSums(t(t(d_expr)*freqs))
		return(c_expr);
	})
	return(clust_expr)
}

source("Setup_autoannotation.R")

## Ranked list
cluster_means <- get_rel_expression(int.obj@assays$RNA@data, int.obj@meta.data$Fine_clusters, int.obj@meta.data$donor)
cluster_detect <- get_rel_expression(int.obj@assays$RNA@counts > 0, int.obj@meta.data$Fine_clusters, int.obj@meta.data$donor)

mark_mean <- my_markers(cluster_means)
mark_detect <- my_markers(cluster_detect);

