# Manual Wilcox DE test

run_wilcox <- function(obj, binned, ident.1=min(binned), ident.2=NULL) {
	if (!is.null(ident.2)) {
		de_out <- t( apply(obj@assays$RNA@data, 1, function(x){my_wilcox(x[binned==ident.1], 
													x[binned==ident.2])} ) )
	} else {
		de_out <- t( apply(obj@assays$RNA@data, 1, function(x){my_wilcox(x[binned==ident.1], 
													x[binned!=ident.1])} ) )
	}
	colnames(de_out) <- c("log2fc", "mean.1", "mean.2", "pct.1", "pct.2", "AUC", "p.value")
	de_out <- data.frame(de_out)
	de_out$q.value <- p.adjust(de_out$p.value, method="fdr")
	return(de_out)
}



my_wilcox <- function(x, y) {
	res <- wilcox.test(x, y)
	AUC <- res$statistic/(length(x)*length(y))
	mu_x <- mean(x)
	mu_y <- mean(y)
	if (mu_x == 0) {
		mu_x <- min(x[x>0])/(length(x)+1)
	}
	if (mu_y == 0) {
		mu_y <- min(y[y>0])/(length(y)+1)
	}
	log2fc <- log2( (exp(mean(x))-0.99)/(exp(mean(y))-0.99) )
	
	pct.x <- sum(x>0)/length(x)
	pct.y <- sum(y>0)/length(y)
	mu_x[pct.x == 0] <- 0
	mu_y[pct.y == 0] <- 0
	log2fc[mu_x == 0] <- 0
	
	return( c(log2fc, mu_x, mu_y, pct.x, pct.y, AUC, res$p.value) )
}

# Wrapper for fgsea that incorporates some aspects of Cytoscape

do_fgsea <- function(scored_genes, pathways=MSigAll, fdr=0.05, nmax=20, jaccard=0.25){
	# scored_genes = named vector of scores, names = gene ID.
	# pathways must be formated as read in by fgsea::gmtPathways
	# fdr = multiple testing correction significance level
	# nmax = maximum number of pathways to plot in the graph figure - equal number up/down
	# jaccard = overlap theshold for linking two terms together in the gragh
	set.seed(2910)
	# Clean up scored genes
	scored_genes[scored_genes > 0 & !is.finite(scored_genes)] <- max(scored_genes[is.finite(scored_genes)])+1
	scored_genes[scored_genes < 0 & !is.finite(scored_genes)] <- min(scored_genes[is.finite(scored_genes)])-1

	# Run fgsea
	res <- fgsea(pathways, scored_genes, minSize=15, maxSize=1000, nperm=100000)

	# Error catching: no significant enrichments
	if (sum(!is.na(res$pval) & res$padj < fdr) == 0) {print("No significant enrichments"); return();}

	# Cleaning up NAs
	res <- res[!is.na(res$pval) & res$padj < fdr,]
	res <- res[order(res$NES),]
	res_full <- res;
	if (nrow(res) > nmax*2) {
		res_pos <- data.frame(res[unlist(res$NES) >0,])
		res_neg <- data.frame(res[unlist(res$NES) <0,])
		res_pos <- res_pos[order(abs(unlist(res_pos$NES)), decreasing=T),]
		res_neg <- res_neg[order(abs(unlist(res_neg$NES)), decreasing=T),]
		res <- rbind(res_pos[1:min(nrow(res_pos), nmax),], res_neg[1:min(nrow(res_neg), nmax),])
		res <- res[order(res$NES),]
		if (nrow(res_neg) == 0) {
			res <- res_pos[1:min(nrow(res_pos), nmax),]
			res <- res[order(res$NES),]
		}
		if (nrow(res_pos) == 0) {
			res <- res_pos[1:min(nrow(res_neg), nmax),]
			res <- res[order(res$NES),]			
		}
	}

	size <- abs(res$NES)
	colour <- sign(res$NES)
	col_palette <- c("dodgerblue", "grey50", "firebrick")
	gene_lists <- res[,"leadingEdge"]
	if (! is.null(dim(gene_lists))) {
		gene_lists_tmp <- list()
		for (i in 1:nrow(gene_lists)) {
			gene_lists_tmp[[i]] <- unlist(gene_lists[i,1])
		}
		gene_lists <- gene_lists_tmp;
	}
	sim_mat <- matrix(0, nrow=length(gene_lists), ncol=length(gene_lists))
	trim <- c();
	
	for (i in 1:length(gene_lists)) {
		for (j in i:length(gene_lists)) {
			int <- length(intersect(unlist(gene_lists[[i]]), unlist(gene_lists[[j]])))
			uni <- length(union(unlist(gene_lists[[i]]), unlist(gene_lists[[j]])))
			sim_mat[i,j] <- int/uni
			sim_mat[j,i] <- int/uni
			if (int/uni > 0.8 & i != j) {
				trim <- c(trim, j);
			}
		}
	}
	colnames(sim_mat) <- unlist(res[,1])
	rownames(sim_mat) <- unlist(res[,1])

	require(igraph)
	G <- simplify(graph_from_adjacency_matrix(sim_mat > jaccard, mode="undirected"))
	plot(G, vertex.color=col_palette[colour+2], vertex.size=size*5, edge.width=2)
	res$cluster <- components(G)$membership
	return(list(rich=res_full, graph=G, vertex_col = col_palette[colour+2], vertex_size = size*5))
}



# These are all tested and debugged.

# Wrapper for rowMeans to ensure using a version appropriate
#   for sparse matrices, and allowing for only one row or column
my_rowMeans <- function(x) {
        if (!is.null(ncol(x))) {
                if (ncol(x) > 1) {
                        return(Matrix::rowMeans(x))
                }
		if (ncol(x) == 0) {
			return(rep(NA, nrow(x)));
		}
        }
        return(x);
}
# Wrapper for rowSums to ensure using a version appropriate
#   for sparse matrices, and allowing for only one row or column
my_rowSums <- function(x) {
        if (!is.null(ncol(x))) {
                if (ncol(x) > 1) {
                        return(Matrix::rowSums(x))
                }
		if (ncol(x) == 0) {
			return(rep(NA, nrow(x)));
		}
        }
        return(x);
}
# Wrapper for colMeans to ensure using a version appropriate
#   for sparse matrices, and allowing for only one row or column
my_colMeans <- function(x) {
        if (!is.null(nrow(x))) {
                if (nrow(x) > 1) {
                        return(Matrix::colMeans(x))
                }
		if (nrow(x) == 0) {
			return(rep(NA, ncol(x)));
		}
        }
        return(x);
}
# Wrapper for colSums to ensure using a version appropriate
#   for sparse matrices, and allowing for only one row or column
my_colSums <- function(x) {
        if (!is.null(nrow(x))) {
                if (nrow(x) > 1) {
                        return(Matrix::colSums(x))
                }
		if (nrow(x) == 0) {
			return(rep(NA, ncol(x)));
		}
        }
        return(x);
}

# Row means or row sums by groups.
group_rowmeans <- function(MAT, group_labs, type=c("mean","sum")) {
        d <- split(seq(ncol(MAT)), group_labs);
	if (type[1] == "mean") {
	        mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
	} else {
	        mus <- sapply(d, function(group) my_rowSums(MAT[,group]))
	} 
        return(mus);
}

# Row var by groups.
group_rowvars <- function(MAT, group_labs) {
        d <- split(seq(ncol(MAT)), group_labs);
	mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
	sq_mus <- sapply(d, function(group) my_rowMeans(MAT[,group]^2))
	vars <- sq_mus - mus^2
        return(vars);
}

# Col means or col sums by groups.
group_colmeans <- function(MAT, group_labs, type=c("mean", "sum")) {
        d <- split(seq(nrow(MAT)), group_labs);
	if (type[1] == "mean") {
        	mus <- sapply(d, function(group) my_colMeans(MAT[group,]))
	} else {
        	mus <- sapply(d, function(group) my_colSums(MAT[group,]))
	}
        return(mus);
}

# Average expression of a matrix by cluster avoiding biases due
# due to different numbers of cells by donor across clusters.
#
# - default weights expression by the overall frequency of donors 
#   across the whole dataset
# setting weight to FALSE gives equal weight to each donor.
get_rel_expression <- function(mat, clusters, donors, weight=TRUE) {
        c <- split(seq(ncol(mat)), clusters);
        donor_freqs <- table(donors)/length(donors)
        # avg expression per donor in this cluster
        clust_expr <- sapply(c, function(clust) {
                d_expr <- group_rowmeans(mat[,clust], donors[clust]);
		if (weight) {
                	# weight by overall frequency of donors
                	freqs <- donor_freqs[match(colnames(d_expr), names(donor_freqs))]
                	freqs <- as.vector(freqs)/sum(freqs)
		} else {
			freqs <- rep(1, ncol(d_expr))
		}
                c_expr <- my_rowSums(t(t(d_expr)*freqs))
                return(c_expr);
        })
        return(clust_expr)
}

# Table of total expression of cells from each donor in each cluster 
#  - for edgeR
get_pseudobulk <- function(mat, clusters, donors) {
        c <- split(seq(ncol(mat)), clusters);
        donor_freqs <- table(donors)/length(donors)
        # avg expression per donor in this cluster
        clust_expr <- lapply(c, function(clust) {
                d_expr <- group_rowmeans(mat[,clust], donors[clust], type="sum");
		if(is.null(dim(d_expr))) {
			l <- sapply(d_expr, length)
			keep <- which(l == nrow(mat))
			d_expr <- matrix(d_expr[[keep]], ncol=length(keep), byrow=FALSE);
			rownames(d_expr) <- rownames(mat);
			colnames(d_expr) <- paste(clusters[clust[1]], levels(donors)[keep], sep="_")
		} else {
			colnames(d_expr) <- paste(clusters[clust[1]], colnames(d_expr), sep="_")
		}
                return(d_expr);
        })
	out <- clust_expr[[1]];
	for (i in 2:length(clust_expr)) {
		c_names <- c(colnames(out), colnames(clust_expr[[i]]))
		out <- cbind(out, clust_expr[[i]]);
		if (is.null(dim(out))){
			out <- matrix(out, ncol=1)
			rownames(out) <- rownames(mat)
		}
		colnames(out) <- c_names
	}
        return(out)
}

# Table of mean expression of cells from each donor in each cluster 
#   - for heatmap
get_pseudobulk_means <- function(mat, clusters, donors) {
        c <- split(seq(ncol(mat)), clusters);
        donor_freqs <- table(donors)/length(donors)
        # avg expression per donor in this cluster
        clust_expr <- lapply(c, function(clust) {
                d_expr <- group_rowmeans(mat[,clust], donors[clust], type="mean");
		colnames(d_expr) <- paste(clusters[clust[1]], colnames(d_expr), sep="_")
                return(d_expr);
        })
	out <- clust_expr[[1]];
	for (i in 2:length(clust_expr)) {
		c_names <- c(colnames(out), colnames(clust_expr[[i]]))
		out <- cbind(out, clust_expr[[i]]);
		colnames(out) <- c_names
	}
        return(out)
}

