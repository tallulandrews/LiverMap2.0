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
my_colMeans <- function(x) {
        if (!is.null(nrow(x))) {
                if (nrow(x) > 1) {
                        return(Matrix::colMeans(x))
                }
        }
        return(x);
}
my_colSums <- function(x) {
        if (!is.null(nrow(x))) {
                if (nrow(x) > 1) {
                        return(Matrix::colSums(x))
                }
        }
        return(x);
}

group_rowmeans <- function(MAT, group_labs, type=c("mean","sum")) {
        d <- split(seq(ncol(MAT)), group_labs);
	if (type[1] == "mean") {
	        mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
	} else {
	        mus <- sapply(d, function(group) my_rowSums(MAT[,group]))
	} 
        return(mus);
}
group_colmeans <- function(MAT, group_labs, type=c("mean", "sum")) {
        d <- split(seq(nrow(MAT)), group_labs);
	if (type[1] == "mean") {
        	mus <- sapply(d, function(group) my_colMeans(MAT[group,]))
	} else {
        	mus <- sapply(d, function(group) my_colSums(MAT[group,]))
	}
        return(mus);
}

# Average expression per cluster per donor, optional- weight by donor freqs
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

# Table of mean expression of cells from each donor in each cluster - for DE
get_pseudobulk <- function(mat, clusters, donors) {
        c <- split(seq(ncol(mat)), clusters);
        donor_freqs <- table(donors)/length(donors)
        # avg expression per donor in this cluster
        clust_expr <- sapply(c, function(clust) {
                d_expr <- group_rowmeans(mat[,clust], donors[clust], type="sum");
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


# means for heatmap?
get_pseudobulk_means <- function(mat, clusters, donors) {
        c <- split(seq(ncol(mat)), clusters);
        donor_freqs <- table(donors)/length(donors)
        # avg expression per donor in this cluster
        clust_expr <- sapply(c, function(clust) {
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


