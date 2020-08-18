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
        clust_expr <- sapply(c, function(clust) {
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

