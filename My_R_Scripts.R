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

# Average expression per cluster per donor, optional- weight by donor freqs
group_rowmeans <- function(MAT, group_labs) {
        d <- split(seq(ncol(MAT)), group_labs);
        mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
        return(mus);
}

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


