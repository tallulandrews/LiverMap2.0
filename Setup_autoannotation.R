require(SingleCellExperiment)
require(scmap)

map1_ref <- readRDS("/home/gelder/MacParlandLabData/human/HumanLiver1.0/scmap_reference.rds")

map1_markers <- read.table("/home/gelder/MacParlandLabData/human/HumanLiver1.0/my_marker_genes.txt", header=T, stringsAsFactors=FALSE)

require(CellTypeProfiles)
my_markers <- function(mat) {
        on_off <- matrix(0, ncol=ncol(mat), nrow=nrow(mat));
        my_split_max_gap <- function(x) {
                x <- sort(x)
                jumps <- diff(x);
                br_pt <- which(jumps == max(jumps))
                return(c(x[br_pt], max(jumps)));
        }
        thresh <- apply(mat, 1, my_split_max_gap);
        on_off <- t(sapply(1:ncol(thresh), function(i) {mat[i,] > thresh[1,i]}))
	rownames(on_off) <- rownames(mat);
	colnames(thresh) <- rownames(mat);
        return(list(score=thresh[2,], on_off=on_off));
}


run_scmap_seurat <- function(myseur, scmap_ref=map1_ref, return_sce=FALSE) {
	myseur@assays$RNA@counts <- myseur@assays$RNA@counts[match(rownames(myseur@assays$RNA@data), rownames(myseur@assays$RNA@counts)),]
	mysce <- SingleCellExperiment(assays=list(counts=myseur@assays$RNA@counts, logcounts=myseur@assays$RNA@data), colData=myseur@meta.data)


#	mysce <- as.SingleCellExperiment(myseur)
	rowData(mysce)$feature_symbol=rownames(mysce);
	mysce <- mysce[!grepl("^MT-", rownames(mysce)),] #remove MT genes.

	# scmap_cluster
	scmap_annotation <- scmapCluster( projection = mysce,
		index_list = list(lm1 = metadata(map1_ref)$scmap_cluster_index), 
		threshold=0.1)
	mysce$scmap_id <- scmap_annotation$scmap_cluster_labs
	mysce$scmap_score <- scmap_annotation$scmap_cluster_siml
	mysce$scmap_id <-  scmap_annotation$scmap_cluster_labs
	mysce$scmap_score <-  scmap_annotation$scmap_cluster_siml


	myseur@meta.data$scmap_cluster_anno <- data.frame(id=mysce$scmap_id, similarity=mysce$scmap_score);

	# scmap_cell
	scmap_cell_res <- scmapCell(mysce, index_list=list(lm1=metadata(map1_ref)$scmap_cell_index));
	getmode <- function(v) {
	   uniqv <- unique(v)
	   uniqv[which.max(tabulate(match(v, uniqv)))]
	}

	cell_anno <-  apply(scmap_cell_res$lm1$cells,2,
        	function(x){
                	anns = map1_ref$cell_type1[x]; 
                	assign = getmode(anns); 
                	if(length(assign) > 1) {return("ambiguous")}
                	else{return(assign)}
        })


	myseur@meta.data$scmap_cell_anno <- cell_anno;
	mysce$scmap_cell_anno <- cell_anno;
	if (return_sce) {
		return(list(seurat=myseur, sce=mysce));
	} else {
		return(myseur)
	}
}

cell_anno_to_cluster_anno <- function(cellids, clusterids) {
	tab <- table(scmap_annotation$scmap_cluster_labs, mysce$seurat_clusters)
	clusterlab <-  apply(tab, 2, function(x){rownames(tab)[which(x==max(x))]})
	return(data.frame(cluster=colnames(tab), lab=clusterlab));
}




Use_markers_for_anno <- function(mat, clusters, ref_markers=map1_markers) {
	# get average expression by cluster
	cluster_means <- my_row_mean_aggregate(mat, clusters);
	# get % detect by cluster
	tmp <- mat;
	tmp[tmp>0] <-1;
	cluster_detect <- my_row_mean_aggregate(tmp, clusters);

	# get markers based on the maximum jump between clusters.
	mark_mean <- my_markers(cluster_means);
	mark_detect <- my_markers(cluster_detect);

	# good marker = change of 0.3 in mean expression or 
	# change of 0.1 in proportion of cells expressing the marker.
	# this is the same as I used for the reference markers.
	good <- mark_mean$score > 0.3 | mark_detect$score > 0.1;
	# mark one cluster or many?
	unique <- rowSums(mark_mean$on_off) == 1 & 
		  rowSums(mark_detect$on_off) == 1
	# detection rate & mean agree perfectly?
	agree <- apply((mark_mean$on_off+mark_detect$on_off), 1,
			function(x){sum(x==1)==0})

	# Those clusters where both methods agree the marker in "on"
	# and only positive markers ('on' in less than half the clusters)
	tab <- mark_mean$on_off & mark_detect$on_off
	tab <- tab[good & rowSums(tab) < ncol(tab)/2,]

	# cross reference with the reference markers
	ref <- ref_markers[ref_markers[,2] != "None",]
	ref <- ref[ref[,1] %in% rownames(tab),]
	ref[,2] <- factor(ref[,2])
	tab <- tab[match(ref[,1],rownames(tab)),]
	
	# use hypergeometric test/fisher's exact test
	# to determine significant enrichments for a set of
	# reference markers.
	result <- vector();
	c_lab <- vector();
	for (lab in unique(ref[,2])) {
		n_lab <- sum(ref[,2] == lab);
		if (n_lab < 3) {next;}
		N <- nrow(ref);
		xs <- colSums(tab[ref[,2] == lab,])
		ks <- colSums(tab);
		ps <- sapply(1:length(ks), function(i){phyper(xs[i], n_lab, N-n_lab, ks[i], lower.tail=FALSE)});
		result <- rbind(result, ps);
		c_lab <- c(c_lab, lab);
	}
	colnames(result) <- colnames(tab);
	rownames(result) <- c_lab;
	
	# assign each novel cluster to its best reference cluster.
	best <- apply(result,2,function(x){
		if (sum(x==min(x))==1) {
			return(rownames(result)[which(x==min(x))])
		} else {
			return("ambiguous")
		}
	})

	return(list(ps=result, 
		cluster_assign=data.frame(cluster=names(best), label=best), 
		cell_assign=best[match(clusters,names(best))]))
}


