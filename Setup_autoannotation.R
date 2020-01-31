require(SingleCellExperiment)
require(scmap)

map1_ref <- readRDS("/home/gelder/MacParlandLabData/human/HumanLiver1.0/scmap_reference.rds")

map1_markers <- read.table("/home/gelder/MacParlandLabData/human/HumanLiver1.0/my_marker_genes.txt", header=T)

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
        return(list(score=thresh[2,], on_off=on_off));
}


run_scmap_seurat <- function(myseur, scmap_ref=map1_ref, return_sce=FALSE) {
	mysce <- as.SingleCellExperiment(myseur)
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

Use_markers_for_anno <- function(mat, clusters, ref_markers) {


}


