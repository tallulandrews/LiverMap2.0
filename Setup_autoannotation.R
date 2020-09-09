require(SingleCellExperiment)
require(scmap)
require(SCINA)

source("/cluster/home/tandrews/scripts/LiverMap2.0/My_R_Scripts.R")
#auto_anno_dir <- "/home/gelder/MacParlandLabData/human/HumanLiver1.0/";
auto_anno_dir <- "/cluster/projects/macparland/TA/AutoAnnotation"

map1_ref <- readRDS(paste(auto_anno_dir,"scmap_reference.rds", sep="/"))

map1_markers <- read.table(paste(auto_anno_dir, "my_marker_genes.txt", sep="/"), header=T, stringsAsFactors=FALSE)

types <- unique(map1_markers[,2])
types <- types[types != "None"]
map1_markers_list <- list();
for (t in types) {
	map1_markers_list[[t]]<-map1_markers[map1_markers[,2] == t,1]
}

run_SCINA <- function(mat, marker_list=map1_markers_list) {
	# SLOW-ish
	c_rate <- max(0.99, 1-10/ncol(mat));
	results <- SCINA::SCINA(mat, marker_list, max_iter=100, convergence_n=4, convergence_rate=c_rate, sensitivity_cutoff=0.05, allow_unknown=1)
	
}

#require(CellTypeProfiles)
my_markers <- function(mat) {
        on_off <- matrix(0, ncol=ncol(mat), nrow=nrow(mat));
        my_split_max_gap <- function(x) {
                x <- sort(x)
                jumps <- diff(x);
                br_pt <- which(jumps == max(jumps))
		if (length(br_pt) > 1) {
			br_pt <- br_pt[length(br_pt)]
		}
                return(c(x[br_pt], max(jumps)));
        }
        thresh <- apply(mat, 1, my_split_max_gap);
        on_off <- t(sapply(1:ncol(thresh), function(i) {mat[i,] > thresh[1,i]}))
	rownames(on_off) <- rownames(mat);
	colnames(thresh) <- rownames(mat);
        return(list(score=thresh[2,], on_off=on_off));
}


run_scmap_seurat <- function(myseur, scmap_ref=map1_ref, return_sce=FALSE) {
	# make sure raw counts and lognormalized matrices match
	myseur@assays$RNA@counts <- myseur@assays$RNA@counts[match(rownames(myseur@assays$RNA@data), rownames(myseur@assays$RNA@counts)),]
	# create SCE
	mysce <- SingleCellExperiment(assays=list(counts=myseur@assays$RNA@counts, logcounts=myseur@assays$RNA@data), colData=myseur@meta.data)


#	mysce <- as.SingleCellExperiment(myseur)
	rowData(mysce)$feature_symbol=rownames(mysce);
	mysce <- mysce[!grepl("^MT-", rownames(mysce)),] #remove MT genes.
	mysce <- mysce[!grepl("^RPS-", rownames(mysce)),] #remove Ribo genes.
	mysce <- mysce[!grepl("^RPL-", rownames(mysce)),] #remove Ribo genes.
	
#	for (i in seq(from=1, to=ncol(mysce), by=1000)) {
#		tmp <- mysce[,seq(from=i, to=min(ncol(mysce), i+1000))]
	# scmap_cluster
	scmap_annotation <- scmapCluster( projection = mysce,
		index_list = list(lm1 = metadata(map1_ref)$scmap_cluster_index), 
		threshold=0.1)
	mysce$scmap_id <- as.vector(scmap_annotation$scmap_cluster_labs)
	mysce$scmap_score <- as.vector(scmap_annotation$scmap_cluster_siml)


	myseur@meta.data$scmap_cluster_anno <- mysce$scmap_id
	myseur@meta.data$scmap_cluster_score <- mysce$scmap_score

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

cell_anno_to_cluster_anno <- function(cellanno, clusterids) {
	tab <- table(cellanno, clusterids)
	clusterlab <-  apply(tab, 2, function(x){
		out <- rownames(tab)[which(x==max(x))]
		if (length(out) > 1) {out <- "ambiguous"}
		return(out);
		})
	return(data.frame(cluster=colnames(tab), lab=clusterlab));
}

Use_markers_for_anno <- function(mat, clusters, ref_markers=map1_markers) {
	# get average expression by cluster
	cluster_means <- group_rowmeans(mat, clusters);
	# get % detect by cluster
	tmp <- mat;
	tmp[tmp>0] <- 1;
	cluster_detect <- group_rowmeans(tmp, clusters);

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
		ps[ks==0] <- 1;
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
