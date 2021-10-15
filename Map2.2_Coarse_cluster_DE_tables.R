source("../../scripts/LiverMap2.0/My_R_Scripts.R")


get_metadata_outcome <- function(obj){
	tmp <- read.delim("/cluster/home/tandrews/scripts/LiverMap2.0/Caudate_recip_data_Dec_3_20.csv", sep=",")
	reject <- tmp[ match(obj@meta.data$donor, tmp[,1]), "Post.LT.Rejection"]
	reject <- factor(reject, levels=c("N", "?", "Y"))
	obj@meta.data$trans.rejected <- reject
	return(obj)
}

sample2meta <- function(obj, meta_col="trans.rejected") {
	tab <- table(obj@meta.data$sample, obj@meta.data[,meta_col])
	res = apply(tab, 1, function(x) {
			out = colnames(tab)[which(x == max(x))];
			if (length(out) > 1) {return(NA)} else {return(out)}
			})
	return(res)
}

## Combining DE results from Seurat

#Fisher's Combined Probabliltiy test https://en.wikipedia.org/wiki/Fisher%27s_method
fisher_combined <- function(p_vals) {
	p_vals[is.na(p_vals)] <- 1;
	chi_sq <- sum(log(p_vals))*-2
	new_p <- pchisq(chi_sq, df=2*length(p_vals), lower.tail=FALSE, log.p=FALSE)
	return(new_p)
}

fix_nas <- function(out) {
	nas <- is.na(out[,1]);
	out[nas,1] <- 1;
	out[nas,2] <- 0;
	out[nas,3] <- 0;
	out[nas,4] <- 0;
	out[nas,5] <- 1;
	return(out);
}

combine_tests_Seurat <- function(out1, out2, n1, n2) {
	all_genes <- sort(unique(c(rownames(out1), rownames(out2))))
	out1 <- fix_nas(out1[match(all_genes, rownames(out1)),])
	out2 <- fix_nas(out2[match(all_genes, rownames(out2)),])
	
	l2fc <- (out1[,"avg_logFC"]*n1 + out2[,"avg_logFC"]*n2)/(n1+n2)
	perc1 <- (out1[,"pct.1"]*n1 + out2[,"pct.1"]*n2)/(n1+n2)
	perc2 <- (out1[,"pct.2"]*n1 + out2[,"pct.2"]*n2)/(n1+n2)
	p_vals <- cbind(out1[,"p_val_adj"], out2[,"p_val_adj"])
	if (abs(n1-n2)/min(n1, n2) < 10) {
		new_pvals <- apply(p_vals, 1, fisher_combined)
	} else {
		if (n1 < n2) { new_pvals <- out2[,"p_val_adj"] }
		if (n2 < n1) { new_pvals <- out1[,"p_val_adj"] }
	}
	new_out <- cbind(l2fc, perc1, perc2, new_pvals);
	rownames(new_out) <- all_genes;
	colnames(new_out) <- c("avg_logFC", "pct.1", "pct.2", "p_val_adj");
	new_out <- new_out[order(new_out[,"p_val_adj"]),]
	return(new_out)
}

combine_tests_edgeR <- function(out1, out2, n1, n2) {
	all_genes <- sort(unique(c(rownames(out1), rownames(out2))))
	out1 <- fix_nas(out1[match(all_genes, rownames(out1)),])
	out2 <- fix_nas(out2[match(all_genes, rownames(out2)),])
	
	l2fc <- (out1[,"avg_logFC"]*n1 + out2[,"avg_logFC"]*n2)/(n1+n2)
	perc1 <- (out1[,"pct.1"]*n1 + out2[,"pct.1"]*n2)/(n1+n2)
	perc2 <- (out1[,"pct.2"]*n1 + out2[,"pct.2"]*n2)/(n1+n2)
	p_vals <- cbind(out1[,"p_val_adj"], out2[,"p_val_adj"])
	if (abs(n1-n2)/min(n1, n2) < 10) {
		new_pvals <- apply(p_vals, 1, fisher_combined)
	} else {
		if (n1 < n2) { new_pvals <- out2[,"p_val_adj"] }
		if (n2 < n1) { new_pvals <- out1[,"p_val_adj"] }
	}
	new_out <- cbind(l2fc, perc1, perc2, new_pvals);
	rownames(new_out) <- all_genes;
	colnames(new_out) <- c("avg_logFC", "pct.1", "pct.2", "p_val_adj");
	new_out <- new_out[order(new_out[,"p_val_adj"]),]
	return(new_out)
}


edgeR_twogroup_pseudobulk_de <- function(pseudobulk, labels, group1=labels[1], group2=labels[2]) {
		require(edgeR);
		keep <- labels %in% c(group1, group2);
		pseudobulk <- pseudobulk[,keep]		
		labels <- labels[keep]
		pseudobulk <- pseudobulk[rowSums(pseudobulk > 0) > 3,]

		edger <- DGEList(counts=pseudobulk, group=labels);
		edger <- estimateDisp(edger)
		et <- exactTest(edger)
		out <- topTags(et, nrow(pseudobulk));
		return(out);
}


require(Seurat)
set.seed(28210)
tag <- "Map2.2_EmptyOnly_Coarse"
obj <- readRDS("Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno_dimreduce.rds")
obj <- get_metadata_outcome(obj)
obj@meta.data$donor_sex <- factor(obj@meta.data$donor_sex)
cluster_column <- "Coarse_clusters"

dat_3pr <- obj[,obj@meta.data$assay_type == "3pr"]
dat_5pr <- obj[,obj@meta.data$assay_type == "5pr"]

# Cluster DE
for (c in unique(obj@meta.data[,cluster_column]) ) {
	if (sum(dat_3pr@meta.data[,cluster_column] == c) > 10 & 
		sum(dat_5pr@meta.data[,cluster_column] == c) > 10) {
		
		out_3pr <- FindMarkers(dat_3pr, group.by=cluster_column, ident.1=c, logfc.threshold=0)
		out_5pr <- FindMarkers(dat_5pr, group.by=cluster_column, ident.1=c, logfc.threshold=0)
	
		out <- combine_tests_Seurat(out_3pr, out_5pr, 
			sum(dat_3pr@meta.data[,cluster_column] == c), 
			sum(dat_5pr@meta.data[,cluster_column] == c))
	} else {
		if (sum(dat_3pr@meta.data[,cluster_column] == c) > 10) {
			out_3pr <- FindMarkers(dat_3pr, group.by=cluster_column, ident.1=c, logfc.threshold=0)
			out <- out_3pr[,-1]
		} else if (sum(dat_5pr@meta.data[,cluster_column] == c) > 10) {
			out_5pr <- FindMarkers(dat_5pr, group.by=cluster_column, ident.1=c, logfc.threshold=0)
			out <- out_5pr[,-1]
		}
	}

	write.table( out, file=paste( "Cluster_DE_Tables", tag, c, "combined_DE_table.csv", sep="_"), 
			sep=",", row.names=T, col.names=T, quote=T )
}
	
## Pseudobulk - edgeR
	
pseudobulk_3pr <- group_rowmeans(dat_3pr@assays$RNA@counts, dat_3pr@meta.data$sample, type="sum")
pseudobulk_5pr <- group_rowmeans(dat_5pr@assays$RNA@counts, dat_5pr@meta.data$sample, type="sum")

keep_3pr <- table(dat_3pr@meta.data$sample) >= 10
keep_5pr <- table(dat_5pr@meta.data$sample) >= 10

labels_3pr <- sample2meta(dat_3pr, meta_col="trans.rejected")[keep_3pr]
labels_5pr <- sample2meta(dat_5pr, meta_col="trans.rejected")[keep_5pr]
edgeR_3pr <- edgeR_twogroup_pseudobulk_de(pseudobulk_3pr[,keep_3pr], labels_3pr, group1="Y", group2="N")
edgeR_5pr <- edgeR_twogroup_pseudobulk_de(pseudobulk_5pr[,keep_5pr], labels_5pr, group1="Y", group2="N")

## Problem: Nothing significant! & extremely inconsistent btw 3' and 5'

## Cell-level - more significant & much more consistent
out_3pr <- FindMarkers(dat_3pr, group.by="trans.rejected", ident.1="Y", ident.2="N", logfc.threshold=0)
out_5pr <- FindMarkers(dat_5pr, group.by="trans.rejected", ident.1="Y", ident.2="N", logfc.threshold=0)
out <- combine_tests_Seurat(out_3pr, out_5pr, 
			sum(dat_3pr@meta.data[,"trans.rejected"] %in% c("Y", "N")), 
			sum(dat_5pr@meta.data[,"trans.rejected"] %in% c("Y", "N")))
write.table( out, file=paste( "SubCluster_DE_Tables", tag, c, "combined_DE_TransReject.csv", sep="_"),
		 sep=",", row.names=T, col.names=T, quote=T )

out_3pr <- FindMarkers(dat_3pr, group.by="donor_sex", ident.1="F", ident.2="M", logfc.threshold=0)
out_5pr <- FindMarkers(dat_5pr, group.by="donor_sex", ident.1="F", ident.2="M", logfc.threshold=0)
out <- combine_tests_Seurat(out_3pr, out_5pr, 
			sum(dat_3pr@meta.data[,"donor_sex"] %in% c("F", "M")), 
			sum(dat_5pr@meta.data[,"donor_sex"] %in% c("F", "M")))
write.table( out, file=paste( "SubCluster_DE_Tables", tag, c, "combined_DE_Sex.csv", sep="_"), 
		sep=",", row.names=T, col.names=T, quote=T )

out_3pr <- FindMarkers(dat_3pr, group.by="donor_age_group", ident.1="young", ident.2="elderly", logfc.threshold=0)
write.table( out_3pr[,-1], file=paste( "Cluster_DE_Tables", tag, c, "3pr_DE_Age.csv", sep="_"), 
		sep=",", row.names=T, col.names=T, quote=T )

out_3pr <- FindMarkers(dat_3pr, group.by="donor_bmi_group", ident.1="normal", ident.2="very obese", logfc.threshold=0)
write.table( out_3pr[,-1], file=paste( "SubCluster_DE_Tables", tag, c, "3pr_DE_Age.csv", sep="_"), 
		sep=",", row.names=T, col.names=T, quote=T )

