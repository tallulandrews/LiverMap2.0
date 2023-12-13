source("../../scripts/LiverMap2.0/My_R_Scripts.R")


files <- c("AntiB_harmony_Subcluster.rds", 
		"Cholangiocyte_harmony_Subcluster.rds", 
		"Endo_harmony_Subcluster.rds", 
		"Hepatocyte1_harmony_Subcluster.rds", 
		"Hepatocyte2_harmony_Subcluster.rds", 
		"Macrophage_harmony_Subcluster.rds", 
		"NKT_harmony_Subcluster.rds", 
		"Stellate_harmony_Subcluster.rds")



get_metadata_outcome <- function(obj){
	tmp <- read.delim("../../Caudate_recip_data_Dec 3_20.csv", sep=",")
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

	pseudobulk <- get_pseudobulk_means(obj@assays$RNA@scale.data, obj@meta.data[,cluster_column], obj@meta.data$sample)
	obj <- get_metadata_outcome(obj)
	Ns <- table(obj@meta.data$sample, obj@meta.data[,cluster_column]) 

	## The sample sizes breaking down by sample and subcluster are too small

	pseudobulk <- group_rowmeans(obj@assays$RNA@scale.data, obj@meta.data$sample, type="mean")
	pseudobulk2 <- group_rowmeans(obj@assays$RNA@data, obj@meta.data$sample, type="mean")


	Ns <- table(obj@meta.data$sample) 

	keep <- Ns >= 10;

	# trans.rejected
	label <- sample2meta(obj, meta_col="trans.rejected")[keep]
	expr_mat <- pseudobulk[,keep]
	#expr_mat <- expr_mat[Matrix::rowSums(obj@assay$RNA@counts)

	de <- t(apply(expr_mat, 1, function(x) {
		res = t.test(x[label=="Y"], x[label=="N"])
		return(c(res$estimate, res$p.value))
		}))
	de <- cbind(de, p.adjust(de[,3], "fdr"))
	colnames(de) <- c("Y", "N", "pval", "qval")

	# Sex
	label <- sample2meta(obj, meta_col="donor.sex")[keep]
	expr_mat <- pseudobulk[,keep]
	
	de <- t(apply(expr_mat, 1, function(x) {
		res = t.test(x[label=="F"], x[label=="M"])
		return(c(res$estimate, res$p.value))
		}))
	de <- cbind(de, p.adjust(de[,3], "fdr"))
	colnames(de) <- c("Y", "N", "pval", "qval")

## Problem: Nothing significant!

###### DE at single cell level #####
expr_mat <- obj@assays$RNA@scale.data
label <- obj@meta.data$trans.rejected

de <- t(apply(expr_mat, 1, function(x) {
		res = t.test(x[label=="Y"], x[label=="N"])
		return(c(res$estimate, res$p.value))
		}))
de <- cbind(de, p.adjust(de[,3], "fdr"))
colnames(de) <- c("Y", "N", "pval", "qval")

write.table(rownames(de)[ de[,4] < 0.05 & de[,1] > de[,2] ], file="tmp_up.txt", row.names=F, col.names=F, quote=F)
write.table(rownames(de)[ de[,4] < 0.05 & de[,1] < de[,2] ], file="tmp_dn.txt", row.names=F, col.names=F, quote=F)


#de_score <- round(de[,1]-de[,2], digits=1)
de_score <- -log(de[,4])*sign(de[,1]-de[,2])