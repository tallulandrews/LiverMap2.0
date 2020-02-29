tidy_Seurat_output <- function(out) {
	out$score <- out$avg_logFC*(out$pct.1-out$pct.2)
	out <- out[order(out$score, decreasing=T),]
	return(out)
}

#args <- commandArgs(trailingOnly=TRUE)
# file
# clustername
# edger_file
# cluster_col_name

require(Seurat)
obj <- readRDS(args[1])

if (!"Use_clusters" %in% colnames(obj@meta.data)) {
	obj@meta.data$Use_clusters <- obj@meta.data[,args[4]]
}
 
out_seur_wilcox <- FindMarkers(obj, ident.1=args[2], group.by="Use_clusters", logfc.threshold=-Inf, min.pct=0.1, min.diff.pct = -Inf, test.use="wilcox")

require(MAST)
out_seur_MAST <- FindMarkers(obj, ident.1=args[2], group.by="Use_clusters", logfc.threshold=-Inf, min.pct=0.1, min.diff.pct = -Inf, test.use="MAST", latent.vars=obj$donor)

## MAST? ## - Slow and can't figure out how to do the contrasts.
tidy_mast_output <- function(output, this_contrast) {
	res <- output$datatable[output$datatable$contrast==this_contrast & output$datatable$component %in% c("C","D"),]
	res$fdr <- p.adjust(unlist(res[,4]), method="fdr")

	con <- res[res$component =="C",]
	dis <- res[res$component =="D",]

	sig <- sign(con$coef) == sign(dis$coef)
	out <- cbind(con[sig,], dis[sig,])
	out$overall <- out[,7]+out[,16]
	out <- out[order(out$overall, decreasing=T),]
	return(out);
}

#require(MAST)
#sca <- FromMatrix(as.matrix(obj@assays$RNA@data))
#ngenes <- colSums(obj@assays$RNA@counts > 0)
#ncounts <- colSums(obj@assays$RNA@counts)
#donor <- obj@meta.data$donor
#cluster <- obj@meta.data$Use_clusters
#mod <- zlm(~cluster + donor + ngenes, sca)

#set.seed(3891)

#res1 <- summary(mod, doLRT=Hypothesis(this_contrast))

#out <- tidy_mast_output(res1, this_contrast)


require("edgeR")
set.seed(9175)

edger_obj <- DGEList(obj@assays$RNA@counts, samples=obj@meta.data, group=obj@meta.data$Use_clusters)
design <- model.matrix(~Use_clusters+donor+nFeature_RNA, data=edger_obj$samples)

coef_names <- colnames(design);

if (!file.exsists(args[3])) {
	edger_obj <- estimateDisp(edger_obj, design)
	fit <- glmQLFit(edger_obj, design)
	saveRDS(fit, args[3]);
} else {
	fit <- readRDS(args[3])
}

contrast_vec <-rep(0, length(coef_names));
this_cluster <- paste("Use_clusters", arg[2], sep="");
contrast_vec[grepl("Use_clusters", coef_names)] <- -1;
contrast_vec[1] <- -1
if (this_cluster %in% coef_names) {
	contrast_vec[coef_names==this_cluster]<- 1;
} else {
	contrast_vec[1] <- 1;
}

this_vs_all <- glmQLFTest(fit, contrast=contrast_vec)
out <- topTags(this_vs_all, nrow(obj));
