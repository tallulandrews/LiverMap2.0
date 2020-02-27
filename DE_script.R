tidy_Seurat_output <- function(out) {
	out$score <- out$avg_logFC*(out$pct.1-out$pct.2)
	out <- out[order(out$score, decreasing=T),]
	return(out)
}


require(Seurat)
obj <- readRDS()
out <- FindMarkers(obj, ident.1="0", group.by="Use_clusters", logfc.threshold=-Inf, min.pct=0.1, min.diff.pct = -Inf, test.use="wilcox")

require(MAST)
out <- FindMarkers(obj, ident.1="0", group.by="Use_clusters", logfc.threshold=-Inf, min.pct=0.1, min.diff.pct = -Inf, test.use="MAST", latent.vars=obj$donor)


require(MAST)
sca <- FromMatrix(as.matrix(obj@assays$RNA@data))
ngenes <- colSums(obj@assays$RNA@counts > 0)
ncounts <- colSums(obj@assays$RNA@counts)
donor <- obj@meta.data$donor
cluster <- obj@meta.data$Use_clusters
mod <- zlm(~cluster + donor + ngenes +0, sca)

set.seed(3891)

res1 <- summary(mod, doLRT=Hypothesis(this_contrast))

out <- tidy_mast_output(res1, this_contrast)

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
