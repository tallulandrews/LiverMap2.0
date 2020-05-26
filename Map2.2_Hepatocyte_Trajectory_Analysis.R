dat <- readRDS("Hep-Traj_Subcluster_analysis.rds")
require("Seurat")
require("destiny")

set.seed(2820)


#transition <- dat[,dat@meta.data$Use_clusters %in% c(0,2,3)]
transition <- dat;
transition <- FindVariableFeatures(transition, nfeatures=1000)
is.Central <- transition@meta.data$consistent_labs=="PericentralHep"
is.Portal <- transition@meta.data$consistent_labs=="PeriportalHep"

v_genes <- VariableFeatures(transition)
transition <- RunPCA(transition, features = VariableFeatures(object = transition))
dat@meta.data$seurat_clusters <- dat@meta.data$Use_clusters


plot_props <- function(pseudotime, is.Central, is.Portal) {
	bin_width = diff(range(pseudotime))/50
	prop_C = sapply(pseudotime, function(x) {
		mean(is.Central[pseudotime < x+bin_width/2 & 
		pseudotime > x-bin_width/2])})
	prop_P = sapply(pseudotime, function(x) {
		mean(is.Portal[pseudotime < x+bin_width/2 & 
		pseudotime > x-bin_width/2])})
	plot(pseudotime[order(pseudotime)], 100*prop_C[order(pseudotime)], 
		type="l", lwd=2.5, xlab="", ylab="% of cells",
		ylim=c(0,1), col="plum")	
	lines(pseudtotime[order(pseudotime)], 100*prop_P[order(pseudotime)], lwd=2.5, col="dodgerblue")
	invisible(cbind(prop_C, prop_P))
}

png("Hepatocyte_transition_pca1.png", width=6, height=6, units="in", res=50)
DimPlot(transition, reduction="pca", group.by="scmap_anno2")
pseudotime_pca <- -1*transition@reductions$pca@cell.embeddings[,2]
dev.off()
png("Hepatocyte_transition_pca2.png", width=6, height=6, units="in", res=50)
plot_props(pseudotime_pca, is.Central, is.Portal)
title(xlab="Pseudotime (PC2)")
dev.off()

transition@meta.data$is.Central <- is.Central
transition@meta.data$is.Portal <- is.Portal
transition@meta.data$pca_pseudotime <- pseudotime_pca

mat <- transition@assays$RNA@data[rownames(transition) %in% v_genes,]
set.seed(30192)
dm <- DiffusionMap(t(as.matrix(mat)), n_pcs=20, n_eigs=3)
pseudotime_dm = dm@eigenvectors[,2]

png("Hepatocyte_transition_dm1.png", width=6, height=6, units="in", res=50)
plot(dm@eigenvectors[,1], dm@eigenvectors[,2], col=c("red", "black")[is.Central+1]) 
dev.off()

png("Hepatocyte_transition_dm2.png", width=6, height=6, units="in", res=50)
plot_props(pseudotime_dm, is.Central)
dev.off();

transition@meta.data$dm_pseudotime <- pseudotime_dm
saveRDS(transition, "Hepatocyte_transition.rds")



# Monocle DE

require("monocle")
#The sm.ns function states that Monocle should fit a natural spline through 
#the expression values to help it describe the changes in expression as a 
#function of progress.

#Monocle assigns each cell a "pseudotime" value, which records its progress 
#through the process in the experiment. The model can test against changes 
#as a function of this value. Monocle uses the VGAM package to model a gene's 
#expression level as a smooth, nonlinear function of pseudotime.

exprs <- transition@assays$RNA@counts
pd <- AnnotatedDataFrame(transition@meta.data)
pd$Pseudotime <- pseudotime_pca
fd <- data.frame(gene_short_name=rownames(transition), 
		is.feature=rownames(transition) %in% v_genes)
rownames(fd) <- rownames(transition)
fd <- AnnotatedDataFrame(fd)

mono_obj <- newCellDataSet(exprs,
                phenoData = pd,
                featureData = fd,
                expressionFamily=negbinomial.size())

mono_obj <- estimateSizeFactors(mono_obj)
mono_obj <- estimateDispersions(mono_obj)

diff_test_res <- differentialGeneTest(mono_obj,
fullModelFormulaStr = "~sm.ns(Pseudotime) + donor", reducedModelFormulaStr = "~donor")

res <- diff_test_res[order(diff_test_res$qval),]

write.table(res, file="Hepatocyte_Trajectory_monocle_de.csv", sep=",")



my_rowMeans <- function(x) {
        if (!is.null(ncol(x))) {
                if (ncol(x) > 1) {
                        return(Matrix::rowMeans(x))
                }
        }
        return(x);
}

bin_width <- range(pseudotime_pca)/30

smoothed_expr <- sapply(pseudotime_pca, function(x) {
	my_rowMeans(transition@assays$RNA@data[,
		pseudotime_pca < x+bin_width/2 & 
		pseudotime_pca > x-bin_width/2])}
	)


smoothed_expr <- smoothed_expr[,order(pseudotime_pca)]
n_columns <- 100
n_each <- ceiling(ncol(smoothed_expr)/n_columns)
thing <- rep(1:100, each=n_each)
thing <- thing[1:ncol(smoothed_expr)]
thing <- split(seq(ncol(smoothed_expr)),factor(thing))
heatmap_dat <- sapply(thing, function(group) my_rowMeans(smoothed_expr[,group]))

heatmap_dat <- heatmap_dat[match(rownames(res), rownames(heatmap_dat)),]
res$dir <- sign(heatmap_dat[,100]-heatmap_dat[,1])
res$score <- log(res$qval)*res$dir

write.table(res, file="Hepatocyte_Trajectory_monocle_de.csv", sep=",")

sig <- res[res$qval < 0.01,]
sig <- sig[order(sig$score),]
to_plot <- heatmap_dat[match(rownames(sig), rownames(heatmap_dat)),]

to_plot <- t(apply(to_plot, 1, scale))
lim <- max(abs(min(to_plot)), max(to_plot))

png("Hepatocyte_transistion_monocle_de_heatmap.png", width=6, height=6, units="in", res=300)
heatmap(to_plot, Colv=NA, Rowv=NA, scale="none",
col=colorRampPalette(c("magenta", "black", "yellow"))(20), 
breaks=seq(from=-lim, to=lim, length=21), xaxt="n", xlab="Pseudotime")
dev.off()

require("gprofiler2")
res_early <- gprofiler2::gost(rownames(sig[sig$dir==1,]), ordered_query=T, correction_method="fdr", sources=c("GP:BP", "KEGG", "REAC"))$result
res_late <- gprofiler2::gost(rev(rownames(sig[sig$dir==-1,])), ordered_query=T, correction_method="fdr", sources=c("GP:BP", "KEGG", "REAC"))$result

write.table(as.matrix(res_early), file="Hepatocyte_transition_monocle_enrichment_early.csv", sep=",")

write.table(as.matrix(res_late), "Hepatocyte_transition_monocle_enrichment_late.csv", sep=",")
