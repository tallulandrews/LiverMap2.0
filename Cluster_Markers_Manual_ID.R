set.seed(1189)

int.obj <- readRDS("integration_harmony_plus_analysis.rds")


# Get auto-annotation for fine & coarse clusters

cluster_assign <- function(x) {
	freqs <- table(x)/length(x);
	return(freqs[freqs==max(freqs)])
}

general_cell_type_ids <- as.character(int.obj@meta.data$scmap_anno2)
general_cell_type_ids[general_cell_type_ids %in% c("CD3abTcells","gdTcells2")] <- "Tcells"
general_cell_type_ids[general_cell_type_ids %in% c("MatureBcells","AntibodysecretingBcells")] <- "Bcells"
general_cell_type_ids[general_cell_type_ids %in% c("PeriportalLSECs","CentralvenousLSECs")] <- "LSECs"
general_cell_type_ids[general_cell_type_ids %in% c("interzonalHep","PeriportalHep", "UnidentifiedHep", "PericentralHep")] <- "Hepatocyte"
general_cell_type_ids[general_cell_type_ids %in% c("Non-inflammatoryMacrophages", "InflamatoryMacrophages")] <- "Macrophage"

assign_specific <- sapply(split(int.obj@meta.data$scmap_anno2, int.obj@meta.data$Fine_clusters), cluster_assign)
assign_general <- sapply(split(general_cell_type_ids, int.obj@meta.data$Fine_clusters), cluster_assign)


Cluster_guesses <- paste(0:(length(assign_general)-1),".",rep("Unknown", length(assign_general)), sep="");
Cluster_guesses[assign_general > 0.5] <- names(assign_general)[assign_general>0.5]
Cluster_guesses[assign_specific > 0.5] <- names(assign_specific)[assign_specific>0.5]

## Figures for Talk ##
agg_coord_by_cluster <- function(coords, clusters) {
	x <- split(seq(nrow(coords)), clusters)
	result <- sapply(x, function(a) apply(coords[a, ],2,median))
	return(result)
}



lab_pos <- agg_coord_by_cluster(int.obj@reductions$tsne@cell.embeddings, int.obj@meta.data$Fine_clusters)
lab_pos_umap <- agg_coord_by_cluster(int.obj@reductions$umap@cell.embeddings, int.obj@meta.data$Fine_clusters)
lab_id <- colnames(lab_pos)
clust_autoanno <- sapply(strsplit(Cluster_guesses, "\\."), function(a){a[[2]]})
clust_autoanno[clust_autoanno == "UnidentifiedHep"] <- "Hepatocyte";
clust_autoanno[clust_autoanno == "CD3abTcells"] <- "Tcells";
clust_autoanno[clust_autoanno == "PericentralHep"] <- "Hepatocyte";
clust_autoanno[clust_autoanno == "Unknown"] <- "Hepatocyte";
clust_autoanno[clust_autoanno == "AntibodysecretingBcells"] <- "Bcells";
clust_autoanno[clust_autoanno == "CentralvenousLSECs"] <- "cvLSEC";
clust_autoanno[clust_autoanno == "InflamatoryMacrophages"] <- "inf-Macrophage";
clust_autoanno[clust_autoanno == "Non-inflammatoryMacrophages"] <- "non-Macrophage";


Colour_scheme <- rbind(
	c("cvLSECs", "goldenrod1"),
	c("LSECs", "violet"),
	c("Bcells", "darkorchid"),
	c("Tcells", "forestgreen"),
	c("Macrophage", "cornflowerblue"),
	c("InflamatoryMacrophages", "navy"),
	c("Non-inflamatoryMacrophages", "dodgerblue"),
	c("Erythoidcells", "darkred"),
	c("Hepatocytes", "firebrick2"))
Colour_scheme <- Colour_scheme[order(Colour_scheme[,1]),]

int.obj@meta.data$autoanno_c <- clust_autoanno[int.obj@meta.data$Fine_clusters]

require("ggplot2")
png("General_autoanno_label_harmony_tsne.png", width=7.5, height=6, units="in", res=300)
DimPlot(int.obj, reduction = "tsne", group.by = "autoanno_c", pt.size = .1)+scale_color_manual(values=Colour_scheme[,2])+annotate("text", x=lab_pos[1,], y=lab_pos[2,], label=lab_id, colour="grey35")
dev.off();
png("General_autoanno_label_harmony_umap.png", width=7.5, height=6, units="in", res=300)
DimPlot(int.obj, reduction = "umap", group.by = "autoanno_c", pt.size = .1)+scale_color_manual(values=Colour_scheme[,2])+annotate("text", x=lab_pos_umap[1,], y=lab_pos_umap[2,], label=lab_id, colour="grey35")
dev.off();


metadata <- read.table("Metadata20LiverMapPlusParams.csv", sep=",", header=T)

pats <- levels(factor(int.obj@meta.data$donor))
sex <- as.character(metadata[match(pats, metadata$Name),"sex"])
age <- as.character(metadata[match(pats, metadata$Name),"age"])
bmi <- as.character(metadata[match(pats, metadata$Name),"BMI"])
sex[!sex %in% c("M", "F")] <- "Unk"
int.obj@meta.data$sex <- sex[factor(int.obj@meta.data$donor)]
int.obj@meta.data$age <- age[factor(int.obj@meta.data$donor)]
int.obj@meta.data$bmi <- bmi[factor(int.obj@meta.data$donor)]


png("General_donor_label_harmony_tsne.png", width=6.5, height=6, units="in", res=300)
DimPlot(int.obj, reduction = "tsne", group.by = "donor", pt.size = .1)
dev.off();

png("General_sex_harmony_tsne.png", width=6.5, height=6, units="in", res=300)
DimPlot(int.obj, reduction = "tsne", group.by = "sex", pt.size = .1)+scale_color_manual(values=c("salmon","cornflowerblue", "grey50"))
dev.off();

png("General_sex_harmony_umap.png", width=6.5, height=6, units="in", res=300)
DimPlot(int.obj, reduction = "umap", group.by = "sex", pt.size = .1)+scale_color_manual(values=c("salmon","cornflowerblue", "grey50"))
dev.off();

png("General_cluster_harmony_umap.png", width=6.5, height=6, units="in", res=300)
DimPlot(int.obj, reduction = "umap", group.by = "Fine_clusters", pt.size = .1)
dev.off();

png("General_cluster_harmony_tsne.png", width=6.5, height=6, units="in", res=300)
DimPlot(int.obj, reduction = "tsne", group.by = "Fine_clusters", pt.size = .1)
dev.off();

## Super General summary & metadata analysis ##

int.obj$general_lab <- int.obj$autoanno_c
int.obj$general_lab[grepl("Macropha", int.obj$general_lab)] <- "Macrophage"
int.obj$general_lab[grepl("LSEC", int.obj$general_lab)] <- "LSECs"

table(int.obj$general_lab)
table(int.obj$general_lab)/sum(table(int.obj$general_lab))

tab <- table(int.obj$general_lab, int.obj$sex)
tab <- t(t(tab)/colSums(tab))
tab <- tab[,-3]


png("General_types_by_gender.png", width=9, height=5, units="in", res=300)
par(mfrow=c(1,2))
pie(tab[,1], col=Colour_scheme[c(1,3,4,6,7,9),2], main="Male")
pie(tab[,2], col=Colour_scheme[c(1,3,4,6,7,9),2], main="Female")
dev.off()


# categorize continuous variables.
int.obj@meta.data$age <- as.numeric(int.obj@meta.data$age)
int.obj@meta.data$bmi <- as.numeric(int.obj@meta.data$bmi)
bins_Age <- c(0,35, 60,300)
names_age <- c("young", "adult", "elderly")
bins_BMI <- c(0, 15,18,25,30,35,40, 100);
names_bmi <- c("very weight", "underweight", "normal", "overweight", "obese", "very obese")

binned <- cut(int.obj@meta.data$age, breaks=bins_Age)
lab <- names_age[binned]
int.obj@meta.data$age_cat <- lab

binned <- cut(int.obj@meta.data$bmi, breaks=bins_BMI)
lab <- names_bmi[binned]
int.obj@meta.data$bmi_cat <- lab


tab <- table(int.obj$general_lab, int.obj$age_cat)
tab <- t(t(tab)/colSums(tab))
tab <- tab[,c(3,1,2)]


png("General_types_by_age.png", width=9*3/2, height=5, units="in", res=300)
par(mfrow=c(1,3))
pie(tab[,1], col=Colour_scheme[c(1,3,4,6,7,9),2], main="Young")
pie(tab[,2], col=Colour_scheme[c(1,3,4,6,7,9),2], main="Adult")
pie(tab[,3], col=Colour_scheme[c(1,3,4,6,7,9),2], main="Elderly")
dev.off()



int.obj$type_anno <- int.obj$scmap_anno2
int.obj$type_anno[grepl("Macropha", int.obj$type_anno)] <- "Macrophage"
int.obj$type_anno[grepl("LSEC", int.obj$type_anno)] <- "LSECs"
int.obj$type_anno[grepl("Tcell", int.obj$type_anno)] <- "Tcells"
int.obj$type_anno[grepl("Bcell", int.obj$type_anno)] <- "Bcells"
int.obj$type_anno[grepl("Hep", int.obj$type_anno)] <- "Hepatocyte"









q();


# Output:

# Auto-anno table:
# % cells/ref cluster - scmap cluster, scmap cell
# marker enrichment for markers from reference clusters.

# Cluster x Annotation table

## Cluster specific Marker genes ###
# Get markers for fine & coarse clusters
# problem: seurat integrated only has 1000 genes....
# harmony one only has ~3,000 - background for GSEA!!!

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

## Significant with MAST
if (!file.exists("harmony_MASTmodel.rds")) {
	set.seed(8817)

	require("MAST")

	sca <- FromMatrix(as.matrix(int.obj@assays$RNA@data)) ##### <---- not sparse this is very Mem heavy!!

	ngenes <- colSums(int.obj@assays$RNA@counts > 0)
	ncounts <- colSums(int.obj@assays$RNA@counts)
	donor <- int.obj@meta.data$donor
	cluster <- int.obj@meta.data$Fine_clusters

	# This is slow with so many cells.... Perhaps move to H4h cluster.
	mod <- zlm(~cluster + donor + ngenes + ncounts, sca)

	saveRDS(mod, "harmony_MASTmodel.rds")
}


#res1 <- summary(mod, doLRT="cluster1")
#res <- res1$datatable[res1$datatable$contrast=="cluster1" & res1$datatable$component %in% c("C","D"),]
#res$fdr <- p.adjust(unlist(res[,4]), method="fdr")

#con <- res[res$component =="C",]
#dis <- res[res$component =="D",]

#sig <- (con$fdr < 0.01 | dis$fdr < 0.01) & sign(con$coef) == sign(dis$coef)
#out <- cbind(con[sig,], dis[sig,])
#out$overall <- out[,7]+out[,16]
#out <- out[order(out$overall, decreasing=T),]
#write.table(out, file=paste(this_contrast,"_MAST_output.txt", sep=""), col.names=T, row.names=F)
 


## Pseudobulks ##
make_pseudobulk_table <- function(mat, clusters, indis) {
	c <- split(seq(ncol(mat)), clusters);
	donor_expr <- function(MAT, ds) {
		d <- split(seq(ncol(MAT)), ds);
		mus <- sapply(d, function(group) my_rowSums(MAT[,group]))
		return(mus);
	}
	pseudobulk_expr <- lapply(c, function(clust) {
		d_expr <- donor_expr(mat[,clust], indis[clust]);
		return(d_expr);
	})
	return(pseudobulk_expr);
}

make_pseudobulk_v2 <- function(mat, clusters, indis) {
	tab <- apply(mat, 1, function(x){tapply(x, list(clusters, indis), sum)})
	c_lab <- rep(levels(as.factor(clusters)), each=length(levels(as.factor(indis))))
	i_lab <- rep(levels(as.factor(indis)), times=length(levels(as.factor(clusters))))
	rownames(tab) <- paste(c_lab, i_lab, sep="_")
	return(list(bulks=t(tab), groups=c_lab, indis=i_lab));
}

pseudobulk <- make_pseudobulk_v2(int.obj@assays$RNA@counts,  int.obj@meta.data$Fine_clusters, int.obj@meta.data$donor)

# subset
keep <- which(!is.na(colSums(pseudobulk$bulks)));
pseudobulk$bulks <- pseudobulk$bulks[,keep]
pseudobulk$groups <- pseudobulk$groups[keep]
pseudobulk$indis <- pseudobulk$indis[keep]


# run edger - No longer super biased but I'm still not happy with the results: l2fc are really high significance is really low.....
cluster = "1"
coef = paste("groups",cluster,sep="")

require("edgeR")
edger.obj <- DGEList(counts=pseudobulk$bulks, group=pseudobulk$groups);
edger.obj <- calcNormFactors(edger.obj)
design <- model.matrix(~pseudobulk$groups + pseudobulk$indis)
edger.obj <- estimateDisp(edger.obj, design)
fit <- glmQLFit(edger.obj, design) #glmQLFit/glmFit
contrast_vec <- rep(-1, ncol(design));
this_group = which(grepl(coef, colnames(design)));
if (length(this_group) > 0) {
	contrast_vec[this_group[1]] <- 1;
}
contrast_vec[1] <- 0
lrt <- glmQLFTest(fit, contrast=contrast_vec) #glmQLFTest/glmLRT
topTags(lrt, 100)


## Ranked list using my method ##

# Matrix of relative expression levels across clusters
# per cluster get means per donor and weight by overall donor freqs

get_rel_expression <- function(mat, clusters, donors) {
	c <- split(seq(ncol(mat)), clusters);
	donor_freqs <- table(donors)/length(donors)
	# avg expression per donor in this cluster
	donor_expr <- function(MAT, ds) {
		d <- split(seq(ncol(MAT)), ds);
		mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
		return(mus);
	}
	clust_expr <- sapply(c, function(clust) {
		d_expr <- donor_expr(mat[,clust], donors[clust]);
		# weight by overall frequency of donors
		freqs <- donor_freqs[match(colnames(d_expr), names(donor_freqs))]
		freqs <- as.vector(freqs)/sum(freqs)
		c_expr <- my_rowSums(t(t(d_expr)*freqs))
		return(c_expr);
	})
	return(clust_expr)
}

source("Setup_autoannotation.R")

## Ranked list
cluster_means <- get_rel_expression(int.obj@assays$RNA@data, int.obj@meta.data$Fine_clusters, int.obj@meta.data$donor)
cluster_detect <- get_rel_expression(int.obj@assays$RNA@counts > 0, int.obj@meta.data$Fine_clusters, int.obj@meta.data$donor)

mark_mean <- my_markers(cluster_means)
mark_detect <- my_markers(cluster_detect);

