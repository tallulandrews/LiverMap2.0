# too look into: Monocle DE, diffcyt : https://bioconductor.org/packages/release/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html

tidy_Seurat_output <- function(out) {
	out$score <- out$avg_logFC*(out$pct.1-out$pct.2)
	out <- out[order(out$score, decreasing=T),]
	return(out)
}

args <- commandArgs(trailingOnly=TRUE)
# file
# clustername
# projname
# cluster_col_name (Use_clusters)
# confounder_col_name (donor)

require(Seurat)
obj <- readRDS(args[1])
source("/cluster/home/tandrews/scripts/LiverMap2.0/Map2.2_add_metadata.R")
meta_data <- get_metadata(obj@meta.data$donor);

obj@meta.data <- cbind(obj@meta.data, meta_data);

if (length(args) >= 4) {
	obj@meta.data$Use_clusters <- obj@meta.data[,args[4]]
}
if (length(args) >= 5) {
	obj@meta.data$donor <- obj@meta.data[,args[5]]
}

exclude1 <- is.na(obj@meta.data$Use_clusters) | obj@meta.data$Use_clusters == "";
exclude2 <- is.na(obj@meta.data$donor) | obj@meta.data$donor == "";
obj <- obj[,!(exclude1 | exclude2)];

out_seur_wilcox <- FindMarkers(obj, ident.1=args[2], group.by="Use_clusters", logfc.threshold=-Inf, min.pct=0.1, min.diff.pct = -Inf, test.use="wilcox")

require(MAST)
out_seur_MAST <- FindMarkers(obj, ident.1=args[2], group.by="Use_clusters", logfc.threshold=-Inf, min.pct=0.1, min.diff.pct = -Inf, test.use="MAST", latent.vars=obj$donor)

all_DE <- list(seur_wilcox=out_seur_wilcox, seur_mast=out_seur_MAST)
saveRDS(all_DE, paste(args[3], args[2], "DE.rds", sep="_"))

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

# too slow!
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


## Pseudobulk ##
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

bulks <- get_pseudobulk(obj@assays$RNA@counts, 
			factor(obj@meta.data$Use_clusters), 
			factor(obj@meta.data$donor))
labs <- strsplit(colnames(bulks), "_")
c_labs <- sapply(labs, function(x){unlist(x[[1]])})
d_labs <- sapply(labs, function(x){unlist(x[[2]])})

bulks_file <- paste(args[3], "_pseudobulk_mat.rds", sep="")
saveRDS(list(mat=bulks, cluster=c_labs, donor=d_labs), bulks_file);


if(!identical(c_labs, d_labs)) {
require("edgeR")
edger_obj <- DGEList(bulks, samples=data.frame(cluster=c_labs, donor=d_labs), group=c_labs)
edger_obj <- calcNormFactors(edger_obj)
design <- model.matrix(~cluster+donor, data=edger_obj$samples)
design <- design[,colSums(design) > 0];

coef_names <- colnames(design);

edger_file <- paste(args[3], "_pseudobulk_edger.rds", sep="")
if (!file.exists(edger_file)) {
        edger_obj <- estimateDisp(edger_obj, design)
        fit <- glmQLFit(edger_obj, design)
        saveRDS(fit, edger_file);
} else {
        fit <- readRDS(edger_file)
}

contrast_vec <-rep(0, length(coef_names));
this_cluster <- paste("cluster", args[2], sep="");
contrast_vec[grepl("cluster", coef_names)] <- -1;
contrast_vec[1] <- -1
if (this_cluster %in% coef_names) {
        contrast_vec[coef_names==this_cluster]<- 1;
} else {
        contrast_vec[1] <- 1;
}

this_vs_all <- glmQLFTest(fit, contrast=contrast_vec)
out <- topTags(this_vs_all, nrow(obj));
all_DE <- list(seur_wilcox=out_seur_wilcox, seur_mast=out_seur_MAST, edger=out)
saveRDS(all_DE, paste(args[3], args[2], "_pseudobulk_DE.rds", sep="_"))
}

#######   END ########


### Diffcyt? ###
#require(diffcyt)

## make input:
#input_list <- list()
#experiment_info <- c()
#for (d in unique(obj@meta.data$donor)) {
#	for (c in unique(obj@meta.data$Use_clusters)) {
#		this <- obj@assays$RNA@counts[,obj@meta.data$donor ==d & obj@meta.data$Use_cluster ==c]
#		sample <- paste(d, c, sep="_")
#		if (sum(obj@meta.data$donor ==d & obj@meta.data$Use_cluster ==c) > 0) {
#			input_list[[sample]] <- as.matrix(t(this));
#			experiment_info <- rbind(experiment_info, c(sample, d, c))
#		}
#	}
#}
#colnames(experiment_info) <- c("sample_id", "patient_id", "cell_type")
#
##input data - list of matrices for each 'sample'
##experiment_info - dataframe of metadata for samples
#
##marker_info - dataframe of metadata for markers
#
##preprocess
#d_se <- prepareData(input_list, data.frame(experiment_info), marker_info=data.frame(gene=rownames(obj), marker_name=rownames(obj)))
##d_se <- transformData(d_se) <- don't do this b/c designed for CyTOF
##d_se <- generateClusters(d_se, seed_clustering = 123) <- FlowSOM clustering
## Calculate cluster cell counts
#d_counts <- calcCounts(d_se)
#
#d_se@elementMetadata$cluster_id <- matched_cell_cell_type
#
## Calculate cluster medians
#d_medians <- calcMedians(d_se)
#
## Create design matrix
## note: selecting columns containing group IDs and patient IDs (for an 
## unpaired dataset, only group IDs would be included)
#design <- createDesignMatrix(
#  data.frame(experiment_info), cols_design = c("cell_type", "patient_id")
#)
## Create contrast matrix
#coef_names <- colnames(design)
#contrast_vec <-rep(0, length(coef_names));
#this_cluster <- paste("cell_type", args[2], sep="");
#contrast_vec[grepl("cell_type", coef_names)] <- -1;
#contrast_vec[1] <- -1
#if (this_cluster %in% coef_names) {
#        contrast_vec[coef_names==this_cluster]<- 1;
#} else {
#        contrast_vec[1] <- 1;
#}
#contrast <- createContrast(contrast_vec)
#
#
## check
#nrow(contrast) == ncol(design)
#
## Test for differential abundance (DA) of clusters
#res_DA <- testDA_edgeR(d_counts, design, contrast)
#
## display table of results for top DA clusters
#topTable(res_DA, format_vals = TRUE)
#
#
