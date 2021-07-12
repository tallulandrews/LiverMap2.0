
require(Seurat)
require(MAST)
require(ggplot2)
require(lme4)

full_metadata_with_subcluster_anno <- Sys.glob("*fullmetadata.rds")
regular_objs <- Sys.glob("*varimax_Subcluster.rds")
full_objs <- Sys.glob("")


args <- commandArgs(trailingOnly=TRUE);

this_i <- as.numeric(args[1]);

cell_type <- regular_objs[this_i]
cell_type <- sapply(strsplit(cell_type, "/"), function(x) { x[[length(x)]] })
cell_type <- sapply(strsplit(cell_type, "_"), function(x) { x[[1]] })


obj <- readRDS(regular_objs[this_i])
obj@meta.data <- readRDS(full_metadata_with_subcluster_anno[this_i])

tmp <- table(obj@meta.data$sample);
exclude <- names(tmp)[tmp < 10]

if (length(exclude) > 0) {
	obj <- obj[, ! (obj@meta.data$sample %in% exclude)]
}

# Subsample large samples - Still way to slow and singular!!!
set.seed(2819)
max_cells <- 1000
if (max(tmp) > max_cells) {
	all_cells_exclude <- c();
	to.subset <- names(tmp)[which(tmp > max_cells)]
	for(sample_i in to.subset) {
		these_cells <- which(obj@meta.data$sample == sample_i)
		exclude <- sample(these_cells, length(these_cells)-max_cells)
		all_cells_exclude <- c(all_cells_exclude, exclude);
	}
	obj <- obj[,-all_cells_exclude]
}
## MAST Linear Mixed Model ##

obj_sce <- as.SingleCellExperiment(obj)
tmp <- as.matrix(obj@assays$RNA@data); rownames(tmp) <- rownames(obj); colnames(tmp) <- colnames(obj)
obj_sca <- FromMatrix(tmp, colData(obj_sce), rowData(obj_sce), class="SingleCellAssay")

# Use PCs instead of individual genes.
pca_mat <- t(obj@reductions$pca@cell.embeddings)
varimax_mat <- t(obj@reductions$varimax@cell.embeddings)
obj_sca <- FromMatrix(varimax_mat, colData(obj_sce), DataFrame(t(obj@reductions$varimax@feature.loadings)), class="SingleCellAssay", check_sanity=FALSE)



#zlm.out2 <- zlm(~donor_sex + donor_age + assay_type + (1|sample) , obj_sca, method="glmer", ebayes=FALSE)
zlm.out2 <- zlm(~donor_sex + donor_age + assay_type + (1|sample) , obj_sca, method="glmer", ebayes=FALSE)

coefAndCI <- as.data.frame(summary(zlm.out2, logFC=FALSE)$datatable)

# SEX #
res_sex <- coefAndCI[coefAndCI[,3] == "donor_sexM" & !is.na(coefAndCI[,6]),]

max_coef.sex <- aggregate(res_sex$coef, by=list(res_sex$primerid), function(x) {x[abs(x) == max(abs(x))]})
zlm.lr.sex <- lrTest(zlm.out2, "donor_sex")
pval.sex <- zlm.lr.sex[,,3][,"hurdle"]

sex.output <- cbind(max_coef.sex, pval.sex[match( max_coef.sex[,1], names(pval.sex))])
sex.output$qval <- p.adjust(sex.output[,3], method="fdr")

print(sum(sex.output$qval < 0.05))
saveRDS(sex.output, paste(cell_type, "_MM_DE_varimax_sex.rds", sep=""))


# AGE #
res_age <- coefAndCI[coefAndCI[,3] == "donor_age" & !is.na(coefAndCI[,6]),]

max_coef.age <- aggregate(res_age$coef, by=list(res_age$primerid), function(x) {x[abs(x) == max(abs(x))]})
zlm.lr.age <- lrTest(zlm.out2, "donor_age")
pval.age <- zlm.lr.age[,,3][,"hurdle"]

age.output <- cbind(max_coef.age, pval.age[match( max_coef.age[,1], names(pval.age))])
age.output$qval <- p.adjust(age.output[,3], method="fdr")

print(sum(age.output$qval < 0.05))
saveRDS(age.output, paste(cell_type, "_MM_DE_varimax_age.rds", sep=""))


# BMI #
#res_bmi <- coefAndCI[coefAndCI[,3] == "donor_bmi" & !is.na(coefAndCI[,6]),]

#max_coef.bmi <- aggregate(res_bmi$coef, by=list(res_bmi$primerid), function(x) {x[abs(x) == max(abs(x))]})
#zlm.lr.bmi <- lrTest(zlm.out2, "donor_bmi")
#pval.bmi <- zlm.lr.bmi[,,3][,"hurdle"]

#bmi.output <- cbind(max_coef.bmi, pval.bmi[match( max_coef.bmi[,1], names(pval.bmi))])
#bmi.output$qval <- p.adjust(bmi.output[,3], method="fdr")

#print(sum(bmi.output$qval < 0.05))
#saveRDS(bmi.output, paste(cell_type, "_MM_DE_bmi.rds", sep=""))

