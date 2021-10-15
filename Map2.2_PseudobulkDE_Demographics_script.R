#require(edgeR)
source("~/scripts/LiverMap2.0/My_R_Scripts.R")
require(matrixStats)

bulks <- readRDS("Merged_EmptyOnly_obj_Map2.2_allgenes_3pr_pseudobulks.rds")

pseudobulk_mat <- bulks$psuedobulk
metadata <- bulks$metadata

args <- commandArgs(trailingOnly=TRUE)

# Issue = Background / Batch effects!
# Solution = scale by each donor.
# Then use glm (not NB anymore).

# Scale & normalize pseudobulks!
pseudobulk_mat <- pseudobulk_mat[!grepl("^MT-", rownames(pseudobulk_mat)),]
pseudobulk_norm <- t(t(pseudobulk_mat)/colSums(pseudobulk_mat) * min(colSums(pseudobulk_mat)))
pseudobulk_scaled <- pseudobulk_norm

avgs <- group_rowmeans(pseudobulk_norm, metadata$donor)
vars <- group_rowvars(pseudobulk_norm, metadata$donor)
vars <- sqrt(vars)
vars[vars == 0] <- 0.0001
for(d in levels(metadata$donor)) {
	d_col <- metadata$donor == d
	d_z_corrected <- (pseudobulk_norm[,d_col]-avgs[,d])/vars[,d]
	pseudobulk_scaled[,d_col] <- d_z_corrected;
}


for(type in unique(metadata$type)) {
	mat <- pseudobulk_scaled[,metadata$type == type]
	mat <- mat[rowMeans(abs(mat)) > 0.1, ];
	mat <- mat[!grepl("^MT-", rownames(mat)),]
	meta <- metadata[metadata$type == type,]

	# Straight Means
	F_mean <- rowMeans(mat[,meta$SEX=="F"])
	F_var <- rowVars(mat[,meta$SEX=="F"])/sqrt(sum(meta$SEX=="F")); 
	names(F_var) <- names(F_mean)
	M_mean <- rowMeans(mat[,meta$SEX=="M"])
	M_var <- rowVars(mat[,meta$SEX=="M"])/sqrt(sum(meta$SEX=="M")); 
	names(M_var) <- names(M_mean)
	young_mean <- rowMeans(mat[,meta$AGE=="young"])
	young_var <- rowVars(mat[,meta$AGE=="young"])/sqrt(sum(meta$AGE=="young")); 
	names(young_var) <- names(young_mean);
	adult_mean <- rowMeans(mat[,meta$AGE=="adult"])
	adult_var <- rowVars(mat[,meta$AGE=="adult"])/sqrt(sum(meta$AGE=="adult")); 
	names(adult_var) <- names(adult_mean);
	elderly_mean <- rowMeans(mat[,meta$AGE=="elderly"])
	elderly_var <- rowVars(mat[,meta$AGE=="elderly"])/sqrt(sum(meta$AGE=="elderly")); 
	names(elderly_var) <- names(elderly_mean);
	
	sex <- meta$SEX
	age <- meta$AGE
	age_c <- meta$age

	# gene-specific DE testing:
	#x <- mat["IL10RA",]
	sex_de <- function(x) {
		test = summary(glm(x~sex+age))
		sex_p_val <- test$coefficients[2,4]
		sex_effect <- test$coefficients[2,1]
		return(c(sex_effect, sex_p_val))
	}
	sex_res <- t(apply(mat, 1, sex_de))
	sex_res <- cbind(sex_res, p.adjust(sex_res[,2], method="fdr"))

	coefficient_contrasts <- function(x) {
		test = glm(x~sex + age_c);
		coefs <- test$coefficients
		#M_vs_F <- coefs[2]
		#F_age <- coefs[3]
		#M_age <- coefs[4]
		
		return(coefs)
	}
	model_res <- t(apply(mat, 1, coefficient_contrasts))
	model_res <- cbind(model_res, F_mean, M_mean, young_mean, adult_mean, elderly_mean)

	# EdgeR - this doesn't work b/c of all the background
	#edger_obj <- DGEList(mat, samples=data.frame(sex=meta$SEX, age=meta$AGE), group=meta$SEX)
	#edger_obj <- calcNormFactors(edger_obj)
	#design <- model.matrix(~sex+age, data=edger_obj$samples)
	#design <- design[,colSums(design) > 0];

	#coef_names <- colnames(design);
        #edger_obj <- estimateDisp(edger_obj, design)
        #fit <- glmQLFit(edger_obj, design)
	#contrast_mat <- rbind(c(1, -1, 0, 0), c(0, 0, 1, -1), c(-1, 0, 1, 0), c(1, 0, 0, -1))
	#contrast_names <- c("F_vs_M", "elderly_vs_young", "elderly_vs_adult", "adult_vs_young")
	#res <- list();
	#for (i in 1:nrow(contrast_mat)) {
	#	this_vs_all <- glmQLFTest(fit, contrast=contrast_mat[i,])
	#	out <- topTags(this_vs_all, nrow(mat));
	#	res[[contrast_names[i]]] <- out;	
	#}
	saveRDS(model_res, paste("Merged_Empty_3pr_CellType", type, "pseudobulk_DE.rds", sep="_"))
}		

