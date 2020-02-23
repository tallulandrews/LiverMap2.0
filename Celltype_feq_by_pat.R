require("Seurat")
set.seed(3918)
int.obj <- int.obj <- readRDS("integration_harmony_plus_analysis.rds")

metadata <- read.table("Metadata20LiverMapPlusParams.csv", sep=",", header=T)

# Add metadata to int.obj

meta.dat.cols <- c("age", "sex", "BMI")
for (i in meta.dat.cols) {
	int.obj@meta.data[,i] <- metadata[match(int.obj@meta.data$donor, metadata$Sample),i]
}

# categorize continuous variables.
bins_Age <- c(0,35, 60,300)
names_age <- c("young", "adult", "elderly")
bins_BMI <- c(0, 15,18,25,30,35,40, 100);
names_bmi <- c("very weight", "underweight", "normal", "overweight", "obese", "very obese")

binned <- cut(int.obj@meta.data$age, breaks=bins_Age)
lab <- names_age[binned]
int.obj@meta.data$age_cat <- lab

binned <- cut(int.obj@meta.data$BMI, breaks=bins_BMI)
lab <- names_bmi[binned]
int.obj@meta.data$bmi_cat <- lab


# Cell-type freqs by donor
counts <- table(int.obj@meta.data$autoanno_c, int.obj@meta.data$donor)
freqs <- t(t(counts)/colSums(counts))
null_freqs <- rowSums(counts)/ sum(counts)

# correct for blood proportion:
types_in_blood_and_liver <- c("Tcells", "Macrophage", "inf-Macrophage", "Bcells","non-Macrophage")
type_in_blood_only <- c("Erythroidcells");



# X = % cells from blood
# ?_blood = proportion of all cells that are cell type ? from blood
# ?_all = proportion of all cells that are cell type ? from either blood or liver
# Mac_blood/Mac_all = X, T_blood/T_all = X, .....
# Mac_blood+T_blood+Erythro_blood.... = X
# X = Erythro_blood + X*Mac_all + X*T_all +....
# (1- Mac_all -T_all ...)*X = Erythro_blood

remove_blood_from_cell_type_freqs <- function(freqs) {
	prop_blood <- freqs[names(freqs) %in% type_in_blood_only]/
			(1-freqs[names(freqs) %in% types_in_blood_and_liver]);
	# Mac_liver = (1-X)*Mac_all
	freqs[names(freqs) %in% types_in_blood_and_liver] <- (1-prop_blood)*freqs[names(freqs) %in% types_in_blood_and_liver]

	freqs <- freqs[!names(freqs) %in% type_in_blood_only]
	freqs <- freqs/sum(freqs);
	return(freqs)
}
corrected_freqs <- apply(freqs, 2, remove_blood_from_cell_type_freqs)
corrected_counts <- round(t(t(corrected_freqs)*colSums(counts[-3,])))

test <- chisq.test(corrected_counts)
diff <- (test$observed-test$expected)/test$expected

meta <- metadata[match(colnames(diff), metadata[,1]),]
