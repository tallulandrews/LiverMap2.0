


clean_sample_names <- function(sample_names) {
	sample_names <- gsub("-NPC$", "", sample_names)
	sample_names <- gsub("_NPC$", "", sample_names)
	return(sample_names);
}


bin_age <- function(age_vec) {
	bins_Age <- c(0,35, 60,300)
	names_age <- c("young", "adult", "elderly")
	binned <- cut(age_vec, breaks=bins_Age)
	lab <- names_age[binned]
	return(lab)
}

bin_BMI <- function(weight_vec) {
	bins_BMI <- c(0, 15,18,25,30,35,40, 100);
	names_bmi <- c("very weight", "underweight", "normal", "overweight", "obese", "very obese")
	binned <- cut(weight_vec, breaks=bins_BMI)
	lab <- names_bmi[binned]
	return(lab)
}

get_metadata <- function(donor_vec) {
	file <- "/cluster/home/tandrews/scripts/LiverMap2.0/LiverMap_SampleProcessingParams.csv";
	donor_vec <- gsub("_NPC", "-NPC", donor_vec)
	metadata <- read.delim(file, header=T, sep=",", stringsAsFactors=FALSE)
	if(sum(grepl("-", donor_vec)) > 0) {
		matching <- match(donor_vec, metadata$Name)
	} else {
		matching <- match(donor_vec, metadata$Sample)
	}
	metadata <- metadata[matching,]
	metadata <- metadata[,c("age", "sex", "BMI")]
	
	metadata$AGE <- factor(bin_age(metadata$age))
	metadata$WEIGHT <- factor(bin_BMI(metadata$BMI))
	metadata$SEX <- factor(metadata$sex)
	return(metadata)
}
