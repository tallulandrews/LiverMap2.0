
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/SubColour_Scheme.R")
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/My_R_Scripts.R")
# Individual Variability in frequency
donor2meta <- function(meta.data, meta_col) {
	tmp <- table(meta.data$donor, meta.data[,meta_col])
	out <- apply(tmp, 1, function(x) {colnames(tmp)[which(x == max(x))]})
	return(out)
}

indivar_freq <- function(cell_type, donor, metadata, donor_colour=rainbow(length(unique(donor))), min.cells=10) {
	# Clean up
	if (is.null(names(donor_colour))) {
		names(donor_colour) <- levels(factor(donor))
	}
		
	freq <- table(donor, cell_type)
	donor_age <- as.numeric(donor2meta(metadata, "donor_age"))
	donor_sex <- donor2meta(metadata, "donor_sex")
	
	
	ns <- rowSums(freq)

	exclude <- ns < min.cells
	freq <- freq[!exclude,]
	donor_age <- donor_age[!exclude]
	donor_sex <- donor_sex[!exclude]

	# colour by sex
	donor_colour <- colours_sex[factor(donor_sex, levels=names(colours_sex))]
	names(donor_colour) <- rownames(freq)

	#donor_colour <- donor_colour[match(rownames(freq), names(donor_colour))]
	#donor_age <- factor(donor_age, levels=c("young", "adult", "elderly"))
	#donor_sex <- factor(donor_sex, levels=c("F", "M"))

	ns <- rowSums(freq)
	freq <- freq/ns
	freq[freq==0] <- 10^-5
	
	freq <- freq[order(donor_age, decreasing=F),]

	xlim = c(0, prod(dim(freq))+ncol(freq))
	xes = c();
	yes = c();
	CI95 = c();
	for(j in (1:ncol(freq))) {
		xes <- c(xes, (1:nrow(freq)) +(nrow(freq)+1)*(j-1))
		yes <- c(yes, freq[,j])
		CI95 <- c(CI95, sqrt(freq[,j]*(1-freq[,j])/ns)*2)
	}
	plot(xes, yes, xaxt="n", xlab="cell type", ylab="frequency", pch=16, col=donor_colour, xlim=xlim)
	v_lines <-(1:xlim[2])[! 1:xlim[2] %in% xes] 
	abline(v= v_lines[1:(length(v_lines)-1)], lty=2)
	arrows(xes, yes-CI95, xes, yes+CI95, length=0)

	label_pos <- seq(from=1, to= max(xlim), length=ncol(freq)+1)
	label_pos <- (label_pos[1:(ncol(freq))]+label_pos[2:(ncol(freq)+1)])/2
	axis(1, at=label_pos, labels=colnames(freq), tick=FALSE)
	return(freq)
}

indi_expr <- function(obj, cluster_col="Subcluster_Manual", donor_col="donor", min.cells=10, tag="test") {
	require(Seurat)
	require(cluster)
	require(gplots)
	set.seed(4302)
	pdf(paste(tag, "indivar_silhouette.pdf", sep="_"), width=8, height=8)

	ns <- table(obj@meta.data[,donor_col])

	mat <- get_pseudobulk_means(obj@assays$RNA@data, obj@meta.data[,cluster_col], obj@meta.data[,donor_col])
	
	mat <- mat[rownames(mat) %in% VariableFeatures(obj),]
	mat <- mat[,ns > min.cells]
	heatmap.2(cor(mat, method="pearson"), Rowv=FALSE, Colv=FALSE, trace="none")
	
	d <- as.dist(1-cor(mat, method="pearson"))
	id <- factor(unlist(lapply(strsplit(colnames(mat), "_"), function(x){x[[1]]})))
	sil <- silhouette(as.numeric(id), dist=d)
	plot(sil)

	scores <- c();
	stderr <- c();
	for (i in unique(id)) {
		scores <- c(scores, mean(as.matrix(1-d)[id==i, id==i]))
		stderr <- c(stderr, sd(as.matrix(1-d)[id==i, id==i])/sqrt(sum(id==i))) # Fix 4 May 2021
	}
	bar_loc <- barplot(scores, name=unique(id), ylab="Cross Donor Correlation (average)", ylim=c(0,1), las=1)
	arrows(bar_loc, scores, bar_loc, scores+2*stderr, angle=90)
	legend("topright", lty=1, c("95% CI"), bty="n")
	dev.off()
	return(list(silhouette=sil, intracor=scores, pseudobulks=mat))
}



indivar_DE_vis <- function(gene_vec, donor, metadata, donor_colour=rainbow(length(unique(donor))), min.cells=10, axis.interval=3, type="sex") {
	# Clean up
	if (is.null(names(donor_colour))) {
		names(donor_colour) <- levels(factor(donor))
	}
	
	# Average across sets of genes
	gene_vec <- my_colMeans(gene_vec)
	
	d <- split(seq(length(gene_vec)), metadata$donor)
	d_means <- sapply(d, function(group){mean(gene_vec[group])})
	d_medians <- sapply(d, function(group){median(gene_vec[group])})
	d_stderr <- sapply(d, function(group){sd(gene_vec[group])/sqrt(length(group))})
	
	donor_age <- as.numeric(donor2meta(metadata, "donor_age"))
	donor_sex <- donor2meta(metadata, "donor_sex")
	donor_rej <- donor2meta(metadata, "trans.rejected")
	
	
	ns <- table(metadata$donor)

	exclude <- ns < min.cells
	d_means <- d_means[!exclude]
	d_medians <- d_medians[!exclude]
	d_stderr <- d_stderr[!exclude]

	donor_age <- donor_age[!exclude]
	donor_sex <- donor_sex[!exclude]
	donor_rej <- donor_rej[!exclude]


	# colour by sex
	donor_colour <- colours_sex[factor(donor_sex, levels=names(colours_sex))]
	names(donor_colour) <- names(d_means)
	if (type == "reject") {
		donor_colour <- colours_reject[factor(donor_rej, levels=names(colours_reject))]
	}


	#donor_colour <- donor_colour[match(rownames(freq), names(donor_colour))]
	#donor_age <- factor(donor_age, levels=c("young", "adult", "elderly"))
	#donor_sex <- factor(donor_sex, levels=c("F", "M"))

	age_order <- order(donor_age, decreasing=F)
	d_means <- d_means[age_order]
	d_medians <- d_medians[age_order]
	d_stderr <- d_stderr[age_order]
	donor_colour <- donor_colour[age_order]
	donor_sex <- donor_sex[age_order]
	donor_age <- donor_age[age_order]


	xlim = c(0, length(d_medians)+1)
	xes = 1:length(d_medians)
	yes = d_medians;
	CI95 = d_stderr*2;

	plot(xes, yes, xaxt="n", xlab="Donor (by age)", ylab="Scaled Expression", pch=16, col=donor_colour, xlim=xlim)
	arrows(xes, yes-CI95, xes, yes+CI95, length=0)
	if (type =="sex") {
		abline(h=mean(d_medians[donor_sex == "F"]), lty=2, col= colours_sex["F"])
		abline(h=mean(d_medians[donor_sex == "M"]), lty=2, col= colours_sex["M"])
	} 
	if (type == "age") {
		abline(v=max(which(donor_age < 40))+0.5, lty=2)
		abline(v=max(which(donor_age < 60))+0.5, lty=2)
	}
	if (type =="reject") {
		abline(h=mean(d_medians[donor_rej == "Y"]), lty=2, col= colours_reject["Y"])
		abline(h=mean(d_medians[donor_rej == "N"]), lty=2, col= colours_reject["N"])
	} 

	x_axis_ats <- seq(from=1, to=length(d_medians), by=axis.interval)
	axis(1, at=x_axis_ats, label=donor_age[x_axis_ats])

	return(yes)
}

