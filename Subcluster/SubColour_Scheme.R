colours_sex = c("pink", "dodgerblue"); names(colours_sex) <- c("F", "M")
colours_age = c("grey65", "grey25", "grey85"); names(colours_age) <- c("adult", "elderly", "young")
colours_reject = c("#7570b3", "grey55", "#1b9e77"); names(colours_reject) <- c("N", "?", "Y")

barplot_freq <- function(obj_in, cluster_col="Coarse_clusters", 
				colours=colours_sex, metadata="donor_sex") {
	obj_in <- obj_in[,obj_in@meta.data[,metadata] != "?"]

	if (metadata == "donor_age_group") {
		obj_in@meta.data[,metadata] <- factor(obj_in@meta.data[,metadata], levels=c("young", "adult", "elderly"))
	} else {
		obj_in@meta.data[,metadata] <- factor(obj_in@meta.data[,metadata])
	}
	tab_sample <- table(obj_in@meta.data[,cluster_col], obj_in@meta.data[,"sample"])
	keep <- colSums(tab_sample) > 50
	sample_meta <- table(obj_in@meta.data[,"sample"], obj_in@meta.data[,metadata])
	sample_lab <- apply(sample_meta, 1, function(x){colnames(sample_meta)[which(x == max(x))]})

	tab_sample <- t(tab_sample[,keep]);
	freq_sample <- t( t(tab_sample)/colSums(tab_sample) )
	n_sample <- rowSums(tab_sample)
	sample_stderr <- sqrt( freq_sample*(1-freq_sample)/n_sample )
	sample_lab <- sample_lab[keep]

	tab_class <- table(obj_in@meta.data[,cluster_col], obj_in@meta.data[,metadata])

	par(mfrow=c(2,1))
	par(mar=c(6,4,1,1))
	# plot 1
	barplot(t(tab_class), col=colours[colnames(tab_class)], ylab="N cells", las=2)
	legend("topright", colnames(tab_class), fill=colours[colnames(tab_class)], bty="n")

	# plot 2
	freq_class <- tab_class/rowSums(tab_class)
	horiz_lines <- colSums(tab_class)
	horiz_lines <- horiz_lines/sum(horiz_lines)
	bar_loc <- barplot(t(freq_class), col=colours[colnames(freq_class)], ylab="Proportion", las=2)
	line_loc = 0;
	lab_order <- colnames(freq_class)
	for(i in 1:(length(horiz_lines)-1)) {
		#for (cluster in 1:ncol(freq_sample)) {
		#	ys <- freq_sample[sample_lab==lab_order[i],cluster]
		#	points(jitter(rep(bar_loc[cluster], length(ys))), ys)
		#}
		line_loc=line_loc+horiz_lines[i];
		abline(h=line_loc, lty=2, lwd=1)
	}
}