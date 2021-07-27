
require(Seurat)
source("/cluster/home/tandrews/scripts/LiverMap2.0/Colour_Scheme.R")
source("/cluster/home/tandrews/scripts/LiverMap2.0/My_R_Scripts.R")

dir = "/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Figures"

metaData <- readRDS("../Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno_MetaDataOnly.rds")

cluster_2_type <- table(metaData$Coarse_clusters, metaData$Coarse_Manual_Anno)
cluster_2_type <- colnames(cluster_2_type)[ apply(cluster_2_type, 1, function(x){which(x == max(x))}) ]
cluster_2_colour <- type_2_colour(cluster_2_type)

# Simpson Index
freq_by_sample <- table(metaData$Coarse_clusters, metaData$sample)
simpson_index <- function(freq_vec) {
	freq_vec <- freq_vec/sum(freq_vec)
	return(sum(freq_vec^2))
}

simpson_by_cluster <- apply(freq_by_sample, 1, simpson_index)
null <- simpson_index(table(metaData$sample));

#pdf("Figure2_Simpson_by_cluster.pdf", width=5, height=5)
png("Figure2_Simpson_by_cluster.png", width=5, height=5, unit="in", res=300)
par(mar=c(4,3,1,1))
barplot(simpson_by_cluster, col=cluster_2_colour, horiz=2, las=1, xlab="Simpson Index (Probabiliity)")
abline(v=null, lty=2)
dev.off()

# Frequency boxplot
# Proportion of sample for each cell-type
# Boxplots with dots

typefreq_by_sample <- table(metaData$sample, metaData$Coarse_Manual_Anno)
typefreq_by_sample <- typefreq_by_sample/rowSums(typefreq_by_sample)*100
bxp <- boxplot(as.vector(typefreq_by_sample)~rep(colnames(typefreq_by_sample), each=nrow(typefreq_by_sample)))
reorder <- rev(order(bxp$stats[3,]))

typefreq_by_sample <- typefreq_by_sample[,reorder]
overall_freq <- table(metaData$Coarse_Manual_Anno)[reorder]/nrow(metaData)

#pdf("Figure2_Freq_by_Sample_boxplot.pdf", width=10, height=6)
png("Figure2_Freq_by_Sample_boxplot.png", width=10, height=6.5, units="in", res=300)
par(mar=c(10,4,1,1))
bxp <- boxplot(as.vector(typefreq_by_sample)~factor(rep(colnames(typefreq_by_sample), each=nrow(typefreq_by_sample)), levels=colnames(typefreq_by_sample)), 
	names=paste(colnames(typefreq_by_sample), " (", round(overall_freq*100, digits=1),"%) ", sep=""), 
	col=type_2_colour(colnames(typefreq_by_sample)), las=2, 
	xlab="", ylab="Frequency (%)", outline=F, ylim=c(0,80))
points(jitter(rep(1:15, each=34)), as.vector(typefreq_by_sample), pch=16)
points(1:15, overall_freq*100, pch=23, bg="goldenrod1", cex=1.5)
legend("topright", "Full Map", bty="n", pch=23, cex=1.5, pt.bg="goldenrod1")
dev.off()

# Test Significance of Change of Proportion! - exclude Eryth for this analysis
typeN_sample <- table(metaData$sample, metaData$Coarse_Manual_Anno)
typeN_sample <- typeN_sample[,reorder]
typeN_sample <- typeN_sample[,-ncol(typeN_sample)]
typefreq_sample_prop <- typefreq_by_sample[,-ncol(typefreq_by_sample)]/100
null <- apply(typefreq_sample_prop,2,median)

# For each sample:
#	Test each cell-type X - is the ratio of X to each other cell-type different from the null?
#	If > 50% of pairwise comparisons are significant (BON multiple testing) in the same direction
#	Then cell-type X is significantly different not due to the difference in another cell-type.
is.diff <- typefreq_sample_prop*0
for (sample_i in 1:nrow(typefreq_sample_prop)) { 
	for (type_i in 1:ncol(typefreq_sample_prop)) {
		Zeds <- c();
		for (type_j in 1:ncol(typefreq_sample_prop)) {
			if (type_j == type_i) {
				Zeds <- c(Zeds, NA);
				next;
			}
			this_prop <- typefreq_sample_prop[sample_i, type_i]/(typefreq_sample_prop[sample_i, type_i]+typefreq_sample_prop[sample_i, type_j])
			null_prop <- null[type_i]/(null[type_i]+null[type_j])
			this_N <- typeN_sample[sample_i, type_i] +  typeN_sample[sample_i, type_j]
			Z <- (this_prop-null_prop)/sqrt((null_prop*(1-null_prop))/this_N)
			Zeds <- c(Zeds, Z)
		}
		threshold <- qnorm(1-0.05/length(Zeds))
		if (sum(Zeds > threshold, na.rm=T) > 0.75*length(Zeds)){ 
			is.diff[sample_i, type_i] <- 1;
		}
		if (sum(Zeds < -1*threshold, na.rm=T) > 0.75*length(Zeds)){ 
			is.diff[sample_i, type_i] <- -1;
		}
	}
}


is.diff <- cbind(is.diff, rep(0, nrow(is.diff)));
#pdf("Figure2_Freq_by_Sample_boxplot.pdf", width=10, height=6)
pt_col <- rep(rgb(0.4,0.4,0.4,0.6), length(as.vector(typefreq_by_sample)));
pt_col[as.vector(is.diff) != 0] <- rgb(0,0,0);

png("Figure2_Freq_by_Sample_boxplot_transparent.png", width=10, height=6.5, units="in", res=300)
par(mar=c(10,4,1,1))
bxp <- boxplot(as.vector(typefreq_by_sample)~factor(rep(colnames(typefreq_by_sample), each=nrow(typefreq_by_sample)), levels=colnames(typefreq_by_sample)), 
	names=paste(colnames(typefreq_by_sample), " (", round(overall_freq*100, digits=1),"%) ", sep=""), 
	col=type_2_colour(colnames(typefreq_by_sample)), las=2, 
	xlab="", ylab="Frequency (%)", outline=F, ylim=c(0,80))
points(jitter(rep(1:15, each=34)), as.vector(typefreq_by_sample), pch=16, col=pt_col)
points(1:15, overall_freq*100, pch=23, bg="goldenrod1", cex=1.5)
legend("topright", c("Full Map", "Differ (q < 0.05)"), bty="n", pch=c(23, 16), cex=1.5, pt.bg=c("goldenrod1", "black"))
dev.off()


# Sample-sample similarity by cluster
# Repeated HVGs only
# SF normalize the pseudobulks
# pearson correlations.

mergedobj <- readRDS("../Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno.rds")

hvgs <- mergedobj@misc$repeated_hvgs

obj_5pr <- mergedobj[rownames(mergedobj) %in% hvgs, mergedobj@meta.data$assay_type == "5pr"]
obj_3pr <- mergedobj[rownames(mergedobj) %in% hvgs, mergedobj@meta.data$assay_type == "3pr"]

dim(mergedobj)
cluster_5pr_pseudo <- get_pseudobulk(obj_5pr@assays$RNA@counts, factor(obj_5pr@meta.data$Coarse_clusters), factor(obj_5pr@meta.data$donor))
cluster_3pr_pseudo <- get_pseudobulk(obj_3pr@assays$RNA@counts, obj_3pr@meta.data$Coarse_clusters, obj_3pr@meta.data$donor)

cluster_5pr_pseudo <- t( t(cluster_5pr_pseudo)/colSums(cluster_5pr_pseudo) * median(colSums(cluster_5pr_pseudo), na.rm=T) )
cluster_3pr_pseudo <- t( t(cluster_3pr_pseudo)/colSums(cluster_3pr_pseudo) * median(colSums(cluster_3pr_pseudo), na.rm=T) )

id <- factor(unlist(lapply(strsplit(colnames(cluster_3pr_pseudo), "_"), function(x){x[[1]]})))
id2 <- factor(unlist(lapply(strsplit(colnames(cluster_5pr_pseudo), "_"), function(x){x[[1]]})))

d_5pr <- as.dist(1-cor(cluster_5pr_pseudo, method="pearson"))
d_3pr <- as.dist(1-cor(cluster_3pr_pseudo, method="pearson"))

scores <- c();
stderr <- c();

cross_scores <- c();

for (i in unique(id)) {
	score_3pr <- mean(as.matrix(1-d_3pr)[id==i, id==i], na.rm=T);
	crossscore_3pr <- mean(as.matrix(1-d_3pr)[id==i, id!=i], na.rm=T);
	score_5pr <- mean(as.matrix(1-d_5pr)[id2==i, id2==i], na.rm=T);
	crossscore_5pr <- mean(as.matrix(1-d_5pr)[id2==i, id2!=i], na.rm=T);
	sd_3pr <- sd(as.matrix(1-d_3pr)[id==i, id==i], na.rm=T)
	sd_5pr <- sd(as.matrix(1-d_5pr)[id2==i, id2==i], na.rm=T)
	n_cors_3pr <- sum(id==i)^2
	n_cors_5pr <- sum(id2==i)^2
	if (is.na(score_5pr)) {
		scores <- c(scores, score_3pr);
		stderr <- c(stderr, sd_3pr/sqrt(n_cors_3pr))
		cross_scores <- c(cross_scores, crossscore_3pr);
	} else {
		scores <- c(scores, (score_3pr+score_5pr)/2)
		stderr <- c(stderr, sqrt(sd_3pr^2+sd_5pr^2)/sqrt(n_cors_3pr+n_cors_5pr))
		cross_scores <- c(cross_scores, (crossscore_3pr+crossscore_5pr)/2)
	}
}

#pdf("Figure2_Similarity_by_cluster.pdf", width=5, height=5)
png("Figure2_Similarity_by_cluster.png", width=5, height=5, units="in", res=300)
par(mar=c(4,3,1,1))
bar_loc <- barplot(scores, name=unique(id), xlab="Cross Donor Correlation (average)", xlim=c(0,1), las=1, horiz=T, col=cluster_2_colour)
arrows(scores, bar_loc, scores+2*stderr,  bar_loc, angle=90, len=0)
legend("topright", lty=1, c("95% CI"), bty="n")
abline(v=mean(cross_scores), lty=2)
dev.off()


# Most variable genes by cell-type
obj_5pr <- mergedobj[ , mergedobj@meta.data$assay_type == "5pr"]
obj_3pr <- mergedobj[ , mergedobj@meta.data$assay_type == "3pr"]

all_cell_types <- sort(unique(obj_3pr@meta.data$Coarse_Manual_Anno))

all_scores=c()
all_means=c()
all_types=c()
all_genes=c();

all_cell_types <- all_cell_types[! all_cell_types %in% c("Eryth", "ErythHep")]

for (i in all_cell_types) {
	# 3pr
	this_3pr <- obj_3pr[,obj_3pr@meta.data$Coarse_Manual_Anno == i]
	means_3pr <- group_rowmeans(this_3pr@assays$RNA@data, this_3pr@meta.data$sample)
	btw_var <- group_rowvars(means_3pr, rep("A", ncol(means_3pr)))
	within_var <- rowMeans(group_rowvars(this_3pr@assays$RNA@data, this_3pr@meta.data$sample))
	scores_3pr <- log2(btw_var[,1]/within_var)
	scores_3pr[within_var == 0] <- 0;
	if (sum(obj_5pr@meta.data$Coarse_Manual_Anno == i) > 10) {
		# 5pr
		this_5pr <- obj_5pr[,obj_5pr@meta.data$Coarse_Manual_Anno == i]
		means_5pr <- group_rowmeans(this_5pr@assays$RNA@data, this_5pr@meta.data$sample)
		btw_var <- group_rowvars(means_5pr, rep("A", ncol(means_5pr)))
		within_var <- rowMeans(group_rowvars(this_5pr@assays$RNA@data, this_5pr@meta.data$sample))
		scores_5pr <- log2(btw_var[,1]/within_var)
		scores_5pr[within_var == 0] <- 0;
		# Score = log2(btw/within)
		this_means <- (my_rowMeans(means_3pr)+my_rowMeans(means_5pr))/2
		this_scores <- (scores_3pr+scores_5pr)/2
	} else {
		this_means <- my_rowMeans(means_3pr)
		this_scores <- scores_3pr
	}
	all_means <- c(all_means, this_means)
	all_scores <- c(all_scores, this_scores)
	all_types <- c(all_types, rep(i, nrow(this_3pr)))
	all_genes <- c(all_genes, rownames(obj_3pr));
	print(i)
	print(head(sort(this_scores, decreasing=T), 10))

	png(paste(i, "cross_variable_genes.png", sep="_"), width=12, height=5, units="in", res=300)
	plot(this_means, this_scores, pch=16, col=type_2_colour(i), xlab="mean", ylab="btw/within variability", ylim=c(-15, 5))
	to.label <- quantile(this_scores, prob=1-10/length(this_scores))
	to.label <- this_scores > to.label
	text(x=this_means[to.label], this_scores[to.label], names(this_scores)[to.label], pos=3)
	dev.off()
}

summary(all_scores)

keep <- is.finite(all_scores)

to.label <- quantile(all_scores[keep], prob=1-30/sum(keep));
to.label <- keep & all_scores > to.label

tmp <- unique(all_types[keep])

#pdf("Figure2_HighlyVariable_scatter.pdf", width=12, height=5)
png("Figure2_HighlyVariable_scatter.png", width=12, height=5, units="in", res=300)
plot(all_means[keep], all_scores[keep], pch=16, col=type_2_colour(all_types[keep]), xlab="mean", ylab="btw/within variability", ylim=c(-15, 5))
text(x=all_means[to.label], all_scores[to.label], names(all_scores)[to.label], pos=rep(c(3,1,4), times=10))
legend("bottomright", tmp, horiz=F, fill=type_2_colour(tmp), bty="n", ncol=5)
dev.off()

#expr_means_3pr <- group_rowmeans(obj_3pr@assays$RNA@data, obj_3pr@meta.data$Coarse_Manual_Anno)
#expr_means_5pr <- group_rowmeans(obj_5pr@assays$RNA@data, obj_5pr@meta.data$Coarse_Manual_Anno)

#var_group_3pr <- group_rowvars(obj_3pr@assays$RNA@data, obj_3pr@meta.data$Coarse_Manual_Anno)
#var_group_5pr <- group_rowvars(obj_5pr@assays$RNA@data, obj_5pr@meta.data$Coarse_Manual_Anno)

#overall_var_3pr <- group_rowvars(obj_3pr@assays$RNA@data, rep("all", ncol(obj_3pr)) )
#overall_var_5pr <- group_rowvars(obj_5pr@assays$RNA@data, rep("all", ncol(obj_5pr)) )

#l2fc_3pr <- log2(var_group_3pr/overall_var_3pr[,1])
#l2fc_5pr <- log2(var_group_5pr/overall_var_5pr[,1])

#l2fc_3pr[!is.finite(l2fc_3pr)] <- 0
#l2fc_5pr[!is.finite(l2fc_5pr)] <- 0

#l2fc_3pr <- l2fc_3pr[,match(colnames(l2fc_5pr), colnames(l2fc_3pr))]

# Consistent direction btw 3pr and 5pr
#disagree <- sign(l2fc_3pr) != sign(l2fc_5pr)
#l2fc_3pr[disagree] <- 0
#l2fc_5pr[disagree] <- 0

#l2fc_total <- l2fc_3pr + l2fc_5pr

#sort(rownames(l2fc_total)[which(apply(l2fc_total, 1, max) > 10)]) # Fibers (Stellate), Immunoglobulins (Bcells), CC (NK cell), CXC genes (cholangiocytes), mucus cholangiocytes


# xaxt = mean expression
# yaxt = log2( across sample var / within sample var)
# label top genes
# colour by cell-type



