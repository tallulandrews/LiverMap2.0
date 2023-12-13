summary_stats_human <- function(seur_obj, assay="RNA") {
	Xist = "XIST";
	a <- read.delim("/cluster/projects/macparland/TA/ExternalData/NCBI_human_chrY_NOT_chrX_2023_05_23.txt")
	Ychr = as.character(a$Symbol);
	print_summary_stats(seur_obj, assay=assay, pattern.mt="^MT-", pattern.ribo="^RP[SL]", Ychr_genes=Ychr, Xist=Xist)
}
	
summary_stats_mouse <- function(seur_obj, assay="RNA") {
	Xist = "Xist";
	a <- read.delim("/cluster/projects/macparland/TA/ExternalData/NCBI_mouse_chrY_NOT_chrX_2023_05_23.txt")
	Ychr = as.character(a$Symbol);
	print_summary_stats(seur_obj, assay=assay, pattern.mt="^mt-", pattern.ribo="^Rp[sl]", Ychr_genes=Ychr, Xist=Xist)
}


print_summary_stats <- function(slice_obj, assay=c("RNA", "Spatial"), pattern.mt="^MT-", pattern.ribo="^RP[SL]", Ychr_genes, Xist="XIST") {
        slice_obj@meta.data$percent.mt <- PercentageFeatureSet(slice_obj, pattern = pattern.mt)
        slice_obj@meta.data$percent.ribo <- PercentageFeatureSet(slice_obj, pattern = pattern.ribo)

        ng = nrow(slice_obj)
        print(paste("n Genes =", ng))
        nc = ncol(slice_obj)
        print(paste("n Spots =", nc))
	slice_obj@meta.data[,paste("nCount", assay, sep="_")] <- Matrix::colSums(slice_obj@assays[[assay[1]]]@counts)
        tc = sum(slice_obj@meta.data[,paste("nCount", assay, sep="_")])
        print(paste("Total UMIs =", tc))
        umi_pC = median(slice_obj@meta.data[,paste("nCount", assay, sep="_")])
        print(paste("UMI/spot =", umi_pC))
        gene_pc = median(slice_obj@meta.data[,paste("nFeature", assay, sep="_")])
        print(paste("Gene/spot =", gene_pc))
        pct_mt = sum(slice_obj@meta.data[,paste("nCount", assay, sep="_")]*slice_obj@meta.data$percent.mt/100)/tc*100
        print(paste("% Mito =", pct_mt))
        pct_ribo = sum(slice_obj@meta.data[,paste("nCount", assay, sep="_")]*slice_obj@meta.data$percent.ribo/100)/tc*100
        print(paste("% Ribo =", pct_ribo))
	if (!( Xist %in% rownames(slice_obj))) {
		Xist <- 0
	} else {
		Xist <- mean(slice_obj@assays[[assay[1]]]@counts[Xist,])
	}
	Ychr <- mean(Matrix::colMeans(slice_obj@assays[[assay[1]]]@counts[rownames(slice_obj) %in% Ychr_genes,]))
	print(paste("Xist = ", Xist, "Ychr = ", Ychr))
}

