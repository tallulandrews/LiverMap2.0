require(Seurat)

C70_obj <- readRDS("/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C70_5pr_reseq_EmptyOnly.rds")

#cite <- as.sparse(read.csv(file = "", sep = ",", 
#    header = TRUE, row.names = 1))
cite <- readRDS("/cluster/projects/macparland/TA/DeepDive/antibody_x_cell_matrices/Ctmp_CITE_ED_lilst_all_Matrix.rds")

common_cells <- colnames(cite)[colnames(cite) %in% colnames(C70_obj)]

cite <- cite[,match(common_cells, colnames(cite))]
C70_obj <- C70_obj[,match(common_cells, colnames(C70_obj))]

cite <- cite[Matrix::rowSums(cite) > 5,]

table(C70_obj@meta.data$seurat_clusters, C70_obj@meta.data$marker_labs)
table(C70_obj@meta.data$seurat_clusters, C70_obj@meta.data$consistent_labs)


new.cluster.ids <- c("CentralHep1", "PortHep1", "PortHep2", "Hep1", "CentralHep2", "Hep2", "CentralHep3", "LSECs", "InfMac")
names(new.cluster.ids) <- as.character(levels(C70_obj@meta.data$seurat_clusters))
C70_obj <- SetIdent(C70_obj, value=C70_obj@meta.data$seurat_clusters)
C70_obj <- RenameIdents(C70_obj, new.cluster.ids)

png("Cluster_plot.png", width=7, height=6, units="in", res=300)
DimPlot(C70_obj, label=TRUE)
dev.off()

# Kruskal Wallace - one-way non-parametric ANOVA
group <- C70_obj@active.ident
ps <- apply(cite, 1, function(x){kruskal.test(x, group)$p.value})
qs <- p.adjust(ps, method="fdr")
sum(qs < 0.05)
sum(ps < 0.05)


signif <- rownames(cite)[ps < 0.05]

rownames(cite) <- paste("CITE", rownames(cite), sep=":")

source("~/scripts/LiverMap2.0/My_R_Scripts.R")
avg <- group_rowmeans(cite, group)

# add CITE data to object
C70_obj[["CITE"]] <- CreateAssayObject(counts = cite)
C70_obj <- NormalizeData(C70_obj, assay = "CITE", normalization.method = "CLR")
C70_obj <- ScaleData(C70_obj, assay = "CITE")

signif2 <- paste("CITE", gsub("_", "-", signif), sep=":")

ncol = 4
nrow = ceiling(length(signif2)/ncol);
pdf("All_the_plots.pdf", width=ncol*5, height=nrow*4.5)
FeaturePlot(C70_obj, features = signif2, min.cutoff = "q05", max.cutoff = "q95", ncol = ncol)
dev.off()

ncol = 3
nrow = ceiling(length(signif2)/ncol);
pdf("All_the_plots2.pdf", width=ncol*6, height=nrow*4.5)
VlnPlot(C70_obj, features = signif2, ncol = ncol)
dev.off()


## Fold-change in LSEC or Macs ##
signif3 <- avg[,"LSECs"]/rowMeans(avg)
signif3 <- names(tail(sort(signif3), 10))


signif4 <- avg[,"InfMac"]/rowMeans(avg)
signif4 <- names(tail(sort(signif4), 10))

signif5 <- c(signif3, signif4)

ncol = 4
nrow = ceiling(length(signif5)/ncol);
pdf("All_the_plots3.pdf", width=ncol*5, height=nrow*4.5)
FeaturePlot(C70_obj, features = signif5, min.cutoff = "q05", max.cutoff = "q95", ncol = ncol)
dev.off()

ncol = 3
nrow = ceiling(length(signif5)/ncol);
pdf("All_the_plots4.pdf", width=ncol*6, height=nrow*4.5)
VlnPlot(C70_obj, features = signif5, ncol = ncol)
dev.off()
