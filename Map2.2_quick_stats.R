obj <- readRDS("C46_reseq_Anno_SeurObj.rds")

dim(obj)

tot_per_cell <- Matrix::colSums(obj@assays$RNA@counts)
median(tot_per_cell)
sum(tot_per_cell)


gene_per_cell <- Matrix::colSums(obj@assays$RNA@counts>0)
median(gene_per_cell)
median(Matrix::colSums(obj@assays$RNA@counts))

mean(Matrix::colSums(obj@assays$RNA@counts[grep("^MT-",rownames(obj)),])/Matrix::colSums(obj@assays$RNA@counts))*100

mean(obj@meta.data$percent.mt)

table(obj@meta.data$seurat_clusters)

table(obj@meta.data$consistent_labs)
table(obj@meta.data$marker_labs)

