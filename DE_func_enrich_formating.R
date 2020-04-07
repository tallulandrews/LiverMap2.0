proj_name = "Harmony_Integrated"
set.seed(82710)
n_clusters = 23
cluster_names <- 1:n_clusters -1

require("gprofiler2")

require("biomaRt")

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
listFilters(ensembl)

gene2desc <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id", "description"), filters="with_hgnc", mart=ensembl, values=TRUE)


for (c in cluster_names) {
## for each cluster
if (!file.exists(paste("DEGeneFiles/", proj_name, "_", c, "_DE.rds", sep=""))) {next;}
obj <- readRDS(paste("DEGeneFiles/", proj_name, "_", c, "_DE.rds", sep=""))

wilcox_res <- obj[["seur_wilcox"]]
wilcox_res$dir <- sign(wilcox_res$avg_logFC)
wilcox_res$desc <- gene2desc[match(rownames(wilcox_res), gene2desc[,1]) ,3]

res <- gprofiler2::gost(rownames(wilcox_res)[wilcox_res$p_val_adj < 0.05 & wilcox_res$dir > 0],
 ordered_query=T, correction_method="fdr", sources=c("GP:BP", "KEGG", "REAC"))$result

#pval, avg_logFC. pct.1, pct.2, p_val_adj, dir, desc

res <- res[order(res$p_value),]

lfc=(res$intersection_size/res$query_size)/(res$term_size/res$effective_domain_size)
rich_reformated <- data.frame(pval=res$p_value, avg_logFC=lfc, 
	pct1=res$intersection_size/res$query_size, 
	pct2=res$term_size/res$effective_domain_size,
	p_val_adj=res$p_value, dir=rep(1, nrow(res)), desc=res$term_name, 
	stringsAsFactors=FALSE)
rownames(rich_reformated) <- res$term_id
colnames(wilcox_res) <- colnames(rich_reformated)

tab <- rbind(wilcox_res, rich_reformated)
write.table(tab, paste(proj_name, c, "DE_for_anno.csv", sep="_"), sep=",", row.names=T, col.names=T)
}