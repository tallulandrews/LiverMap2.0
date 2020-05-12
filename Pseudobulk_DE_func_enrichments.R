

# Get Gene descriptions from biomart
require("biomaRt")

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
listFilters(ensembl)

gene2desc <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id", "description"), filters="with_hgnc", mart=ensembl, values=TRUE)

require("gprofiler2")

files <- Sys.glob("*_pseudobulk_DE.rds")
proj_name = "20livers_v2"
set.seed(82710)

max_de <- 2000

for (i in 1:length(files)) {
	f <- files[i];
	bits <- unlist(strsplit(f, "_");)
	cluster_id <- bits[3];

	obj <- readRDS(f)

	pseudobulk_res <- obj[["edger"]]$table;
	pseudobulk_res$dir <- sign(pseudobulk_res$logFC);
	pseudobulk_res$des <- gene2desc[match(rownames(pseudobulk_res), gene2desc[,1]),3]

	de_quantile_thresh <- quantile(pseudobulk_res$PValue, max_de/nrow(pseudobulk_res))
	de_gene_list <- rownames(de_quantile_thresh)[
		pseudobulk_res$FDR < 0.05 & 
		pseudobulk$PValue < de_quantile_thresh & 
		pseudobulk$dir > 0]

res <- gprofiler2::gost(de_gene_list, ordered_query=T, correction_method="fdr", sources=c("GP:BP", "KEGG", "REAC"))$result

#pval, avg_logFC. pct.1, pct.2, p_val_adj, dir, desc

res <- res[order(res$p_value),]

lfc=(res$intersection_size/res$query_size)/(res$term_size/res$effective_domain_size)
rich_reformated <- data.frame(pval=res$p_value, avg_logFC=lfc, 
	pct1=res$intersection_size/res$query_size, 
	pct2=res$term_size/res$effective_domain_size,
	p_val_adj=res$p_value, dir=rep(1, nrow(res)), desc=res$term_name, 
	stringsAsFactors=FALSE)
rownames(rich_reformated) <- res$term_id
colnames(pseudobulk_res) <- colnames(rich_reformated)

tab <- rbind(pseudobulk_res, rich_reformated)
write.table(tab, paste(proj_name, c, "DE_for_anno.csv", sep="_"), sep=",", row.names=T, col.names=T)
}
