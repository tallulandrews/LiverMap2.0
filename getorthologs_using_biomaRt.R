
# Human-mouse

require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attributes = c("ensembl_gene_id","mmusculus_homolog_ensembl_gene","mmusculus_homolog_perc_id_r1")

attributes=c(attributes,"mmusculus_homolog_orthology_type", "mmusculus_homolog_subtype", "mmusculus_homolog_perc_id")

orth.mouse = getBM( attributes,filters="with_mmusculus_homolog",values=TRUE, mart = human, bmHeader=FALSE)

write.table(orth.mouse, col.names=T, row.names=F, file="human-mouse_orthologs.txt")

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
h_ensg_symbol <- getBM(c("ensembl_gene_id", "hgnc_symbol"), mart=human, bmHeader=FALSE, values=TRUE)
write.table(h_ensg_symbol, col.names=T, row.names=F, file="Hsap_ensg2symbol.txt")

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
m_ensg_symbol <- getBM(c("ensembl_gene_id", "hgnc_symbol"), mart=mouse, bmHeader=FALSE, values=TRUE)
write.table(m_ensg_symbol, col.names=T, row.names=F, file="Mmus_ensg2symbol.txt")
