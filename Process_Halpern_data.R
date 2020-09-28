
####### Halpern Spatial Stuff
Halpern_data <- read.delim("/cluster/projects/macparland/TA/ExternalData/HalpernSpatialLocHep/NIHMS70855-supplement-Supplementary_Table_3.csv", sep=",", header=T)

# Remap gene IDs
gene_ids <- Halpern_data[,1]
Hepatocyte_spatial_profiles <- Halpern_data[,2:10]
exclude <- Halpern_data$q.values > 0.05 | is.na(Halpern_data$q.values);
Hepatocyte_spatial_profiles <- Hepatocyte_spatial_profiles[!exclude,];
gene_ids <- gene_ids[!exclude];
Hepatocyte_spatial_profiles <- t(apply(Hepatocyte_spatial_profiles, 1, scale));
colnames(Hepatocyte_spatial_profiles) <- paste("Layer", 1:9, sep="")


id_list <- strsplit( as.character(gene_ids), ";")

n_ids <- unlist(lapply(id_list, length))
max_n_ids <- max(unlist(lapply(id_list, length)))

source("~/R-Scripts/Ensembl_Stuff.R")


map_id <- function(ids) {
        h_ids <- General_Map(ids, in.org="Mmus", out.org="Hsap", in.name="symbol", out.name="symbol")
        h_ids <- h_ids[h_ids != ""]
        if (length(h_ids) == 0) {return("")}
        if (length(h_ids) == 1) {return(h_ids)}
        tab <- table(h_ids)
        h_ids <- names(tab)[which(tab ==max(tab))]
        return(h_ids[1]);
}

h_genes_simple <- General_Map(unlist(id_list[n_ids==1]), in.org="Mmus", out.org="Hsap", in.name="symbol", out.name="symbol")
h_genes_complex <- sapply(id_list[n_ids > 1], map_id)

h_genes <- as.character(gene_ids);
h_genes[n_ids==1] <- h_genes_simple;
h_genes[n_ids>1] <- h_genes_complex;

# deal with duplicated genes
q.values <- Halpern_data$q.values[!exclude]
dups <- unique(h_genes[duplicated(h_genes)])
dups <- dups[dups != ""]

for (g in dups) {
        qs <- q.values[h_genes==g];
        best <- which(qs == min(qs));
        best <- best[1]
        keep <- which(h_genes==g)[best]
        h_genes[h_genes == g] <- ""
        h_genes[keep] <- g
}

Hepatocyte_spatial_profiles <- Hepatocyte_spatial_profiles[h_genes != "",]
h_genes <- h_genes[h_genes != ""]
rownames(Hepatocyte_spatial_profiles) <- h_genes
saveRDS(Hepatocyte_spatial_profiles, file="/cluster/projects/macparland/TA/ExternalData/HalpernSpatialLocHep/Human_ortho_sig_profiles.rds");
####### Halpern Done
